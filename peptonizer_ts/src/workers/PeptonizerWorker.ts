/**
 * This worker contains all instructions to run the different steps that are required by the Peptonizer. All of these
 * functions are implemented in the same Worker instance and are loaded at the same time into memory instead of into
 * smaller separate workers.
 *
 * @author Pieter Verschaffelt
 */

import {
    ClusterEffectsTaskData,
    ClusterEffectsTaskDataResult, ComputeGoodnessDataResult, ComputeGoodnessTaskData,
    ExecutePepgmTaskData,
    ExecutePepgmTaskDataResult,
    FetchUnipeptEffectTaskData,
    FetchUnipeptEffectTaskResult,
    GenerateGraphTaskData,
    GenerateGraphTaskDataResult,
    InputEventData,
    OutputEventData,
    PerformEffectsWeighingTaskData,
    PerformEffectsWeighingTaskResult,
    ResultType,
    WorkerTask
} from "./PeptonizerWorkerTypes.ts";
import init, { 
    perform_effects_weighing_wasm, 
    execute_pepgm_wasm, 
    fetch_unipept_taxa_wasm, 
    generate_pepgm_graph_wasm, 
    cluster_effects_wasm,
    compute_goodness_wasm
} from "../../pkg/peptonizer_rust.js";

interface DedicatedWorkerGlobalScope {
    postMessage: (message: OutputEventData) => void;
    submitPepgmProgress: (progressType: "graph" | "max_residual" | "iteration", currentValue: number, maxValue: number, workerId: number) => void;
}

declare const self: DedicatedWorkerGlobalScope & typeof globalThis;

async function fetchUnipeptEffectInformation(data: FetchUnipeptEffectTaskData): Promise<FetchUnipeptEffectTaskResult> {
    console.time("Execution time fetching Unipept information");
    
    let score_keys = [...data.peptidesScores.keys()];
    let peptidesScores = JSON.stringify(score_keys);
    let effectQuery = JSON.stringify(data.effectQuery);

    const unipeptJson = await fetch_unipept_taxa_wasm(peptidesScores, data.rank, effectQuery);
    
    console.timeEnd("Execution time fetching Unipept information");

    return { unipeptJson };
}

async function performEffectsWeighing(data: PerformEffectsWeighingTaskData): Promise<PerformEffectsWeighingTaskResult> {
    console.time("Execution time effects weiging");
    
    let peptidesEffects = JSON.stringify(Object.fromEntries(data.peptidesEffects));
    let peptidesScores = JSON.stringify(Object.fromEntries(data.peptidesScores));
    let peptidesCounts = JSON.stringify(Object.fromEntries(data.peptidesCounts));

    const [sequenceScoresCsv, effectsWeightsCsv] = await perform_effects_weighing_wasm(
        peptidesEffects,
        peptidesScores,
        peptidesCounts,
        data.effectsInGraph,
        data.rank
    );

    console.timeEnd("Execution time effects weiging");
    return {
        sequenceScoresCsv,
        effectsWeightsCsv
    };

    
}

async function generateGraph(data: GenerateGraphTaskData): Promise<GenerateGraphTaskDataResult> {
    console.time("Execution time generating graph");

    const factor_graph_bytes = generate_pepgm_graph_wasm(data.sequenceScoresCsv);
    
    console.timeEnd("Execution time generating graph");

    return {
        factor_graph_bytes
    };
}


async function executePepgm(data: ExecutePepgmTaskData): Promise<ExecutePepgmTaskDataResult> {
    console.time("Execution time Nori");

    const effectScoresJson = execute_pepgm_wasm(data.factor_graph_bytes, data.alpha, data.beta, true, data.prior);

    console.timeEnd("Execution time Nori");
    return {
        effectScoresJson
    };
}

async function clusterEffects(data: ClusterEffectsTaskData): Promise<ClusterEffectsTaskDataResult> {
    console.time("Execution time clustering effects");

    const effectClusterHeadsCsv = cluster_effects_wasm(data.sequenceScoresCsv, data.effectsWeightsCsv, data.similarityThreshold)

    console.timeEnd("Execution time clustering effects");
    return {
        effectClusterHeadsCsv
    };
}

async function computeGoodness(data: ComputeGoodnessTaskData): Promise<ComputeGoodnessDataResult> {
    console.time("Execution time computing goodness");
    
    let peptonizerResults = JSON.stringify(Object.fromEntries(data.peptonizerResults));
    const goodness = compute_goodness_wasm(data.effectClusterHeadsCsv, peptonizerResults);

    console.timeEnd("Execution time computing goodness");
    return {
        goodness
    }
}

self.submitPepgmProgress = function(
    progressType: "graph" | "max_residual" | "iteration",
    currentValue: number,
    maxValue: number,
    workerId: number
) {
    const resultMessage: OutputEventData = {
        resultType: ResultType.PROGRESS,
        task: WorkerTask.EXECUTE_PEPGM,
        workerId: workerId,
        progressUpdate: {
            progressType,
            currentValue,
            maxValue
        }
    }

    self.postMessage(resultMessage);
}


self.onmessage = async (event: MessageEvent<InputEventData>): Promise<void> => {
    try {
        // Make sure loading is done
        await init();

        // Destructure the data from the event
        const eventData = event.data;

        let output: FetchUnipeptEffectTaskResult | PerformEffectsWeighingTaskResult | GenerateGraphTaskDataResult | ExecutePepgmTaskDataResult | ClusterEffectsTaskDataResult | ComputeGoodnessDataResult | undefined;

        if (eventData.task === WorkerTask.FETCH_UNIPEPT_TAXON) {
            output = await fetchUnipeptEffectInformation(eventData.input);
        } else if (eventData.task === WorkerTask.PERFORM_TAXA_WEIGHING) {
            output = await performEffectsWeighing(eventData.input);
        } else if (eventData.task === WorkerTask.GENERATE_GRAPH) {
            output = await generateGraph(eventData.input);
        } else if (eventData.task === WorkerTask.EXECUTE_PEPGM) {
            output = await executePepgm(eventData.input);
        } else if (eventData.task === WorkerTask.CLUSTER_TAXA) {
            output = await clusterEffects(eventData.input);
        } else if (eventData.task === WorkerTask.COMPUTE_GOODNESS) {
            output = await computeGoodness(eventData.input);
        } else {
            throw new Error("Unknown task type passed to worker!");
        }

        if (!output) {
            throw new Error("No valid output was generated by worker!");
        }

        self.postMessage({
            resultType: ResultType.SUCCESSFUL,
            workerId: eventData.workerId,
            task: eventData.task,
            output
        });
    } catch (error: any) {
        self.postMessage({
            resultType: ResultType.FAILED,
            workerId: event.data.workerId,
            error: error.toString()
        });
    }
};
