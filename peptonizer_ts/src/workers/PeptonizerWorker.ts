/**
 * This worker contains all instructions to run the different steps that are required by the Peptonizer. All of these
 * functions are implemented in the same Worker instance and are loaded at the same time into memory instead of into
 * smaller separate workers.
 *
 * @author Pieter Verschaffelt
 */

import {
    ClusterTaxaTaskData,
    ClusterTaxaTaskDataResult, ComputeGoodnessDataResult, ComputeGoodnessTaskData,
    ExecutePepgmTaskData,
    ExecutePepgmTaskDataResult,
    FetchUnipeptTaxonTaskData,
    FetchUnipeptTaxonTaskResult,
    GenerateGraphTaskData,
    GenerateGraphTaskDataResult,
    InputEventData,
    OutputEventData,
    PerformTaxaWeighingTaskData,
    PerformTaxaWeighingTaskResult,
    ResultType,
    WorkerTask
} from "./PeptonizerWorkerTypes.ts";
import init, { 
    perform_taxa_weighing_wasm, 
    execute_pepgm_wasm, 
    fetch_unipept_taxa_wasm, 
    generate_pepgm_graph_wasm, 
    cluster_taxa_wasm,
    compute_goodness_wasm
} from "../../pkg/peptonizer_rust.js";

interface DedicatedWorkerGlobalScope {
    postMessage: (message: OutputEventData) => void;
    submitPepgmProgress: (progressType: "graph" | "max_residual" | "iteration", currentValue: number, maxValue: number, workerId: number) => void;
}

declare const self: DedicatedWorkerGlobalScope & typeof globalThis;

async function fetchUnipeptTaxonInformation(data: FetchUnipeptTaxonTaskData): Promise<FetchUnipeptTaxonTaskResult> {
    console.time("Execution time fetching Unipept information");
    
    const score_keys = [...data.peptidesScores.keys()];
    const peptidesScores = JSON.stringify(score_keys);
    const taxonQuery = JSON.stringify(data.taxonQuery);

    const unipeptJson = await fetch_unipept_taxa_wasm(peptidesScores, data.rank, taxonQuery);
    
    console.timeEnd("Execution time fetching Unipept information");

    return { unipeptJson };
}

async function performTaxaWeighing(data: PerformTaxaWeighingTaskData): Promise<PerformTaxaWeighingTaskResult> {
    console.time("Execution time taxa weiging");
    
    const peptidesTaxa = JSON.stringify(Object.fromEntries(data.peptidesTaxa));
    const peptidesScores = JSON.stringify(Object.fromEntries(data.peptidesScores));
    const peptidesCounts = JSON.stringify(Object.fromEntries(data.peptidesCounts));

    const [sequenceScoresCsv, taxaWeightsCsv] = await perform_taxa_weighing_wasm(peptidesTaxa, peptidesScores, peptidesCounts, data.taxaInGraph);

    console.timeEnd("Execution time taxa weiging");
    return {
        sequenceScoresCsv,
        taxaWeightsCsv
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

    const taxonScoresJson = execute_pepgm_wasm(data.factor_graph_bytes, data.alpha, data.beta, true, data.prior);

    console.timeEnd("Execution time Nori");
    return {
        taxonScoresJson
    };
}

async function clusterTaxa(data: ClusterTaxaTaskData): Promise<ClusterTaxaTaskDataResult> {
    console.time("Execution time clustering taxa");

    const clusteredTaxaWeightsCsv = cluster_taxa_wasm(data.sequenceScoresCsv, data.taxaWeightsCsv, data.similarityThreshold)

    console.timeEnd("Execution time clustering taxa");
    return {
        clusteredTaxaWeightsCsv
    };
}

async function computeGoodness(data: ComputeGoodnessTaskData): Promise<ComputeGoodnessDataResult> {
    console.time("Execution time computing goodness");
    
    const peptonizerResults = JSON.stringify(Object.fromEntries(data.peptonizerResults));
    const goodness = compute_goodness_wasm(data.clusteredTaxaWeightsCsv, peptonizerResults);

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

        let output: FetchUnipeptTaxonTaskResult | PerformTaxaWeighingTaskResult | GenerateGraphTaskDataResult | ExecutePepgmTaskDataResult | ClusterTaxaTaskDataResult | ComputeGoodnessDataResult | undefined;

        if (eventData.task === WorkerTask.FETCH_UNIPEPT_TAXON) {
            output = await fetchUnipeptTaxonInformation(eventData.input);
        } else if (eventData.task === WorkerTask.PERFORM_TAXA_WEIGHING) {
            output = await performTaxaWeighing(eventData.input);
        } else if (eventData.task === WorkerTask.GENERATE_GRAPH) {
            output = await generateGraph(eventData.input);
        } else if (eventData.task === WorkerTask.EXECUTE_PEPGM) {
            output = await executePepgm(eventData.input);
        } else if (eventData.task === WorkerTask.CLUSTER_TAXA) {
            output = await clusterTaxa(eventData.input);
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
    } catch (error) {
        self.postMessage({
            resultType: ResultType.FAILED,
            workerId: event.data.workerId,
            error: String(error)
        });
    }
};
