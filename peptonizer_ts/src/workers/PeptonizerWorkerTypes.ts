// Define a specific type of inputs that are expected for each task that can be performed by this worker.
enum WorkerTask {
    FETCH_UNIPEPT_TAXON,
    PERFORM_TAXA_WEIGHING,
    GENERATE_GRAPH,
    EXECUTE_PEPGM,
    CLUSTER_TAXA,
    COMPUTE_GOODNESS
}

interface FetchUnipeptEffectTaskData {
    peptidesScores: Map<string, number>;
    rank: string;
    effectQuery: number[];
}

interface PerformEffectsWeighingTaskData {
    peptidesEffects: Map<string, number[]>;
    peptidesScores: Map<string, number>;
    peptidesCounts: Map<string, number>;
    rank?: string;
    effectsInGraph: number;
}

interface GenerateGraphTaskData {
    sequenceScoresCsv: string;
}

interface ExecutePepgmTaskData {
    factor_graph_bytes: Uint8Array,
    alpha: number,
    beta: number,
    prior: number,
}

interface ClusterEffectsTaskData {
    sequenceScoresCsv: string,
    effectsWeightsCsv: string,
    similarityThreshold: number
}

interface ComputeGoodnessTaskData {
    effectClusterHeadsCsv: string,
    peptonizerResults: Map<string, number>
}

type SpecificInputEventData =
    { task: WorkerTask.FETCH_UNIPEPT_TAXON, input: FetchUnipeptEffectTaskData} |
    { task: WorkerTask.PERFORM_TAXA_WEIGHING, input: PerformEffectsWeighingTaskData } |
    { task: WorkerTask.GENERATE_GRAPH, input: GenerateGraphTaskData } |
    { task: WorkerTask.EXECUTE_PEPGM, input: ExecutePepgmTaskData } |
    { task: WorkerTask.CLUSTER_TAXA, input: ClusterEffectsTaskData } |
    { task: WorkerTask.COMPUTE_GOODNESS, input: ComputeGoodnessTaskData };

type CommonInputEventData = { workerId: number };

type InputEventData = SpecificInputEventData & CommonInputEventData;

interface FetchUnipeptEffectTaskResult {
    unipeptJson: string,
}

interface PerformEffectsWeighingTaskResult {
    sequenceScoresCsv: string,
    effectsWeightsCsv: string
}

interface GenerateGraphTaskDataResult {
    factor_graph_bytes: Uint8Array
}

interface ExecutePepgmTaskDataResult {
    effectScoresJson: string
}

interface PepgmProgressUpdate {
    progressType: "graph" | "max_residual" | "iteration",
    currentValue: number,
    maxValue: number
}

interface ClusterEffectsTaskDataResult {
    effectClusterHeadsCsv: string
}

interface ComputeGoodnessDataResult {
    goodness: number
}

enum ResultType {
    SUCCESSFUL,
    PROGRESS,
    FAILED,
    CANCELLED
}

type SpecificOutputEventData = { resultType: ResultType.SUCCESSFUL } & (
    { task: WorkerTask.FETCH_UNIPEPT_TAXON, output: FetchUnipeptEffectTaskResult } |
    { task: WorkerTask.PERFORM_TAXA_WEIGHING, output: PerformEffectsWeighingTaskResult } |
    { task: WorkerTask.GENERATE_GRAPH, output: GenerateGraphTaskDataResult } |
    { task: WorkerTask.EXECUTE_PEPGM, output: ExecutePepgmTaskDataResult } |
    { task: WorkerTask.CLUSTER_TAXA, output: ClusterEffectsTaskDataResult } |
    { task: WorkerTask.COMPUTE_GOODNESS, output: ComputeGoodnessDataResult });

type CommonOutputEventData = { workerId: number };

type ErrorOutputEvent = { resultType: ResultType.FAILED, error: string };
type ProgressOutputEvent = { resultType: ResultType.PROGRESS, task: WorkerTask.EXECUTE_PEPGM, progressUpdate: PepgmProgressUpdate };

type OutputEventData = (SpecificOutputEventData | ErrorOutputEvent | ProgressOutputEvent) & CommonOutputEventData;

export {
    WorkerTask,
    ResultType
};

export type {
    FetchUnipeptEffectTaskData, 
    PerformEffectsWeighingTaskData,
    GenerateGraphTaskData,
    ExecutePepgmTaskData,
    ClusterEffectsTaskData,
    ComputeGoodnessTaskData,
    SpecificInputEventData,
    InputEventData,
    FetchUnipeptEffectTaskResult,
    PerformEffectsWeighingTaskResult,
    GenerateGraphTaskDataResult,
    ExecutePepgmTaskDataResult,
    ClusterEffectsTaskDataResult,
    ComputeGoodnessDataResult,
    OutputEventData,
    PepgmProgressUpdate,
    ProgressOutputEvent
};
