import PeptonizerWorker from './PeptonizerWorker.ts?worker&inline';
import {
    ClusterTaxaTaskData,
    ComputeGoodnessTaskData,
    ExecutePepgmTaskData,
    FetchUnipeptTaxonTaskData,
    GenerateGraphTaskData,
    InputEventData,
    OutputEventData,
    PepgmProgressUpdate,
    PerformTaxaWeighingTaskData,
    ResultType,
    SpecificInputEventData,
    WorkerTask
} from "./PeptonizerWorkerTypes.ts";
import async, { QueueObject } from "async";
import { PeptonizerResult } from "../Peptonizer.ts";
import { PeptonizerProgressListener } from "../PeptonizerProgressListener.ts";

/**
 * A worker pool that can be used to generate factor graphs for PepGM that can, in turn, be send as input to the belief
 * propagation algorithm.
 */
class WorkerPool {
    private workers: [Worker, number][] = [];
    private allWorkers: Worker[] = [];
    private queue: QueueObject<{ queueInput: SpecificInputEventData, progressListener?: PeptonizerProgressListener }>;
    private isCancelled: boolean = false;

    constructor(workerCount: number = 1) {
        for (let i = 0; i < workerCount; i++) {
            const worker = new PeptonizerWorker();
            this.allWorkers.push(worker);
            this.workers.push([worker, i]);
        }

        this.queue = async.queue(async(
            queueData
        ) => {
            // Retrieve worker from the pool.
            const [worker, workerId] = this.workers.pop()!;

            if (queueData.queueInput.task === WorkerTask.EXECUTE_PEPGM && queueData.progressListener) {
                const parameterSet = queueData.queueInput.input;
                queueData.progressListener.taskStarted(
                    {
                        alpha: parameterSet.alpha,
                        beta: parameterSet.beta,
                        prior: parameterSet.prior
                    },
                    workerId
                );
            }

            // eslint-disable-next-line @typescript-eslint/no-explicit-any -- queue result type varies per task, narrowed by each public method's return type
            const result = await new Promise<any>((resolve, reject) => {
                worker.onmessage = this.handleWorkerMessages(resolve, reject, queueData.progressListener);

                const workerTask: InputEventData = {
                    ...queueData.queueInput,
                    workerId
                }

                worker.postMessage(workerTask);
            });

            // Add worker back to the pool
            this.workers.push([worker, workerId]);

            if (queueData.queueInput.task === WorkerTask.EXECUTE_PEPGM && queueData.progressListener) {
                const parameterSet = queueData.queueInput.input;
                queueData.progressListener.taskFinished(
                    {
                        alpha: parameterSet.alpha,
                        beta: parameterSet.beta,
                        prior: parameterSet.prior
                    },
                    workerId
                );
            }

            return result;
        }, workerCount);
    }

    public async fetchUnipeptTaxonInfo(peptidesScores: Map<string, number>, rank: string, taxonQuery: number[]): Promise<string> {
        const eventData: FetchUnipeptTaxonTaskData = {
            peptidesScores,
            rank,
            taxonQuery
        };

        return await this.queue.pushAsync({ queueInput: { task: WorkerTask.FETCH_UNIPEPT_TAXON, input: eventData }, progressListener: undefined });
    }

    /**
     * Generates a CSV-file representing a dataframe with all the taxa weights required for the Peptonizer. These
     * taxa weights will be used in a subsequent step of the Peptonizer to generate the factor graph.
     *
     * @param peptidesTaxa Mapping between peptides and the associated taxa. If a filtering by taxa (or another
     * criterium) is required, this needs to be done before passing this mapping to this function.
     * @param peptidesScores Mapping between peptide sequences that need to be considered by the peptonizer and a
     * scoring value assigned to each sequence by prior steps (e.g. search engines).
     * @param peptidesCounts Mapping between peptide sequences and their occurrences in the input file.
     * @param taxaInGraph How many taxa are being used in the graphical model?
     * @return A CSV-representation of a dataframe with taxon weights.
     */
    public async performTaxaWeighing(
        peptidesTaxa: Map<string, number[]>,
        peptidesScores: Map<string, number>,
        peptidesCounts: Map<string, number>,
        taxaInGraph: number,
    ): Promise<[string, string]> {
        if (this.isCancelled) {
            throw new Error("Workerpool is no longer active. Cancel has been called on this pool before.");
        }

        const eventData: PerformTaxaWeighingTaskData = {
            peptidesTaxa,
            peptidesScores,
            peptidesCounts,
            taxaInGraph
        };

        return await this.queue.pushAsync({ queueInput: { task: WorkerTask.PERFORM_TAXA_WEIGHING, input: eventData }, progressListener: undefined });
    }

    public async generateGraph(
        sequenceScoresCsv: string
    ): Promise<Uint8Array> {
        if (this.isCancelled) {
            throw new Error("Workerpool is no longer active. Cancel has been called on this pool before.");
        }

        const eventData: GenerateGraphTaskData = {
            sequenceScoresCsv
        };

        return await this.queue.pushAsync({ queueInput: { task: WorkerTask.GENERATE_GRAPH, input: eventData }, progressListener: undefined });
    }

    public async executePepgm(
        factor_graph_bytes: Uint8Array,
        alpha: number,
        beta: number,
        prior: number,
        progressListener?: PeptonizerProgressListener
    ): Promise<PeptonizerResult> {
        if (this.isCancelled) {
            throw new Error("Workerpool is no longer active. Cancel has been called on this pool before.");
        }

        const eventData: ExecutePepgmTaskData = {
            factor_graph_bytes,
            alpha,
            beta,
            prior
        };

        return await this.queue.pushAsync({ queueInput: { task: WorkerTask.EXECUTE_PEPGM, input: eventData }, progressListener });
    }

    public async clusterTaxa(
        sequenceScoresCsv: string,
        taxaWeightsCsv: string,
        similarityThreshold: number = 0.9
    ): Promise<string> {
        if (this.isCancelled) {
            throw new Error("Workerpool is no longer active. Cancel has been called on this pool before.");
        }

        const eventData: ClusterTaxaTaskData = {
            sequenceScoresCsv,
            taxaWeightsCsv,
            similarityThreshold
        }

        return await this.queue.pushAsync({ queueInput: { task: WorkerTask.CLUSTER_TAXA, input: eventData }, progressListener: undefined });
    }

    public async computeGoodness(
        clusteredTaxaWeightsCsv: string,
        peptonizerResults: Map<string, number>
    ): Promise<number> {
        if (this.isCancelled) {
            throw new Error("Workerpool is no longer active. Cancel has been called on this pool before.");
        }

        const eventData: ComputeGoodnessTaskData = {
            clusteredTaxaWeightsCsv,
            peptonizerResults
        };

        return await this.queue.pushAsync({ queueInput: { task: WorkerTask.COMPUTE_GOODNESS, input: eventData }, progressListener: undefined });
    }

    /**
     * Stop execution of all currently running task and clean up all data that was left behind.
     */
    public close(): void {
        this.isCancelled = true;
        while (this.allWorkers.length > 0) {
            const worker = this.allWorkers.pop();
            worker?.terminate();
        }
    }

    /**
     * This function takes care of the results that are returned by the worker and converts them into something
     * usable for the rest of the framework.
     *
     * @param resolve
     * @param reject
     * @param progressListener
     * @private
     */
    private handleWorkerMessages(
        // eslint-disable-next-line @typescript-eslint/no-explicit-any -- resolves with whichever task-specific output type the caller expects
        resolve: (x: any) => void,
        // eslint-disable-next-line @typescript-eslint/no-explicit-any -- standard Promise executor reject signature
        reject: (reason?: any) => void,
        progressListener?: PeptonizerProgressListener
    ): (event: MessageEvent<OutputEventData>) => void {
        return (event: MessageEvent<OutputEventData>) => {
            const eventData = event.data;

            if (eventData.resultType === ResultType.SUCCESSFUL) {
                if (eventData.task === WorkerTask.FETCH_UNIPEPT_TAXON) {
                    resolve(eventData.output.unipeptJson)
                } else if (eventData.task === WorkerTask.PERFORM_TAXA_WEIGHING) {
                    resolve([eventData.output.sequenceScoresCsv, eventData.output.taxaWeightsCsv]);
                } else if (eventData.task === WorkerTask.GENERATE_GRAPH) {
                    resolve(eventData.output.factor_graph_bytes);
                } else if (eventData.task === WorkerTask.EXECUTE_PEPGM) {
                    const peptonizerResult: PeptonizerResult = new Map();
                    for (const [key, value] of Object.entries(JSON.parse(eventData.output.taxonScoresJson))) {
                        peptonizerResult.set(key, value as number);
                    }
                    resolve(peptonizerResult);
                } else if (eventData.task === WorkerTask.CLUSTER_TAXA) {
                    resolve(eventData.output.clusteredTaxaWeightsCsv);
                } else if (eventData.task === WorkerTask.COMPUTE_GOODNESS) {
                    resolve(eventData.output.goodness);
                }
            } else if (eventData.resultType === ResultType.PROGRESS) {
                if (progressListener && eventData.task === WorkerTask.EXECUTE_PEPGM) {
                    this.notifyPepgmProgressListener(eventData.progressUpdate, eventData.workerId, progressListener);
                }
            } else if (eventData.resultType === ResultType.FAILED) {
                reject(eventData.error);
            }
        }
    }

    private notifyPepgmProgressListener(
        progressUpdate: PepgmProgressUpdate,
        workerId: number,
        progressListener: PeptonizerProgressListener
    ): void {
        const currentValue = progressUpdate.currentValue;
        const maxValue = progressUpdate.maxValue;

        if (progressUpdate.progressType === "graph") {
            progressListener.graphsUpdated(currentValue, maxValue, workerId);
        } else if (progressUpdate.progressType === "max_residual") {
            progressListener.maxResidualUpdated(currentValue, maxValue, workerId);
        } else if (progressUpdate.progressType === "iteration") {
            progressListener.iterationsUpdated(currentValue, maxValue, workerId);
        }
    }
}

export { WorkerPool };
