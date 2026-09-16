import {PeptonizerParameterSet, PeptonizerProgressListener} from "./PeptonizerProgressListener.ts";
import {WorkerPool} from "./workers/WorkerPool.ts";

type PeptonizerResult = Map<string, number>;

class Peptonizer {
    private isCancelled: boolean = false;
    private workerPool!: WorkerPool;

    /**
     * Start the peptonizer. This function takes in a PSM-file that has been read in earlier (so no file paths here). The
     * PSMS including their intensities are then used as input to the Peptonizer-pipeline. This pipeline will finally
     * return a Map in which NCBI effect IDs are mapped onto their probabilities (as computed by the Peptonizer).
     *
     * @param peptidesEffects Mapping between peptides and the associated effects. If a filtering by effects (or another
     * criterium) is required, this needs to be done before passing this mapping to this function.
     * @param peptidesScores Mapping between the peptide sequences that should be present in the peptonizer and a scoring
     * metric (derived from a prior search engine step) for each peptide.
     * @param peptidesCounts Mapping between peptide sequences and the amount of times they occur in the input.
     * @param alphas An array of possible values for the alpha parameter. This parameter indicates the probability that an
     * observed effect also indicates the presence of a peptide.
     * @param betas An array of possible values for the beta parameter. This parameter indicates the probability of
     * detecting a peptide at random.
     * @param priors An array of possible values for the gamma (or prior) parameter. Gamma indicates the prior probability
     * of a effect being present.
     * @param rank At which NCBI effect rank should the Peptonizer perform the effect inference.
     * @param effectsInGraph How many effects are being used in the graphical model?
     * @param progressListener Is called everytime the progress of the belief propagation algorithm has been updated.
     * @param workers The amount of Web Workers that can be spawned and used simultaneously to run the Peptonizer.
     * @return Mapping between NCBI effect IDs (integer, > 0) and probabilities (float in [0, 1]). If execution of the
     * Peptonizer was cancelled before it was completed, it returns an undefined result set.
     */
    async peptonize(
        peptidesEffects: Map<string, number[]>,
        peptidesScores: Map<string, number>,
        peptidesCounts: Map<string, number>,
        alphas: number[],
        betas: number[],
        priors: number[],
        rank: string = "species",
        effectsInGraph: number = 100,
        progressListener?: PeptonizerProgressListener,
        workers: number = 2
    ): Promise<PeptonizerResult | undefined> {
        return this.runPipeline(
            peptidesEffects,
            peptidesScores,
            peptidesCounts,
            alphas,
            betas,
            priors,
            rank,
            effectsInGraph,
            progressListener,
            workers
        );
    }

    /**
     * Start the peptonizer in functional-analysis mode. Identical pipeline to `peptonize()`, except that
     * `peptidesFunctions` maps peptides onto functional annotation IDs instead of taxonomic effect IDs, and no
     * NCBI rank normalization is performed (functional annotations have no such rank).
     *
     * @param peptidesFunctions Mapping between peptides and the associated functional annotation IDs.
     * @param peptidesScores Mapping between the peptide sequences that should be present in the peptonizer and a scoring
     * metric (derived from a prior search engine step) for each peptide.
     * @param peptidesCounts Mapping between peptide sequences and the amount of times they occur in the input.
     * @param alphas An array of possible values for the alpha parameter. This parameter indicates the probability that an
     * observed effect also indicates the presence of a peptide.
     * @param betas An array of possible values for the beta parameter. This parameter indicates the probability of
     * detecting a peptide at random.
     * @param priors An array of possible values for the gamma (or prior) parameter. Gamma indicates the prior probability
     * of a effect being present.
     * @param functionsInGraph How many functions are being used in the graphical model?
     * @param progressListener Is called everytime the progress of the belief propagation algorithm has been updated.
     * @param workers The amount of Web Workers that can be spawned and used simultaneously to run the Peptonizer.
     * @return Mapping between functional annotation IDs and probabilities (float in [0, 1]). If execution of the
     * Peptonizer was cancelled before it was completed, it returns an undefined result set.
     */
    async functionalAnalysis(
        peptidesFunctions: Map<string, number[]>,
        peptidesScores: Map<string, number>,
        peptidesCounts: Map<string, number>,
        alphas: number[],
        betas: number[],
        priors: number[],
        functionsInGraph: number = 100,
        progressListener?: PeptonizerProgressListener,
        workers: number = 2
    ): Promise<PeptonizerResult | undefined> {
        return this.runPipeline(
            peptidesFunctions,
            peptidesScores,
            peptidesCounts,
            alphas,
            betas,
            priors,
            undefined,
            functionsInGraph,
            progressListener,
            workers
        );
    }

    /**
     * Shared pipeline behind both `peptonize()` and `functionalAnalysis()`: effect weighing, factor graph
     * generation, a belief-propagation grid search over every parameter set, clustering, and goodness-based
     * selection of the best result. `rank` is forwarded to the effect-weighing step for taxonomic analysis;
     * `functionalAnalysis()` passes `undefined` since functional annotations have no NCBI rank to normalize to.
     */
    private async runPipeline(
        peptidesEffects: Map<string, number[]>,
        peptidesScores: Map<string, number>,
        peptidesCounts: Map<string, number>,
        alphas: number[],
        betas: number[],
        priors: number[],
        rank: string | undefined,
        effectsInGraph: number,
        progressListener?: PeptonizerProgressListener,
        workers: number = 2
    ): Promise<PeptonizerResult | undefined> {
        try {
            this.workerPool = new WorkerPool(workers);

            // Compute the total amount of tasks and which tasks to perform.
            const parameterSets: PeptonizerParameterSet[] = [];

            for (const alpha of alphas) {
                for (const beta of betas) {
                    for (const prior of priors) {
                        parameterSets.push({
                            alpha,
                            beta,
                            prior
                        });
                    }
                }
            }

            // Notify any listeners that the Peptonizer did start running (and report which set of parameters will be tuned)
            progressListener?.peptonizerStarted(parameterSets.length, parameterSets);

            // const unipept_json = await this.workerPool.fetchUnipeptEffectInfo(peptidesScores, rank, effectQuery);

            const effectWeighingResult = await this.workerPool.performEffectsWeighing(
                peptidesEffects,
                peptidesScores,
                peptidesCounts,
                rank,
                effectsInGraph
            );

            if (this.isCancelled) {
                return;
            }
            const [sequenceScoresCsv, effectWeightsCsv] = effectWeighingResult;

            const factor_graph_bytes = await this.workerPool.generateGraph(sequenceScoresCsv);

            const pepgmPromises: Promise<PeptonizerResult>[] = [];

            if (this.isCancelled) {
                return;
            }

            console.time("Execution time pepgm total");
            for (const paramSet of parameterSets) {
                pepgmPromises.push(
                    this.workerPool.executePepgm(factor_graph_bytes, paramSet.alpha, paramSet.beta, paramSet.prior, progressListener)
                );
            }

            // Wait until all parameter sets have been tuned...
            const peptonizerResults = await Promise.all(pepgmPromises);
            console.timeEnd("Execution time pepgm total");

            // Now that we have all the results generated by the peptonizer, we need to figure out which one yields the best
            // results

            if (this.isCancelled) {
                return;
            }

            // First compute the clustered effects weights
            const effectClusterHeadsCsv = await this.workerPool.clusterEffects(sequenceScoresCsv, effectWeightsCsv);

            if (this.isCancelled) {
                return;
            }

            let bestGoodness = -1;
            let bestResult: PeptonizerResult | undefined;
            for (const result of peptonizerResults) {
                const goodness = await this.workerPool.computeGoodness(effectClusterHeadsCsv, result);

                if (goodness > bestGoodness) {
                    bestGoodness = goodness;
                    bestResult = result;
                }
            }

            if (this.isCancelled) {
                return undefined;
            }

            if (!bestResult) {
                throw new Error("No results found!");
            }

            this.workerPool.close();
            progressListener?.peptonizerFinished();

            return bestResult;
        } catch (error) {
            throw error;
        } finally {
            this.workerPool.close();
        }
    }

    public cancel(): void {
        this.isCancelled = true;
        this.workerPool.close();
    }
}

export { Peptonizer };
export type { PeptonizerResult }
