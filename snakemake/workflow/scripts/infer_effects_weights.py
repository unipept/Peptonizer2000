import argparse
import gzip

from peptonizer_rust import perform_effects_weighing_py, parse_input_peptides_py

parser = argparse.ArgumentParser()

parser.add_argument(
    "--number-of-effects",
    type=int,
    required=True,
    help="Number of effects to include in the final Peptonizer2000 output.",
)
parser.add_argument(
    "--sequence-scores-dataframe-file",
    type=str,
    required=False,
    help="Output: path to a CSV-file that will contain all computed sequence scores.",
)
parser.add_argument(
    "--effects-weights-dataframe-file",
    type=str,
    required=False,
    help="Output: path to a CSV-file that will contain all computed effects weights.",
)
parser.add_argument(
    "--unipept-response-file",
    type=str,
    required=True,
)
parser.add_argument(
    "--effect-rank",
    type=str,
    required=False,
    default="species",
    help="Taxonomic rank at which you want the Peptonizer2000 results to be resolved.",
)
parser.add_argument(
    "--pout-file",
    type=str,
    required=True,
    help="Input: path to percolator (ms2rescore) '.pout' file.",
)

args = parser.parse_args()

with open(args.pout_file, 'rt', encoding='utf-8') as file:
    file_contents = file.read()

# Parse the input MS2Rescore file
pep_score, pep_psm_counts = parse_input_peptides_py(file_contents)


# Read the Unipept response file
unipept_responses = ""
with open(args.unipept_response_file, "r") as file:
    unipept_responses = file.read()

sequence_scores, effects_weights = perform_effects_weighing_py(
    unipept_responses,
    pep_score,
    pep_psm_counts,
    args.number_of_effects,
    args.effect_rank
)

print("Started dumping produced results to CSV-files...")
with open(args.sequence_scores_dataframe_file, 'w') as sequences_file, open(args.effects_weights_dataframe_file, 'w') as weights_file:
    sequences_file.write(sequence_scores)
    weights_file.write(effects_weights)
