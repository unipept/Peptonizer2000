import argparse

from peptonizer_rust import cluster_effects_py

parser = argparse.ArgumentParser(description = 'cluster Effects based on peptidome similarity and weight attributed')

parser.add_argument(
    '--sequence-scores-dataframe-file',
    type = str,
    help = 'Path to the sequence scores dataframe file.'
)
parser.add_argument(
    '--effects-weights-dataframe-file',
    type = str,
    help = 'Path to file with weighted effects computed in the effects weighing step.'
)
parser.add_argument(
    '--similarity-threshold',
    type = float,
    help = 'Threshold for the peptidome similarity at which two effects should belong to the same cluster.'
)
parser.add_argument(
    '--out',
    type = str,
    help= 'Path to clustered effects output csv file.'
)

args = parser.parse_args()

with open(args.sequence_scores_dataframe_file, 'r') as sequence_scores_file, open(args.effects_weights_dataframe_file, 'r') as effects_weights_file:
    sequence_scores_csv = sequence_scores_file.read()
    effects_weights_csv = effects_weights_file.read()

effect_cluster_heads_csv = cluster_effects_py(
    sequence_scores_csv,
    effects_weights_csv,
    args.similarity_threshold
)

with open(args.out, 'w') as clustered_effects_file:
    clustered_effects_file.write(effect_cluster_heads_csv)
