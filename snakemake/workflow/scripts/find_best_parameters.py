import argparse
import os
import pandas as pd
import re
import shutil
from os import path

from peptonizer_rust import compute_goodness_py, clean_csv_py

parser = argparse.ArgumentParser()


parser.add_argument(
    "--effect-cluster-heads-file",
    type=str,
    required=True,
    help="Input: path to a CSV-file that contains effect cluster heads.",
)
parser.add_argument(
    "--results-folder",
    type=str,
    required=True,
    help="Path to a folder containing JSON-files with all the results from a prior PepGM analysis.",
)
parser.add_argument(
    "--best-params-file",
    type=str,
    required=True,
    help="Path to the output file where the best suited parameter set should be stored in."
)
parser.add_argument(
    "--best-params-json",
    type=str,
    required=True,
    help="Path to the output file where the results of the best Peptonizer run in JSON format should be stored."
)
parser.add_argument(
    "--best-params-png",
    type=str,
    required=True,
    help="Path to the output file where the results of the best Peptonizer run in PNG format should be stored."
)

args = parser.parse_args()

def find_json_files(folder_path):
    json_files = []

    # Walk through the directory and subdirectories
    for root, dirs, files in os.walk(folder_path):
        for file in files:
            if file.endswith('.json') and file.find("pepgm_results") >= 0:
                json_files.append(os.path.join(root, file))

    return json_files

def extract_parameters(filename):
    # Regular expression to find the patterns 'aX', 'bX', 'pX' where X is a float
    pattern = r'_a([0-9.]+)_b([0-9.]+)_p([0-9.]+)\.'

    match = re.search(pattern, filename)

    if match:
        a = float(match.group(1))
        b = float(match.group(2))
        p = float(match.group(3))
        return a, b, p
    else:
        raise ValueError("The filename does not contain valid 'a', 'b', and 'p' parameters.")

# Get all effect cluster heads required to compute the goodness metric for each results file
effect_cluster_heads_csv = ""
with open(args.effect_cluster_heads_file, 'r') as effect_cluster_heads_file:
    effect_cluster_heads_csv = effect_cluster_heads_file.read()

# Store all result dataframes and the corresponding parameter sets in this list that will be used to finally find the
# best parameter set.
best_param_set = (0, 0, 0)
best_goodness = 0.0
for result_file in find_json_files(args.results_folder):
    alpha, beta, prior = extract_parameters(result_file)
    with open(result_file, "r") as f:
        peptonizer_result = f.read()
        goodness = compute_goodness_py(effect_cluster_heads_csv, peptonizer_result)
        if goodness > best_goodness:
            best_goodness = goodness
            best_param_set = (alpha, beta, prior)

# Write out the best parameters to a CSV file for future reference
(alpha, beta, prior) = best_param_set
with open(args.best_params_file, "w") as f:
    f.write("alpha,beta,prior\n")
    f.write(f"{alpha},{beta},{prior}\n")

# Clean the CSV for the best parameters and write it to the final output directory
best_json_path = path.join(args.results_folder, f"prior{prior}", f"pepgm_results_a{alpha}_b{beta}_p{prior}.json")
with open(best_json_path, "r") as in_file:
    clean_effects_json = clean_csv_py(in_file.read())

    with open(args.best_params_json, "w") as out_file:
        out_file.write(clean_effects_json)

# Copy the plots with the best parameters to the final output directory
shutil.copy(
    path.join(args.results_folder, f"prior{prior}", f"pepgm_results_a{alpha}_b{beta}_p{prior}.png"),
    args.best_params_png
)
