import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import json
import argparse

from peptonizer_rust import get_names_for_taxa_py

"""
Script that takes PepGM .json output, translates taxIDS to scientific names, and barplots the *number of results* highest
scoring effects.
"""

def plot_peptonizer_results(input_file: str, output_file: str, number_of_effects: int = 25):
    """
    Read the results of a Peptonizer run from a JSON-file (denoted by the input_file argument) and write bar charts
    representing these results to a PNG-file.
    """
    assert input_file.lower().endswith(".json"), "Input file should be a JSON file."
    assert output_file.lower().endswith(".png"), "Output file should be a PNG file."

    # Read JSON file
    with open(input_file, "r") as f:
        data = json.load(f)
    effect_scores: Dict[int, float] = {
        int(k): float(v)
        for k, v in data.items()
    }

    # Get top N effects by score
    top_effects = sorted(
        effect_scores.items(),
        key=lambda x: x[1],
        reverse=True
    )[:number_of_effects]
    effect_ids = [effect_id for effect_id, _ in top_effects]
    effect_scores = [score for _, score in top_effects]
    effect_names_dict = json.loads(get_names_for_taxa_py(effect_ids))
    effect_names_dict: Dict[int, str] = {
        int(k): str(v)
        for k, v in effect_names_dict.items()
    }
    effect_names = [effect_names_dict[effect_id] for effect_id in effect_ids]

    # make the barplot
    fig, ax = plt.subplots()
    fig.set_size_inches(30, 15)
    bars = ax.barh(
        range(len(effect_names)),
        effect_scores,
        color="#283593",
    )

    ax.set_yticks(range(len(effect_names)))
    ax.set_yticklabels(effect_names, fontsize=24, color="#283593", fontweight="bold")
    ax.tick_params(axis='y', which='major', pad=15)
    plt.xlim((0, 1))
    plt.xlabel("Probability score", fontsize=35, fontweight="bold")
    ax.xaxis.set_ticks(np.arange(0, 1.2, 0.2))
    ax.xaxis.set_ticklabels([0, 0.2, 0.4, 0.6, 0.8, 1.0], fontsize=35)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["bottom"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.bar_label(bars, fmt='{:,.3f}', fontsize=24, fontweight='bold', color='black', padding=20)

    fig.tight_layout()

    plt.savefig(output_file)
    plt.close()

matplotlib.use("Agg")

parser = argparse.ArgumentParser(description="Generate BarPlot of PepGM results.")

parser.add_argument(
    "--results-file", type=str, help="Path(s) to your PepGM results JSON."
)
parser.add_argument(
    "--number-of-results",
    type=int,
    default=12,
    help="How many effects you want to show up on the results plot.",
)
parser.add_argument(
    "--out",
    type=str,
    help="Path(s) to where the generated BarPlot PNG's should be stored.",
)

args = parser.parse_args()
plot_peptonizer_results(args.results_file, args.results_file.replace(".json", ".png"), args.number_of_results)
