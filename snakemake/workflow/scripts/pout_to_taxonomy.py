import argparse
import json

from peptonizer_rust import fetch_unipept_taxa_py, parse_unique_peptides_py


parser = argparse.ArgumentParser()

parser.add_argument(
    "--taxonomy-query",
    type=str,
    required=True,
    help="Taxa that should be used to query in Unipept. If querying all taxa, put [1].",
)
parser.add_argument(
    "--input-file",
    type=str,
    required=True,
    help="Input: path to peptide and score '.tsv' file.",
)
parser.add_argument(
    "--unipept-response-file",
    type=str,
    required=True,
    help="Path to output file that contains all queried peptide counts (which should be used in the next step)."
)
parser.add_argument(
    "--taxon-rank",
    type=str,
    required=False,
    default="species",
    help="Taxonomic rank at which you want the Peptonizer2000 results to be resolved.",
)

args = parser.parse_args()

with open(args.input_file, 'rt', encoding='utf-8') as file:
    file_contents = file.read()

# Parse the input TSV file
peptides = parse_unique_peptides_py(file_contents)

unipept_response = fetch_unipept_taxa_py(
    peptides,
    args.taxon_rank,
    args.taxonomy_query
)

with open(args.unipept_response_file, "w") as f:
    f.write(unipept_response)
