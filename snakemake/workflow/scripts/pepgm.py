import argparse

from peptonizer_rust import execute_pepgm_py


def str_to_bool(value: str) -> bool:
    """
    Convert a command-line string to a bool. Needed because argparse's built-in
    `type=bool` calls `bool(str)`, which is True for any non-empty string
    (so "--regularized False" would incorrectly parse as True).
    """
    value = value.strip().lower()
    if value in ("true", "1", "yes"):
        return True
    if value in ("false", "0", "no"):
        return False
    raise argparse.ArgumentTypeError(f"Invalid boolean value: {value}")


parser = argparse.ArgumentParser(
    description="Run the PepGM algorithm from command line"
)

parser.add_argument(
    "--communities-graph-bytes-path",
    type=str,
    required=True,
    help="Path to where the binary file of the factor graph (using Louvain communities) is stored.",
)
parser.add_argument(
    "--max-iter",
    nargs="?",
    type=int,
    default=10000,
    help="Max. number of iterations the belief propagation algo will go through.",
)
parser.add_argument(
    "--tol",
    nargs="?",
    type=float,
    default=0.006,
    help="Residual error allowed for the BP algorithm.",
)
parser.add_argument(
    "--out",
    type=str,
    required=True,
    help="Path to the file you want to save your results as.",
)
parser.add_argument(
    "--alpha",
    type=float,
    required=True,
    help="Detection probability of a peptide for the noisy-OR model.",
)
parser.add_argument(
    "--beta", type=float, required=True, help="Probability of wrong detection."
)
parser.add_argument(
    "--prior", type=float, required=True, help="Prior assigned to all effects."
)
parser.add_argument(
    "--regularized",
    type=str_to_bool,
    default=False,
    help="If True, the regularized version of the noisy-OR model is used.",
)

args = parser.parse_args()

with open(args.communities_graph_bytes_path, 'rb') as in_file:
    json_content = execute_pepgm_py(
        in_file.read(),
        args.alpha,
        args.beta,
        args.regularized,
        args.prior,
        args.max_iter,
        args.tol
    )

    with open(args.out, 'w') as out_file:
        out_file.write(json_content)
