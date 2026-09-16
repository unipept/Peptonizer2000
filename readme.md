<div id="top"></div>


<!-- PROJECT SHIELDS -->
<!--
*** I'm using markdown "reference style" links for readability.
*** Reference links are enclosed in brackets [ ] instead of parentheses ( ).
*** See the bottom of this document for the declaration of the reference variables
*** for contributors-url, forks-url, etc. This is an optional, concise syntax you may use.
*** https://www.markdownguide.org/basic-syntax/#reference-style-links
-->
<!-- PROJECT LOGO -->
<br />
<div align="center">
  <a href=https://git.bam.de/tholstei/pepgm/>
    <img src="https://raw.githubusercontent.com/compomics/Peptonizer2000/refs/heads/master/peptonizer_logo.jpg" alt="Logo"  height="300">
  </a>

<h3 align="center">The Peptonizer 2000</h3>

  <p align="center">
    Integrating PepGM and Unipept for probability-based effect inference of metaproteomic samples
    <br />
  </p>
</div>


<!-- TABLE OF CONTENTS -->
<details>
  <summary>Table of Contents</summary>
  <ol>
    <li>
      <a href="#about-the-project">About The Project</a>
      <ul>
      </ul>
    </li>
    <li><a href="#input">Input</a></li>
    <li>
      <a href="#getting-started">Getting Started</a>
      <ul>
        <li><a href="#prerequisites">Prerequisites</a></li>
        <li><a href="#installation">Installation</a></li>
        <li><a href="#preparation">Preparation</a></li>
      </ul>
    </li>
    <li><a href="#usage">Usage</a></li>
    <li><a href="#running-an-analysis-natively-without-snakemake-or-the-website">Running an analysis natively, without Snakemake or the website</a></li>
    <li><a href="#roadmap">Roadmap</a></li>
    <li><a href="#contributing">Contributing</a></li>
    <li><a href="#license">License</a></li>
    <li><a href="#contact">Contact</a></li>
  </ol>
</details>



<!-- ABOUT THE PROJECT -->
## About The Project

Introducing the Peptonizer2000 - a tool that combines the capabilities of Unipept and PepGM to analyze
metaproteomic mass spectrometry-based samples. Originally designed for taxonomic inference of viral
mass spectrometry-based samples, we've extended PepGM's functionality to analyze metaproteomic samples by
retrieving taxonomic information from the Unipept database. The pipeline can also be used for other peptide-effect 
relations, such as functional analysis.

PepGM is a probabilistic graphical model developed by Tanja Holstein et al. that uses belief propagation to infer the taxonomic origin of peptides and taxa in viral samples.
You can learn more about PepGM at [GitHub](https://github.com/BAMeScience/PepGM) page.

Unipept, on the other hand, is a web-based metaproteomics analysis tool that provides taxonomic information for
identified peptides. To make it work seamlessly with PepGM, we've extended Unipept with new functionalities that
restrict the taxa queried and provide all potential taxonomic origins of the peptides queried. Check out more
information about Unipept [here](https://unipept.ugent.be/).

With the Peptonizer2000, you can look forward to a comprehensive and streamlined workflow that simplifies
the process of identifying peptides and their taxonomic origins in metaproteomic samples.

The Peptonizer2000 workflow is comprised of the following steps:

1. Query all identified peptides, provided by the user in a .tsv file, in the Unipept API,
   and restrict the taxonomic range queried based on any prior knowledge of the sample.
2. Assemble the peptide-effect associations provided by Unipept into a bipartite graph,
   where peptides and effects are represented by different nodes, and an edge is drawn between a peptide and a effect
   if they are related.
3. Transform the bipartite graph into a factor graph using convolution trees and conditional probability table
   factors (CPD).
4. Run the belief propagation algorithm multiple times with different sets of CPD parameters until convergence,
   to obtain posterior probabilities of candidate effects.
5. Use an empirically deduced metric to determine the ideal graph parameter set.
6. Output the top scoring effects as a results barchart. The results are also available as comma-separated files
   for further downstream analysis or visualizations.


<div align="center">
    <img src="https://raw.githubusercontent.com/compomics/Peptonizer2000/refs/heads/master/peptonizer_workflow.png" alt="workflow scheme" width="500">
</div>

<br>



<p align="right">(<a href="#top">back to top</a>)</p>

<!-- INPUT -->

## Input

* A .tsv file of your peptides output from any protoemic peptide search method. The first column should be the peptide, the second column it's score attributed by the search engine. An example is provided in test files. <br>
* A config file with your parameters for the peptonizer2000. A more detailed description of the configuration file can be found below. Additionally, an exemplary config file is provided in this repository.

<p align="right">(<a href="#top">back to top</a>)</p>

<!-- GETTING STARTED -->
## Getting Started

### Prerequisites

## Rust implementation
The core algorithm is implemented in Rust in the `peptonizer_rust` folder. Wheels are created for all major platforms so users can use the package on supported systems.

### Running as snakemake workflow
In order to run the Peptonizer2000 on your own system, you should install Conda, Mamba and all of its dependencies.
Follow the installation instructions step-by-step for an explanation of what you should do.

* Make sure that Conda and Mamba are installed. If these are not yet present on your system, you can follow the instructions on their [README](https://github.com/conda-forge/miniforge).
* Go to the "workflow" directory by executing `cd snakemake/workflow` from the terminal.
* Run `conda env create -f env.yml` (make sure to run this command from the workflow directory) in order to install all dependencies and create a new conda environment (which is named "peptonizer" by default).
* Run `mamba install -c conda-forge -c bioconda -n peptonizer snakemake` to install snakemake which is required to run the whole workflow from start-to-finish.
* Run `conda activate peptonizer` to switch the current Conda environment to the peptonizer environment you created earlier.
* Start the peptonizer with the command `snakemake --use-conda --cores 1`. If you have sufficient CPU and memory power available to your system, you can increase the amount of cores in order to speed up the workflow.

If you see the following error while installing dependencies:

```
ERROR: Could not find a version that satisfies the requirement peptonizer_rust (from versions: none)
ERROR: No matching distribution found for peptonizer_rust
```

then the workflow could not find a wheel for `peptonizer_rust` for your platform. To create one manually, expand the instructions below.

<details>
  <summary>Show manual wheel creation instructions</summary>

  - Change to the Rust package directory:

     ```bash
     cd peptonizer_rust
     ```

  - Download the rustup installer:

    ```bash
    curl -sSf -o rustup-init.sh https://sh.rustup.rs
    ```

  - Run the installer with defaults:

    ```bash
    sh rustup-init.sh -y
    ```

  - Make `cargo` available in this shell session:

    ```bash
    source "$HOME/.cargo/env"
    ```

  - Update the Rust toolchain to stable:

    ```bash
    rustup update
    ```

  - Install Linux build tools (Debian/Ubuntu example):

    ```bash
    sudo apt-get update
    sudo apt-get install -y build-essential
    ```

  - Install `maturin` into the active Conda environment (run after `conda activate peptonizer`):

    ```bash
    python -m pip install --upgrade pip setuptools wheel
    python -m pip install maturin
    ```

  - Build the wheel for Python 3.12 (adjust `-i` if using a different Python):

    ```bash
    maturin build --release --out dist -i python3.12
    ```

  The command writes one or more `.whl` files to `peptonizer_rust/dist`.

</details>

### Configuration file

The Peptonizer2000 relies on a configuration file in `yaml` format to set up the workflow.
An example configuration file is provided in `config/config.yaml`. <br>
Do not change the config file location.

<details> 
   <details > <summary> Directory parameters </summary>
   <ul>
      <li>data_dir: relative path to output files </li>
      <li>input_file: relative path to input .tsv </li> 
      <li>log_dir: relative path to log directory</li>
   </ul>
   </details>

   <details > <summary> Analysis specific parameter </summary>
   <ul>
      <li>effects_in_graph: # of inferred effects that appear in the barplot that is created of the results csv</li>
      <li>effects_in_plot: number of effects reported in bar plot</li>
      <li>alpha: grid search increments for alpha (list) </li>
      <li>beta: grid search increments for beta (list) </li>
      <li>prior: grid search increments for prior (list) </li>
      <li>regularized: boolean. If True, the probability for the number of parents effects of a peptide is regularized to be inversely proportional to the number of parents </li>
   </ul>
   </details>
   <details > <summary> UniPept query parameters </summary>
   <ul>
       <li>taxon_rank: NCBI rank at which taxonomic results will be reported </li>
       <li>taxon_query: taxa comprised in the Unipept query. If querying all of Unipept, use 1 (list)</li>
   </ul> 
   </details>
</details>

### Output files

All Peptonizer2000 output files are saved into the results folder and include the following: <br>

Main results: <br>

- peptonizer_results.csv: table with values ID, score, type (contains all taxids under 'ID' and all probabilities under 'score' <br>
- peptonizer_results.png: bar plot of the peptonizer results showing the scores for the #'effects_in_plot' (see config parameters) highest scoring effects
  <br>

Additional files: <br>
- Intermediate results folders sorted by their prior value for all possible grid search parameter combinations
- effects_weights_dataframe.csv: csv file of all taxids that had at least one peptide map to them and their weight 
- pepgm_graph.graphml: graphml file of the graphical model (without convolution tree factors). Useful to visualize the graph structure and peptide-effect connections <br>
- sequence_scores_dataframe.csv: dataframe with petides, effects and scores used to create the graph <br>
- best_parameter.csv: file with best parameter <br>
- unipept_responses.json: response of unipept queries <br>
- effect_cluster_heads_dataframe: additional .csv file resulting from the clustering of effects by peptidome used for rbo<br>


<p align="right">(<a href="#top">back to top</a>)</p>


## Testing the Peptonizer
<!-- Testing -->

To test the Peptonizer2000 and see if it is set up correctly on your machine, we provide a test file under resources/test_files. This should be dowloaded automatically if you follow the installation instructions above. There are several test files from different metaproteomic samples. These are: <br>
- the samples S03, S05 and S11 of the [CAMPI study](https://www.nature.com/articles/s41467-021-27542-8) searched against a sample specific database using X!Tandem and MS2Rescore. The original files are available through [PRIDE under PXD023217](https://www.ebi.ac.uk/pride/archive/projects/PXD023217/). 
- the sample U1 of uneven communities from a [metaproteomic benchmark study by Kleiner](https://www.nature.com/articles/s41467-017-01544-x) searched against a sample specific database. The original files are available through [PRIDE under PXD006118](https://www.ebi.ac.uk/pride/archive/projects/PXD006118)
- the sample F07, a fecal sample, of the [CAMPI study](https://www.nature.com/articles/s41467-021-27542-8) searched against the integrated gene catalog for the human gut using X!Tandem and MS2Rescore. The original files are available through [PRIDE under PXD023217](https://www.ebi.ac.uk/pride/archive/projects/PXD023217/). 

To execute a test run of the Peptonizer2000 using the provided files: 
 
 1. Follow the installation instructions above
 2. In the config file, make sure to point to the test sample you want to use. By default, this is S03
 3. Start to peptonize with the command `snakemake --use-conda --cores 1`. If you have sufficient CPU and memory power available to your system, you can increase the amount of cores in order to speed up the workflow.

<p align="right">(<a href="#top">back to top</a>)</p>

## Running an analysis natively, without Snakemake or the website

The `peptonizer_analysis` folder provides three standalone Rust command-line tools that run the full
Peptonizer2000 pipeline (effect weighing, factor graph construction, a belief-propagation grid search, and
best-parameter selection) directly against TSV files, without installing Conda/Snakemake and without a browser.
You supply the peptide-effect relationships yourself, rather than starting from a raw peptide list — `protein_inference`
and `functional_analysis` then make no Unipept queries at all, since protein and function IDs are used as-is:

- `taxonomic_analysis` — infers taxonomic origin from peptide-to-taxon relationships (taxon IDs are normalized
  to species rank before weighing via a Unipept API call, so this tool needs network access)
- `protein_inference` — infers the source protein from peptide-to-protein relationships (no Unipept queries)
- `functional_analysis` — infers functional annotations from peptide-to-function relationships (no Unipept
  queries)

### Building

```bash
cd peptonizer_analysis
cargo build --release
```

The compiled binaries are written to `target/release/`.

### Input files

Each tool takes three tab-separated files with **no header row**:

| File | Flag | Columns |
|---|---|---|
| Relationships | `--peptide-taxa` / `--peptide-proteins` / `--peptide-functions` | `peptide`, `id` (an integer taxon/function ID, or a protein name for `--peptide-proteins`) |
| Scores | `--peptide-scores` | `peptide`, `score` (a float, e.g. from your search engine) |
| Counts | `--peptide-counts` | `peptide`, `count` (an integer PSM count) |

A peptide can appear on multiple rows of the relationships file to associate it with more than one ID.

### Usage

```bash
./target/release/taxonomic_analysis \
  --peptide-taxa peptide_taxa.tsv \
  --peptide-scores peptide_scores.tsv \
  --peptide-counts peptide_counts.tsv \
  --output taxonomic_analysis_results.tsv

./target/release/protein_inference \
  --peptide-proteins peptide_protein.tsv \
  --peptide-scores peptide_scores.tsv \
  --peptide-counts peptide_counts.tsv \
  --output protein_inference_results.tsv

./target/release/functional_analysis \
  --peptide-functions peptide_functions.tsv \
  --peptide-scores peptide_scores.tsv \
  --peptide-counts peptide_counts.tsv \
  --output functional_analysis_results.tsv
```

`--output` is optional and defaults to `<tool_name>_results.tsv` in the current directory. Run any tool with
`--help`/`-h` for a one-line usage reminder. Unlike the Snakemake workflow's `config.yaml`, the alpha/beta/prior
grid-search ranges are currently fixed per tool rather than user-configurable.

### Output

A TSV with a header row and one row per ID, sorted by descending posterior probability, e.g. for
`taxonomic_analysis`:

```
taxon_id	probability
2	0.98
816	0.42
```

(`protein_id`/`function_id` for the other two tools; `protein_inference` reports the original protein name
instead of an internal ID.)

<p align="right">(<a href="#top">back to top</a>)</p>

<!-- LICENSE -->
## License

Distributed under the Apache 2.0 License. See `LICENSE.txt` for more information.

<p align="right">(<a href="#top">back to top</a>)</p>


<!-- CONTACT -->
## Contact

Tanja Holstein - [@HolsteinTanja](https://twitter.com/HolsteinTanja) - tanja.holstein@ugent.be <br>
Pieter Verschaffelt - pieter.verschaffelt@ugent.be

<div align="center">
  <img src="https://raw.githubusercontent.com/compomics/Peptonizer2000/refs/heads/master/peptonizer_developers.jpeg" alt="Logo"  height="300">
</div>

<p align="right">(<a href="#top">back to top</a>)</p>


<!-- MARKDOWN LINKS & IMAGES -->
<!-- https://www.markdownguide.org/basic-syntax/#reference-style-links -->
[contributors-shield]: https://img.shields.io/github/contributors/BAMeScience/repo_name.svg?style=for-the-badge
[contributors-url]: https://github.com/BAMeScience/repo_name/graphs/contributors
[forks-shield]: https://img.shields.io/github/forks/BAMeScience/repo_name.svg?style=for-the-badge
[forks-url]: https://github.com/BAMeScience/repo_name/network/members
[stars-shield]: https://img.shields.io/github/stars/BAMeScience/repo_name.svg?style=for-the-badge
[stars-url]: https://github.com/BAMeScience/repo_name/stargazers
[issues-shield]: https://img.shields.io/github/issues/BAMeScience/repo_name.svg?style=for-the-badge
[issues-url]: https://github.com/BAMeScience/repo_name/issues
[license-shield]: https://img.shields.io/github/license/BAMeScience/repo_name.svg?style=for-the-badge
[license-url]: https://github.com/BAMeScience/repo_name/blob/master/LICENSE.txt
[linkedin-shield]: https://img.shields.io/badge/-LinkedIn-black.svg?style=for-the-badge&logo=linkedin&colorB=555
[linkedin-url]: https://linkedin.com/in/linkedin_username
[product-screenshot]: images/screenshot.png
[Next.js]: https://img.shields.io/badge/next.js-000000?style=for-the-badge&logo=nextdotjs&logoColor=white
[Next-url]: https://nextjs.org/
[React.js]: https://img.shields.io/badge/React-20232A?style=for-the-badge&logo=react&logoColor=61DAFB
[React-url]: https://reactjs.org/
[Vue.js]: https://img.shields.io/badge/Vue.js-35495E?style=for-the-badge&logo=vuedotjs&logoColor=4FC08D
[Vue-url]: https://vuejs.org/
[Angular.io]: https://img.shields.io/badge/Angular-DD0031?style=for-the-badge&logo=angular&logoColor=white
[Angular-url]: https://angular.io/
[Svelte.dev]: https://img.shields.io/badge/Svelte-4A4A55?style=for-the-badge&logo=svelte&logoColor=FF3E00
[Svelte-url]: https://svelte.dev/
[Laravel.com]: https://img.shields.io/badge/Laravel-FF2D20?style=for-the-badge&logo=laravel&logoColor=white
[Laravel-url]: https://laravel.com
[Bootstrap.com]: https://img.shields.io/badge/Bootstrap-563D7C?style=for-the-badge&logo=bootstrap&logoColor=white
[Bootstrap-url]: https://getbootstrap.com
[JQuery.com]: https://img.shields.io/badge/jQuery-0769AD?style=for-the-badge&logo=jquery&logoColor=white
[JQuery-url]: https://jquery.com 

