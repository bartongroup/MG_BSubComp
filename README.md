# Transcriptomic Changes in *ΔcomA* and *ΔrapA* Mutants in *Bacillus subtilis*

**Maintainer**: Marek Gierlinski
**Collaborators**: Jonathan Griffin, Nicola Stanley-Wall

## Usage

The software in this repository is designed to be run in two stages. The first stage, including FASTQ quality control, cleaning, and mapping, is run on a computing cluster. The second stage, including all downstream analysis, is designed to be run in R on a local laptop or desktop.

### On a Linux Cluster

Copy the contents of the `rna_seq` folder to the cluster filesystem. Change to the `rna_seq` directory and create and activate a `conda` environment:

```bash
cd rna_seq
conda create --name bsub_comp --file spec-file.txt
conda activate bsub_comp
```

Ensure that the FASTQ files are in the `./fastq` subdirectory.

Run Snakemake:

```bash
./run_snake.sh
```

This will trim adapters, perform quality control, download genome files, map reads to the reference, and count reads per gene.
**Note**: This shell script contains a call to *snakemake* with arguments configured for our Sun Grid Engine. If you're using a different cluster setup, you will need to modify the arguments accordingly.

### On a Local Laptop/Desktop

Once Snakemake has finished, we recommend using Positron or RStudio for downstream analysis. If you’re working on a different machine (e.g., I run Positron on a laptop), you will need to copy some data using:

```bash
./scripts/get_data.sh
```

**Note**: Before using this script, edit it to change the remote location, as it currently points to my directory.

Once the mapping results are available on your local machine, open an R console in Positron or RStudio and create the environment using `renv`:

```r
install.packages("renv")
renv::restore()
```

This will install all the required packages. If anything is missing, install additional packages locally using `renv::install()`.

Once the environment is set up, run the `targets` pipeline:

```r
targets::tar_make()
```

This will perform all calculations, generate data objects, figures, and tables (as [targets](https://books.ropensci.org/targets/) objects), and output CSV files in the `./tab` directory. An HTML report will be generated in the `./doc` directory.
