# Read-depth and Product Annotation Plotting Tools

This repository contains three scripts for generating plots using Python: `dplot`, `aplot`, and `daplot`. These tools are designed to create read depth plots, CDS product annotation plots, and combined depth and annotation plots, respectively.

**Note:** This program has experienced a change of developer and is now maintained by Kristian Stevens and Henry Li at [Foundation Plant Services](https://fps.ucdavis.edu/index.cfm) (FPS) at [UC Davis](https://www.ucdavis.edu/) following its initial development by Ruihan Yuan. Additionally, one of the dependencies, [BCBio](https://github.com/bcbio/bcbio-nextgen), stopped receiving support as of 08-16-2024 and is no longer installable from pip. Therefore, the installation instructions now include steps to install dependencies using `conda`.

Jump to [Installation](#installation), [Preparation](#preparation), or [Usage](#usage).

## Overview

In the realm of Next Generation Sequencing (NGS), depth, or coverage, refers to the number of times a specific sequence has been read and it plays a pivotal role in determining the reliability of genomic data analysis. A higher depth at particular genome positions suggests that those regions have been more thoroughly sequenced. Having depth information visually placed against the genome annotation not only shows the reliability of specific areas but also helps derive more accurate conclusions by offering a detailed view of sequence coverage across various regions. 

The goal of this tool is to generate depth and annotation plots that visually enhance such interpretations by vertically placing the plots together, showcasing the number of reads and their corresponding positions on the genome at the same time. Such visualizations are crucial for assessing the uniformity of coverage and pinpointing areas with exceptionally high or low sequencing depth, which may help in the identification of well-sequenced regions versus areas that may require additional attention. All the graphs generated are in SVG format, which should fit any downstream usages.

## Installation

The scripts are written in **Python**, which you will need to run them. Managing Python installations is never easier nowadays with **Anaconda**. We will also use Anaconda to install and manage all of the dependancies. 

To install Anaconda on your machine, please follow the [Anaconda Installation Documentation](https://docs.anaconda.com/anaconda/install/).

After installing Anaconda, to verify you have installed it correctly, use:

```bash
conda --version
```

and if your machine has Anaconda installed, you should see its version info.

To set up the rest of the dependencies, please follow the following steps:

Create a new conda environment:

```bash
conda create -n daplot python=3.12
```

Activate the conda environment:

```bash
conda activate daplot
```

Install the required dependencies:

```bash
conda install -c bioconda -c conda-forge bcbio-gff bwa dna-features-viewer matplotlib pyyaml samtools seaborn
```

After following the installation, you should see all of the dependencies installed using:

```bash
conda list
```

Now, you can clone this project onto any directory on your machine and navigate to the bin folder, where you should see the scripts.

To test all the scrips are ready, use:

```bash
./daplot --version
```

and the current version number should be displayed.

## Preparation

To prepare your raw reads into a depth file, we need to map it to the reference genome and sort it. We recommend the widely-used tools **BWA** for the alignment and **Samtools** for further processing. This process is designed with a focus on a single contig genome. Follow the steps below to prepare your data:

Index your reference genome using BWA:
   
```bash
bwa index {REFERENCE_GENOME.FASTA}
```

Align reads to the reference genome:

```bash
bwa mem -t {NUMBER_OF_THREADS_INTEGER} {REFERENCE_GENOME.FASTA} {READS_INPUT} | samtools sort -o {ALIGNMENT}.bam
```

Generate a depth file using samtools:
   
```bash
samtools depth {ALIGNMENT}.bam > {GENOME}.depth
```

**Note**: Samtools `depth` command counts positions starting from 1, which might become relevant when interpreting your results.

Now, your depth file is ready as inputs for our scripts.

## Usage

For all three scripts, use:

```bash
./scriptname -h
```

to see help menu and detailed argument info.

### dplot

The `dplot` script generates depth plots from a depth file.

**Command:**

```bash
./dplot -i /path/to/sample.depth -o /path/to/save/depth.svg -t "PLOT_TITLE" -n -c {CUTOFF_INTEGER}
```

**Arguments:**

- `-i` or `--input`: Path to the input depth file (required).
- `-o` or `--output`: Path to save the output SVG file (required).
- `-t` or `--title`: Title of the plot (default: “Depth Plot”).
- `-n` or `--normalize`: Normalize the depth values (optional).
- `-c` or `--cutoff`: Depth cutoff value (default: 20).
- `-v` or `--version`: Display the version of the script.

### aplot

The `aplot` script generates color-coded annotation plots from a GFF file.

**Command:**

```bash
./aplot -i /path/to/sample.gff3 -o /path/to/save/annotation.svg -y /path/to/color_settings.yml
```

**Arguments:**

- `-i` or `--input`: Path to the input GFF file (required).
- `-o` or `--output`: Path to save the output SVG file (required).
- `-y` or `--color_settings`: Path to the color settings YAML file (optional, defaults to program settings if left blank).
- `-v` or `--version`: Display the version of the script.

To make changes to the color map for each annotated product, please do so in the **color_setting.yml** file in the bin folder and pass it in to override the default settings using the -c argument. Supports hex color codes but make sure to place them in quotes.

Here is an example color_setting.yml file you could use:

```yml
color_mapping:
  P0 protein: '#88B04B'
  RNA-dependent RNA polymerase: '#92A8D1'
  P1 protein: '#6B5B95'
  coat protein: '#FF6F61'
  aphid transmission protein: '#F7CAC9'

default_color: '#9F9F9F'
```

GFF3 files for your genome can be found and downloaded from [NCBI](https://www.ncbi.nlm.nih.gov/). Read more about GFF3 [here](https://www.ncbi.nlm.nih.gov/datasets/docs/v1/reference-docs/file-formats/about-ncbi-gff3/).

### daplot

The `daplot` script generates combined depth and annotation plots sharing a single x axis from a depth file and a GFF file.

**Command:**

```bash
./daplot -id /path/to/sample.depth -ig /path/to/sample.gff3 -o /path/to/save/combined.svg -n -c {CUTOFF_INTEGER} -y /path/to/color_settings.yml
```

**Arguments:**

- `-id` or `--input_depth`: Path to the input depth file (required).
- `-ig` or `--input_gff`: Path to the input GFF file (required).
- `-o` or `--output`: Path to save the output SVG file (required).
- `-n` or `--normalize`: Normalize the depth values (optional).
- `-c` or `--cutoff`: Depth cutoff value (default: 20).
- `-y` or `--color_settings`: Path to the color settings YAML file (optional, defaults to program settings if left blank).
- `-v` or `--version`: Display the version of the script.

## License

This project is licensed under the MIT License.
