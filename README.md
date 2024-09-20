# VirusProject
Data visualization and analysis using python

## Depth Analysis ##

### Overview ###
In the realm of Next Generation Sequencing (NGS), depth, or coverage, plays a pivotal role in determining the reliability of genomic data analysis.
**Depth** refers to the number of times a specific sequence has been read. A higher depth at particular genome positions suggests that those regions have been more thoroughly sequenced.
This not only implies greater reliability in those areas but also significantly enhances the overall accuracy of genome analysis by offering a detailed view of sequence coverage across various regions.

### Purpose of the Tool ###
The primary goal of this tool is to generate insightful plots that visually represent the relationship between the number of reads and their corresponding positions on the genome.
Such visualizations are crucial for assessing the uniformity of coverage and pinpointing areas with exceptionally high or low sequencing depth.

### Key Features ###

- **Visualization of Sequencing Depth**: Generates comprehensive plots that illustrate the distribution and variance of sequencing depth across the genome.
This facilitates a deeper understanding of coverage patterns, aiding in the identification of well-sequenced regions versus areas that may require additional attention or reads.

- **Enhanced Accuracy and Clarity**: By providing a visual representation of sequence coverage, this tool aids researchers in making informed decisions about their data analysis strategies, 
ensuring that genomic analyses are based on accurately sequenced regions.

### Technical Foundation ###

The methodology and approach implemented in this tool were inspired by and adapted from a tutorial available on [One Stop Data Analysis](https://onestopdataanalysis.com/depth-plot/). 
This foundational guide to depth plot creation offered the initial framework, which was further refined and customized to suit the specific needs of NGS data analysis.

### Getting Started ###

To utilize this tool in your NGS data analysis projects, please refer to the following steps:

1. **Installation**: Instructions on setting up the Python environment and installing necessary libraries.
2. **Usage**: Detailed guide on how to run the program, including command-line arguments and input file specifications.
3. **Examples**: Sample input files along with corresponding output plots to demonstrate the tool's capabilities and expected results.

For more detailed information, please refer to the [Installation](#installation) and [Usage](#usage) sections below.


## Installation ##

### Prerequisites ###

Before using this tool, you'll need to have Python installed on your system. This program is compatible with Python 3.x versions. If you don't have Python installed, follow the steps below to set it up. 

### Installing Python ###

1. **Download Python**: Visit the official Python website at [python.org](https://www.python.org/) and download the latest version of Python 3.x for your operating system (Windows, macOS, or Linux/Unix).

2. **Run the Installer**: Launch the downloaded installer. Ensure to check the box that says **"Add Python 3.x to PATH"** to make Python accessible from the command line.

3. **Verify Installation**: Open your command line interface (CLI) and type `python --version` (or `python3 --version` on some Linux distributions). 
You should see the Python version you installed displayed, confirming the successful installation.

### Dependencies ###
With Python step up, the next step is to install the necessary libraries and dependencies fo the tool.
Such dependencies can be installed using 'pip', Python's package installer.
```
pip install seaborn matplotlib numpy
```

### Usage ###

#### Generating Genome Depth File ####

To generate the input depth file required by the plot generation program, we employ the widely-used tools **BWA** for aligning reads to a reference genome and **Samtools** for processing sequence alignments. This process is designed with a focus on a single contig genome. Follow the steps below to prepare your data:

1. **Format the Reference Genome with BWA**:
   
   First, you need to index your reference genome using BWA. This step prepares the genome for alignment.
   
   ```
   bwa index {REFERENCE_GENOME.FASTA}
   ```

2. **Align Reads to the Reference Genome**:
   
   Next, align your sequencing reads to the reference genome. This step uses BWA to perform the alignment and Samtools to sort the resulting BAM file. Adjust `-t 16` according to the number of threads you wish to use for the alignment process.
   
   ```
   bwa mem -t 16 {REFERENCE_GENOME.FASTA} {READS_INPUT} | samtools sort -o {ALIGNMENT}.bam
   ```

3. **Generate Depth File with Samtools**:
   
   With the alignments in a sorted BAM file, you can now generate the depth report using Samtools. This report will serve as the input for the plot generation program.
   
   ```
   samtools depth {ALIGNMENT}.bam > {GENOME}.depth
   ```

**Note**: Samtools `depth` command counts positions starting from 1. This is important to remember when interpreting your results.

4. **Plot Generation**

   Now, we are ready to use the python script to generate the depth plot. Please refer to the `depth_plot_final.py` in code file.

#### Generating genome Annotation Plot ####

   To enhance our understanding of genomic structures, the DNA annotation plot visualizes different genomic features, such as genes, exons, and regulatory regions, each marked according to its biological function. This visualization is crucial for tasks ranging from gene discovery to regulatory element identification.
   
   **Overview of DNAFeaturesViewer**
   
   We use **DNAFeaturesViewer**, a versatile Python library designed for plotting DNA features. It can interpret data from GenBank, GFF files, or Biopython SeqRecords, allowing for a broad application in genomic studies. Detailed setup instructions are available [here](https://edinburgh-genome-foundry.github.io/DnaFeaturesViewer/).
   
   **Why Annotation Plots Matter**
   
   Annotation plots are essential for visualizing the physical layout of genes and other features on a genome. They are especially useful in comparative genomics, gene editing studies, and in the annotation of newly sequenced genomes. Such plots help researchers quickly identify regions of interest and assess the completeness and accuracy of genomic annotations.
   
   **Advantages of Combining Annotation with Depth Plots**
   
   By integrating annotation plots with depth analysis, researchers can not only visualize the functional layout of the genome but also see how thoroughly each segment has been sequenced. This combined visualization approach provides a comprehensive view of both genomic structure and sequencing coverage. It is particularly valuable for identifying less sequenced or missing functional regions that require further attention. Thus, combined plots serve as an excellent tool for refining research areas within the genome, offering a robust method for investigating and addressing gaps in functional regions.

### Dependencies ###

```
pip install BCBio dna_features_viewer
```

1. **Input File Download**
   
   We use GFF3 file as input in the annotation plot generation step. GFF3 files for desired genome can be found and downloaded from ncbi.
   
2. **Customizing the Plot**

   You can easily customize your plots by modifying the color_mapping dictionary in the script. This allows you to highlight different features in distinct colors, enhancing the plot's readability and focusing on features of interest.
   ```python
   # Customizing colors for genomic features
     color_mapping = {
       'gene': 'blue',
       'exon': 'green',
       'promoter': 'red'
     }
   
3. **Accessing the Scripts**

   The scripts for generating these plots are located in the code folder. You can refer to the following files for the specific implementations:

   * Annotation Plot: Check out [annotation_plot.py](./code/annotation_plot.py) for the code responsible for generating annotation plots.
   * Combined Plot: For the combined annotation and depth plots, see [combined_plot.py](./code/combined_plot.py).
