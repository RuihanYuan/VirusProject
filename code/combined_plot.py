
from dna_features_viewer import BiopythonTranslator, GraphicRecord, GraphicFeature
import seaborn as sns
import numpy as np
from matplotlib import pyplot as plt, lines
from BCBio import GFF


def parse_depth(depth_input):
    #initialize empty list to depth value at each position in the genome
    depth = []
    #to store unique genome ID, set() ensure each only being stored once
    references = set()
    with open(depth_input) as depth_object:
        for row in depth_object:
            genome_id, position, depth_count = row.split()
            references.add(genome_id)
            #if multiple genome IDs has been found, raise error since we only supposed to deal with 1 genome
            if len(references) > 1:
                raise Exception('This script only handles one genome - contig.')
            #python index starts at 1, so convert to int and subtract 1 to accommodate
            position = int(position) - 1
            #num of reads covers a specific position in the genome
            depth_count = int(depth_count)
            if position >= len(depth):
                #extend the litst if position is outside the range
                #calculation of the extended position: the desired position - length of current list + 1(python index starts at 0)
                depth.extend([0] * (position - len(depth) + 1))
                #update the depth count at the current position
            depth[position] = depth_count
    return depth

def plot_depth(ax,depth_report, normalize=False, depth_cut_off=20):
    data = parse_depth(depth_report)
    y_label = "Normalized Coverage" if normalize else "Coverage"
    data = [xx / max(data) for xx in data] if normalize else data
    sns.set(color_codes=True)
    sns.lineplot(x=range(len(data)), y=data,ax=ax)
    ax.set(xlabel='Genome Position (bp)', ylabel=y_label)
    if not normalize:
        ax.add_line(lines.Line2D([0, len(data) + 1], [depth_cut_off], color="r"))
    print("Plot added to the subplot.")

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(6, 3.5), sharex=True, gridspec_kw={"height_ratios": [3, 1]},) 
# Open the GFF3 file and parse it
with open("/path/to/file/.gff3", "r") as gff_input:
    for gff_record in GFF.parse(gff_input):
        # Create a list of GraphicFeature objects from the CDS features in the SeqRecord
        graphic_features = [
            GraphicFeature(
                start=int(feature.location.start),
                end=int(feature.location.end),
                strand=feature.location.strand,
                color="lightgreen",
                label=feature.qualifiers.get('product', [''])[0]
            )
            for feature in gff_record.features if feature.type == "CDS"
        ]

        # Create a GraphicRecord with the list of GraphicFeature objects
        graphic_record = GraphicRecord(sequence_length=len(gff_record), features=graphic_features)

        # Plot the graphic record on the first subplot
        graphic_record.plot(ax=ax1)



plot_depth(ax2, '/path/to/file/.depth', normalize=False, depth_cut_off=20)
plt.tight_layout()
plt.show()
fig.savefig("combined_plot_gff.png")