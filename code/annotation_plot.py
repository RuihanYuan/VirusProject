# Import necessary libraries
from BCBio import GFF
from dna_features_viewer import GraphicFeature, GraphicRecord

# Define a color mapping for different products
color_mapping = {
    "P0 protein": "lightpink",
    "RNA-dependent RNA polymerase": "magenta",
    "P1 protein": "lightgreen",
    "coat protein": "yellow",
    "aphid transmission protein":"orange",
    # Add more mappings as needed
}

# Default color if the product is not in the mapping
default_color = "grey"

# Initialize an empty list to hold GraphicFeature objects
graphic_features = []

# Parse the GFF3 file and create GraphicFeature objects
with open("/path/to/file.gff3", "r") as gff_input:
    for gff_record in GFF.parse(gff_input):
        for feature in gff_record.features:
            if feature.type == "CDS":
                product = feature.qualifiers.get('product', [''])[0]
                color = color_mapping.get(product, default_color)
                graphic_feature = GraphicFeature(
                    start=int(feature.location.start),
                    end=int(feature.location.end),
                    strand=feature.location.strand,
                    color=color,
                    label=product
                )
                graphic_features.append(graphic_feature)

# Create a GraphicRecord outside the loop, after all features have been added to the list
graphic_record = GraphicRecord(sequence_length=len(gff_record.seq), features=graphic_features)

# Plot the graphic record outside the loop
ax, _ = graphic_record.plot(figure_width=10)
ax.figure.tight_layout()
ax.figure.savefig("annotations_plot1.png")  # Save the figure to a file

