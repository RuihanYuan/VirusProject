# !/usr/bin/env python3
# -*- coding: utf-8 -*-

import seaborn as sns
import numpy as np
from matplotlib import pyplot as plt, lines

def parse_depth(depth_input):
    """
    Parses a depth file to extract sequencing depth information for a single genome or contig.
    
    This function reads a depth file where each row contains a genome identifier, position, 
    and depth count. It accumulates depth counts into a list, ensuring each position is 
    accurately represented. The function is designed to handle files pertaining to a single 
    genome or contig; it will raise an exception if multiple genome identifiers are detected.
    
    Parameters:
    - depth_input (str): The path to the input depth file. Each line of the file should follow 
      the format: 'genome_id position depth_count'.
    
    Returns:
    - list: A list of integers representing the sequencing depth at each position of the genome.
            The list index corresponds to the genome position (0-indexed), and the value at each 
            index is the depth count for that position.
    
    Raises:
    - Exception: If the input file contains data for more than one genome or contig, indicating 
      it's not supported by this script.
    """
    depth = []
    references = set()
    with open(depth_input) as depth_object:
        for row in depth_object:
            genome_id, position, depth_count = row.split()
            references.add(genome_id)
            if len(references) > 1:
                raise Exception('This script only handles one genome - contig.')
            position = int(position) - 1
            depth_count = int(depth_count)
            if position >= len(depth):
                depth.extend([0] * (position - len(depth) + 1))
            depth[position] = depth_count
    return depth

def plot_depth(depth_report, output_name, plot_title, normalize=False, depth_cut_off=20):
    data = parse_depth(depth_report)
    y_label = "Normalized Depth" if normalize else "Depth"
    data = [xx / max(data) for xx in data] if normalize else data
    sns.set(color_codes=True)
    plt.title(plot_title)
    ax = plt.subplot(111)
    sns_plot = sns.lineplot(x=range(len(data)), y=data)
    sns_plot.set(xlabel='Genome Position (bp)', ylabel=y_label)
    if not normalize:
        ax.add_line(lines.Line2D([0, len(data) + 1], [depth_cut_off], color="r"))
    plt.savefig(output_name, bbox_inches='tight', dpi=400)
    plt.close()
    print("Done :)")

def main():
    plot_depth("/path/to/sample.depth", "/path/to/save/depth_raw.png", "plot title", False, 40)
    plot_depth("/path/to/sample.depth", "/path/to/save/depth_normalized.png", "plot title", True, 40)
    print("Plot generation complete!")

if __name__ == '__main__':
    main()
