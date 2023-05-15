import os
import sys
import markov_clustering as mc
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import numpy as np
from utils import *

def load_mat(mat_file):
    mat = None
    with open(mat_file, "r") as mat_fd:
        lines = mat_fd.readlines()
        len_entries = len(lines[0].strip().split(","))
        mat = np.zeros((len_entries, len_entries), dtype=np.int32)
        for i, line in enumerate(lines):
            for j, e in enumerate(line.strip().split(",")):
                  mat[i][j] = int(e)
            
        mat_fd.close()
    return mat

def draw_hic_map(mat, utg_ids):
    print("draw hic map")
    fig, ax = plt.subplots(figsize=(10,10))
    mat_scaled = np.log(mat)
    #  im = plt.imshow(mat, cmap='hot', interpolation='nearest', norm=LogNorm(vmin=0.01, vmax=1))
    im = plt.imshow(mat_scaled, cmap='hot', norm=LogNorm())
    ax.set_xticks(ticks=np.arange(0, len(utg_ids), 1), labels=utg_ids, rotation=90)
    ax.set_yticks(ticks=np.arange(0, len(utg_ids), 1), labels=utg_ids)

    fig.colorbar(im)
    plt.title("HiC contact map")
    plt.savefig(mat_file + ".png")

def draw_cluster(mat, utg_ids):
    print("draw cluster")
    print(mat)
    result = mc.run_mcl(mat, inflation=2.0)
    clusters = mc.get_clusters(result)
    # pos=positions, 
    mc.draw_graph(mat, clusters, node_size=50, with_labels=False, edge_color="silver")
    plt.show()


if __name__ == "__main__":
    if len(sys.argv) != 4:
        print(f"{sys.argv[0]} <mat.txt> <cluster.txt> <contigs.fasta>")
        sys.exit(0)
    
    _, mat_file, cluster_file, contig_file = sys.argv

    mat = load_mat(mat_file)
    utgs = get_edges(contig_file)
    utg_ids = list(utgs.keys())

    draw_hic_map(mat, utg_ids)
    # draw_cluster(mat, utg_ids)

    sys.exit(0)