import sys, os

from utils import *


def get_contig_assignment(mapped_paf: str, utgs: list):
    assign_dict = {}
    with open(mapped_paf, "r") as pfd:
        for Line in pfd:
            splited = Line.strip().split("\t")
            seg_no = str(splited[0])
            assert seg_no in utgs
            ref_no = str(splited[5])
            assign_dict[seg_no] = ref_no
        pfd.close()
    return assign_dict

def check_cluster(cluster_file: str, assign_dict: dict, ref_ids: list):
    clusters = []
    with open(cluster_file, "r") as cfd:
        for line in cfd:
            cluster = line.strip().split()
            clusters.append(cluster)
        cfd.close()
    
    for i, cluster in enumerate(clusters):
        print("-----------------------------------------")
        print("-----------------------------------------")
        print(f"-Current checking {i}th cluster with {len(cluster)} entries")
        ref_assigns = dict.fromkeys(ref_ids, 0)
        ref_assigns['na'] = 0
        for cid in cluster:
            if cid in assign_dict:
                ref_assigns[assign_dict[cid]] += 1
            else:
                ref_assigns['na'] += 1
        for rid, rcount in ref_assigns.items():
            if rcount != 0:
                print(f"ref: {rid}, count: {rcount}")
        (max_rid, max_rcount) = max(ref_assigns.items(), key=lambda tple: tple[1])
        print(f"majority ref: {max_rid}, ({max_rcount}/{len(cluster)})")
    return

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print(f"{sys.argv[0]} <contigs .fasta> <cluster .txt> <reference .fasta>")
        sys.exit(0)
    
    _, contig_file, cluster_file, ref_file = sys.argv

    refs = get_edges(ref_file)
    ref_ids = list(refs.keys())
    utgs = get_edges(contig_file)
    utg_ids = list(utgs.keys())

    mapped_contig_paf = "mapped.contigs.paf"
    print("getting contig-ref alignment")
    if not os.path.exists(mapped_contig_paf):
        System("minimap2 --secondary=no {0} {1} > {2}".format(ref_file, contig_file, mapped_contig_paf))
    else:
        print("mapped contigs paf exists")

    assign_dict = get_contig_assignment(mapped_contig_paf, utgs)
    check_cluster(cluster_file, assign_dict, ref_ids)
    sys.exit(0)


    