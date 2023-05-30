import sys, os

from utils import *


def get_clusters(cluster_file: str):
    clusters = []
    with open(cluster_file, "r") as cfd:
        for line in cfd:
            cluster = line.strip().split()
            clusters.append(cluster)
        cfd.close()
    return clusters

"""
    generate a K*S matrix M with K number of clusters and S number of references
    with entry M[k][s] denotes the relative score be assigned to s-reference on k-th cluster
"""
def gen_matrix(clusters: list, ref_ids: list, assign_dict: dict):
    num_clusters = len(clusters)
    num_ref = len(ref_ids)
    mat = dict.fromkeys(range(num_clusters), None)
    non_sense = []
    for i, cids in enumerate(clusters):
        mat[i] = dict.fromkeys(ref_ids, 0)
        for cid in cids:
            if cid not in assign_dict:
                # contig cannot be mapped to any ground-truth
                non_sense.append(cid)
            else:
                ref_no, seg_len, ref_len, matching, mapped, quality = assign_dict[cid]
                mat[i][ref_no] += 1
        print(f">>>{i}-th cluster\n{mat[i].keys()}\n{mat[i].values()}")

    return mat, non_sense

def eval_cluster(mat: dict, ref_ids: list, non_sense: list):
    num_non_sense = len(non_sense)
    # compute precision
    total_utgs = sum([sum(row.values()) for row in mat.values()])
    max_utgs_sum = sum([max(row.values()) for row in mat.values()])
    precision = float(max_utgs_sum) / total_utgs
    print("Precision: ", round(precision, 3))
    # compute recall
    recall_num = 0
    for ref_id in ref_ids:
        col_max = 0
        for i, row in mat.items():
            col_max = max(col_max, row[ref_id])
        recall_num += col_max
    recall = float(recall_num) / (total_utgs + num_non_sense)
    print("Recall: ", round(precision, 3))

    # compute F1-score
    F1 = 2 * (float(precision * recall) / (precision + recall))
    print("F1: ", round(precision, 3))

    # compute Adjusted Rand Index (ARI)
    # measure of similarity between the binning result and its actual grouping. It is calculated as follows.
    return

def get_contig_assignment(mapped_paf: str, utgs: list):
    assign_dict = {}
    with open(mapped_paf, "r") as pfd:
        for Line in pfd:
            splited = Line.strip().split("\t")
            seg_no = str(splited[0])
            seg_len = int(splited[1])
            assert seg_no in utgs
            ref_no = str(splited[5])
            ref_len = int(splited[6])
            matching = int(splited[9])
            mapped = int(splited[10])
            quality = int(splited[11])
            assign_dict[seg_no] = (ref_no, seg_len, ref_len, matching, mapped, quality)
        pfd.close()
    return assign_dict

# def check_cluster(cluster_file: str, assign_dict: dict, refs: dict):
#     clusters = get_clusters(cluster_file)
    
#     for i, cluster in enumerate(clusters):
#         print("-----------------------------------------")
#         print("-----------------------------------------")
#         print(f"-Current checking {i}th cluster with {len(cluster)} entries")
#         ref_assigns = {}
#         for ref_id in refs.keys():
#             ref_assigns[ref_id] = []
#         ref_assigns['na'] = []
#         for cid in cluster:
#             if cid in assign_dict:
#                 (ref_no, seg_len, ref_len, matching, mapped, quality) = assign_dict[cid]
#                 ref_assigns[ref_no].append((cid, ref_no, seg_len, ref_len, matching, mapped, quality))
#             else:
#                 ref_assigns['na'].append(cid)
#         for rid, rstatus in ref_assigns.items():
#             if rstatus != []:
#                 if rid != 'na':
#                     print(f"*ref: {rid} with size: {len(refs[rid])}, contig count: {len(rstatus)}")
#                     for (cid, ref_no, seg_len, ref_len, matching, mapped, quality) in rstatus:
#                         print(f"-->contig: {cid} with size: {seg_len}, match/total mapping: {matching} / {mapped}, qscore: {quality}")
#                 else:
#                     print(f"*NA, contig count: {len(rstatus)}")
#                     for cid in rstatus:
#                         print(f"-->contig: {cid}")
#         (max_rid, max_status) = max(ref_assigns.items(), key=lambda tple: len(tple[1]))
#         print(f"majority ref: {max_rid}, ({len(max_status)}/{len(cluster)})")
#     return

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print(f"{sys.argv[0]} <contigs .fasta> <cluster .txt> <reference .fasta>")
        sys.exit(0)
    
    _, contig_file, cluster_file, ref_file = sys.argv

    refs = get_utgs(ref_file)
    ref_ids = list(refs.keys())
    utgs = get_utgs(contig_file)
    utg_ids = list(utgs.keys())

    mapped_contig_paf = "mapped.contigs.paf"
    print("getting contig-ref alignment")
    if not os.path.exists(mapped_contig_paf):
        System("minimap2 --secondary=no {0} {1} > {2}".format(ref_file, contig_file, mapped_contig_paf))
    else:
        print("mapped contigs paf exists")

    clusters = get_clusters(cluster_file)
    assign_dict = get_contig_assignment(mapped_contig_paf, utgs)
    mat, non_sense = gen_matrix(clusters, ref_ids, assign_dict)
    eval_cluster(mat, ref_ids, non_sense)

    # check_cluster(cluster_file, assign_dict, refs)
    sys.exit(0)


    