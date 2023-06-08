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
def gen_matrix(cluster_file: list, ref_ids: list, assign_dict: dict):
    clusters = get_clusters(cluster_file)
    num_clusters = len(clusters)
    full_ref = list(ref_ids)
    full_ref.append("NA")
    mat = dict.fromkeys(range(num_clusters), None)
    for i, cids in enumerate(clusters):
        mat[i] = dict.fromkeys(full_ref, 0)
        for cid in cids:
            if cid not in assign_dict:
                # contig cannot be mapped to any ground-truth
                mat[i]["NA"] += 1
            else:
                ref_no, seg_len, ref_len, matching, mapped, quality = assign_dict[cid]
                mat[i][ref_no] += 1
        print(f">>>{i}-th cluster\n{mat[i].keys()}\n{mat[i].values()}")

    return mat

def mat2csv(mat: dict, ref_ids: list, prefix: str):
    csv_file = prefix + "mat.csv"
    System("echo "" > " + csv_file)
    with open(csv_file, "w") as csv_fd:
        # header
        s = "cluster vs reference," + ",".join(ref_ids)+",NA\n"
        csv_fd.write(s)
        for i in range(len(mat)):
            row = mat[i]
            arr = [str(i)+";"+str(sum(row.values()))]
            for ref_id in ref_ids:
                arr.append(str(row[ref_id]))
            arr.append(str(row["NA"]))
            csv_fd.write(",".join(arr) + "\n")
        csv_fd.close()
    return

def eval_cluster(mat: dict, ref_ids: list, total_elem: int):
    # compute precision
    num_non_sense = sum([row["NA"] for row in mat.values()])
    classified = sum([sum(row.values()) for row in mat.values()])
    max_utgs_sum = sum([max(row.values()) for row in mat.values()])
    precision = float(max_utgs_sum) / classified
    print(max_utgs_sum)
    print(classified)
    print("Precision: ", round(precision, 3))
    # compute recall
    recall_num = 0
    for ref_id in ref_ids:
        col_max = 0
        for i, row in mat.items():
            col_max = max(col_max, row[ref_id])
        recall_num += col_max
    recall = float(recall_num) / total_elem
    print("Recall: ", round(recall, 3))

    # compute F1-score
    f1 = 2 * (float(precision * recall) / (precision + recall))
    print("F1: ", round(f1, 3))

    # compute Adjusted Rand Index (ARI)
    # measure of similarity between the binning result and its actual grouping. It is calculated as follows.
    return precision, recall, f1

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

def eval_wrapper(contig_file, cluster_file, ref_file, prefix):
    refs = get_utgs(ref_file)
    ref_ids = list(refs.keys())
    utgs = get_utgs(contig_file)

    mapped_contig_paf = "mapped.contigs.paf"
    print("getting contig-ref alignment")
    if not os.path.exists(mapped_contig_paf):
        System("minimap2 --secondary=no {0} {1} > {2}".format(ref_file, contig_file, mapped_contig_paf))
    else:
        print("mapped contigs paf exists")

    assign_dict = get_contig_assignment(mapped_contig_paf, utgs)
    mat = gen_matrix(cluster_file, ref_ids, assign_dict)
    mat2csv(mat, ref_ids, prefix)
    precision, recall, f1 = eval_cluster(mat, ref_ids, len(utgs))

    return precision, recall, f1

# if __name__ == "__main__":
#     if len(sys.argv) != 4:
#         print(f"{sys.argv[0]} <contigs .fasta> <cluster .txt> <reference .fasta>")
#         sys.exit(0)
    
#     _, contig_file, cluster_file, ref_file = sys.argv

#     refs = get_utgs(ref_file)
#     ref_ids = list(refs.keys())
#     utgs = get_utgs(contig_file)
#     utg_ids = list(utgs.keys())

#     mapped_contig_paf = "mapped.contigs.paf"
#     print("getting contig-ref alignment")
#     if not os.path.exists(mapped_contig_paf):
#         System("minimap2 --secondary=no {0} {1} > {2}".format(ref_file, contig_file, mapped_contig_paf))
#     else:
#         print("mapped contigs paf exists")

#     assign_dict = get_contig_assignment(mapped_contig_paf, utgs)
#     mat = gen_matrix(cluster_file, ref_ids, assign_dict)
#     eval_cluster(mat, ref_ids)

#     # check_cluster(cluster_file, assign_dict, refs)
#     sys.exit(0)


    