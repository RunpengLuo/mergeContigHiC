import sys, os
from utils import *

from verify_cluster import eval_wrapper

if __name__ == "__main__":
    exe_prefix = os.path.dirname(os.path.realpath(__file__))
    main_py = exe_prefix + "/main.py"
    verify_py = exe_prefix + "/verify_cluster.py"
    visual_py = exe_prefix + "/visual.py"
    empty = "tmp"
    print("base executable directory: " + exe_prefix)

    if len(sys.argv) != 5:        
        print(f"{sys.argv[0]} <contigs> <hic_forward> <hic_reverse> <reference>")
        sys.exit(1)
    _, contig_file, fwd_file, rve_file, ref_file = sys.argv

    mapq_list = [0, 10, 20, 30, 40]
    inflation_arr = [1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0]
    # inflation_arr = [2.0]

    for mapq in mapq_list:
        for inflation_rate in inflation_arr:
            print(f"----->Experiment on MAPQ: {mapq}, inflation rate: {inflation_rate}")
            cluster_file = f"mcl_graph_Q{mapq}_I{inflation_rate}.cluster.txt"
            if os.path.exists(cluster_file):
                print("skip alignment & clustering for current params, done.")
            else:
                # construct results
                print("* constructing model & clustering..")
                exec_str = f"python {main_py} {contig_file} {fwd_file} {rve_file} {mapq} {inflation_rate} > exp_Q{mapq}_I{inflation_rate}.1.log"
                System(exec_str)

    
    precision_file = "precision.csv"
    System("echo "" > " + precision_file)
    recall_file = "recall.csv"
    System("echo "" > " + recall_file)
    f1_file = "f1.csv"
    System("echo "" > " + f1_file)
    prec_fd = open(precision_file, "w")
    rec_fd= open(recall_file, "w")
    f1_fd = open(f1_file, "w")

    header = "inflation_rate vs mapq," + ",".join([str(m) for m in mapq_list]) + "\n"
    prec_fd.write(header)
    rec_fd.write(header)
    f1_fd.write(header)
    # row: inflation rate
    # column mapq
    for inflation_rate in inflation_arr:
        s_prec = [str(inflation_rate)]
        s_rec = [str(inflation_rate)]
        s_f1 = [str(inflation_rate)]
        for mapq in mapq_list:
            cluster_file = f"mcl_graph_Q{mapq}_I{inflation_rate}.cluster.txt"
            prefix = f"eval_Q{mapq}_I{inflation_rate}"
            print(f">>>>> verifying the cluster om MAPQ: {mapq}, inflation rate: {inflation_rate}")
            (precision, recall, f1) = eval_wrapper(contig_file, cluster_file, ref_file, prefix)
            print("Precision: ", precision)
            print("Recall: ", recall)
            print("F1 score: ", f1)

            s_prec.append(str(precision))
            s_rec.append(str(recall))
            s_f1.append(str(f1))

        prec_fd.write(",".join(s_prec) + "\n")
        rec_fd.write(",".join(s_rec) + '\n')
        f1_fd.write(",".join(s_f1) + "\n")
    
    prec_fd.close()
    rec_fd.close()
    f1_fd.close()
    # heatmap plot
    for mapq in mapq_list:
        mat_file = f"mat_mapq_{mapq}.txt"
        # visualize the heatmap only if #entry leq 500
        if os.path.exists(mat_file):
            print("visualizing the heatmap")
            exec_str = f"python {visual_py} {mat_file} {empty} {contig_file} > exp_Q{mapq}.3.log"
            System(exec_str)
        else:
            print("skip heatmap construction, entries too large")



