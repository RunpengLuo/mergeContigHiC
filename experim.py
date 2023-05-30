import sys, os
from utils import *


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
    _, contig_file, fwd_file, rve_file, reference = sys.argv

    mapq_list = [0, 10, 20, 30, 40]
    inflation_arr = [1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0]
    # inflation_arr = [2.0]

    for mapq in mapq_list:
        for inflation_rate in inflation_arr:
            print(f"----->Experiment on MAPQ: {mapq}, inflation rate: {inflation_rate}")
            # construct results
            print("* constructing model & clustering..")
            exec_str = f"python {main_py} {contig_file} {fwd_file} {rve_file} {mapq} {inflation_rate} > exp_Q{mapq}_I{inflation_rate}.1.log"
            System(exec_str)

            cluster_file = f"mcl_graph_Q{mapq}_I{inflation_rate}.cluster.txt"

            print("verifying the cluster")
            exec_str = f"python {verify_py} {contig_file} {cluster_file} {reference} > exp_Q{mapq}_I{inflation_rate}.2.log"
            System(exec_str)
        
        mat_file = f"mat_mapq_{mapq}.txt"
        # visualize the heatmap only if #entry leq 500
        if os.path.exists(mat_file):
            print("visualizing the heatmap")
            exec_str = f"python {visual_py} {mat_file} {empty} {contig_file} > exp_Q{mapq}.3.log"
            System(exec_str)
        else:
            print("skip heatmap construction, entries too large")



