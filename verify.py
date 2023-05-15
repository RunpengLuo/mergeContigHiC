import os
import sys

import numpy

def System(command):
    return os.system(command)

def list_to_string(ids: list, s=""):
    string = s + " - " if s != "" else ""
    for id in ids:
        string += str(id) + ", "
    return string[:-2] if len(string) >= 2 else ""

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("{0} <contigs> <ref>".format(sys.argv[0]))
        sys.exit(1)
    
    _, contig_file, ref_file = sys.argv
    System("minimap2 -x map-ont --secondary=no {0} {1} > contig2ref.paf".format(
        ref_file, contig_file
    ))

    ref_dict = {}
    utg_scores = {}
    utg_mapq = {}
    with open("contig2ref.paf", "r") as paf_fd:
        for Line in paf_fd:
            splited = Line.strip().split("\t")
            seg_no = str(splited[0])
            seg_len = int(splited[1])
            ref_no = str(splited[5])
            nmatches = int(splited[9])
            match_region = int(splited[10])
            mapq = int(splited[11])
            match_rate = float(nmatches)/match_region
            # tags = splited[12:]
            # nm = int([s.split(":")[2] for s in tags if s.startswith("NM")][0])
            if ref_no not in ref_dict:
                ref_dict[ref_no] = []
            if seg_no not in utg_scores:
                utg_scores[seg_no] = match_rate
            else:
                utg_scores[seg_no] = max(match_rate, utg_scores[seg_no])
            if seg_no not in utg_mapq:
                utg_mapq[seg_no] = mapq
            else:
                utg_mapq[seg_no] = max(mapq, utg_mapq[seg_no])

            # less than 10% mismatch
            if match_rate >= 0.95:
                ref_dict[ref_no].append(seg_no)
            if seg_no == "utig4-876" and ref_no == "CM033309.1":
                print(splited)
                print(match_rate)
                print(match_rate >= 0.95)
                sys.exit(0)

    for ref_id, nodes in ref_dict.items():
        print(list_to_string(nodes, ref_id))

    print(">>>match rate distributions")
    regions, bins = numpy.histogram(a=list(utg_scores.values()))
    print(regions)
    print(bins)

    print(">>>mapq distributions")
    regions, bins = numpy.histogram(a=list(utg_mapq.values()))
    print(regions)
    print(bins)
    sys.exit(0)