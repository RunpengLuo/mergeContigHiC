import sys, os
import pysam

from utils import *


def alignment_bwa(ref_file: str, fwd_file: str, rve_file: str):
    print("bwa mem alignment")
    mapped_fwd_bam = "forward.bam"
    mapped_rve_bam = "reverse.bam"

    if not os.path.exists(mapped_fwd_bam) or not os.path.exists(mapped_rve_bam):
        System("bwa index {0}".format(ref_file))
        System("bwa mem -A1 -B4 -E50 -L0 -t 16 {0} {1} 2>>fwd.log | samtools view -Shb -F 8 > {2}".format(ref_file, fwd_file, mapped_fwd_bam))
        System("bwa mem -A1 -B4 -E50 -L0 -t 16 {0} {1} 2>>rve.log | samtools view -Shb -F 8 > {2}".format(ref_file, rve_file, mapped_rve_bam))
        # "bwa mem -5SP -t8 {0} {1} {2} | samtools view -S -h -b -F 2316 > {3}"
    else:
        print("alignment exists")
    print("done")
    return mapped_fwd_bam, mapped_rve_bam

def build_mat_bam(fwd_bam: str, rve_bam: str, utgs: list, mapq: int):
    utg2idx = {}
    for idx, utg in enumerate(utgs):
        utg2idx[utg] = idx
    utg_idx_file = "utg.idx"
    if not os.path.exists(utg_idx_file):
        System("echo "" > " + utg_idx_file + "\n")
        with open(utg_idx_file, "w") as idx_fd:
            for idx, utg in enumerate(utgs):
                idx_fd.write(f"{idx}\t{utg}")
            idx_fd.close()

    hic2utg_dict = {}
    count_diff = 0
    count_pair = 0
    bam_fwd_reader = pysam.AlignmentFile(fwd_bam, 'rb')  # BAM file reader.
    for read in bam_fwd_reader:
        if read.mapping_quality >= mapq:
            hic2utg_dict[read.query_name] = Pair.init_pair_fwd(read.reference_name)

    bam_rve_reader = pysam.AlignmentFile(rve_bam, 'rb')  # BAM file reader.
    for read in bam_rve_reader:
        if read.query_name in hic2utg_dict and read.mapping_quality >= mapq:
            hic2utg_dict[read.query_name].assign(read.reference_name, False)
            if hic2utg_dict[read.query_name].F != hic2utg_dict[read.query_name].R:
                count_diff += 1
            count_pair += 1

    print("Total aligned pair: ", count_pair)
    print("Total rdiff pair: ", count_diff)

    return c2s_mat(hic2utg_dict, utgs, utg2idx, mapq)


def alignment_minimap2(ref_file: str, fwd_file: str, rve_file: str):
    print("minimap2 alignment")
    mapped_fwd_paf = "mapped.forward.paf"
    mapped_rve_paf = "mapped.reverse.paf"
    if not os.path.exists(mapped_fwd_paf) or not os.path.exists(mapped_rve_paf):
        System("minimap2 --secondary=no {0} {1} > {2}".format(ref_file, fwd_file, mapped_fwd_paf))
        System("minimap2 --secondary=no {0} {1} > {2}".format(ref_file, rve_file, mapped_rve_paf))
    else:
    #    print("sorted mapped bam exists")
        print("alignment exists")
    print("done")
    return mapped_fwd_paf, mapped_rve_paf

# between contig alignment
def build_mat_minimap2(fwd_paf: str, rve_paf: str, utgs: list, qscore):
    utg2idx = {}
    for idx, utg in enumerate(utgs):
        utg2idx[utg] = idx
    
    hic2utg_dict = {}
    fwd_fd = open(fwd_paf, "r")
    rve_fd = open(rve_paf, "r")

    for Line in fwd_fd:
        splited = Line.strip().split("\t")
        seg_no = str(splited[0])
        ref_no = str(splited[5])
        hic2utg_dict[seg_no] = Pair.init_pair_fwd(ref_no)

    count_diff = 0
    count_pair = 0
    for Line in rve_fd:
        splited = Line.strip().split("\t")
        seg_no = str(splited[0])
        ref_no = str(splited[5])
        if seg_no in hic2utg_dict:
            hic2utg_dict[seg_no].assign(ref_no, False)
            if hic2utg_dict[seg_no].F != hic2utg_dict[seg_no].R:
                count_diff += 1
            count_pair += 1
    
    print("Total aligned pair: ", count_pair)
    print("Total rdiff pair: ", count_diff)

    return c2s_mat(hic2utg_dict, utgs, utg2idx)

def mcl_clustering(mat: dict, utgs: list, inflation_rate: float, mapq: int, num_threads: int):
    print("run mcl clustering")
    print(f"inflation_rate: {inflation_rate}")
    print(f"mapq: {mapq}")
    # params
    print("Convert graph to abc format")

    file_prefix="mcl_graph_Q" + str(mapq) + "_I" + str(inflation_rate)
    if not os.path.exists(file_prefix + ".abc"):
        System("echo "" > {0}.abc".format(file_prefix))
        with open(file_prefix + ".abc", "w") as fd:
            for (idx1, idx2), val in mat.items():
                fd.write(f"{utgs[idx1]}\t{utgs[idx2]}\t{val}\n")
            fd.close()
        print("done")
    else:
        print("abc format file exists, skip")

    print("Convert abc format to binary mci & tab format")
    if not os.path.exists(file_prefix + ".mci"):
        System("mcxload -abc {0}.abc --stream-mirror --write-binary -o {0}.mci -write-tab {0}.tab".format(file_prefix))
        print("done")
    else:
        print("mci format file exists, skip")

    print("Solve mcl clusterig")
    if not os.path.exists(file_prefix + ".icl"):
        System("mcl {0}.mci -te {1} -I {2} -o {0}.icl".format(
            file_prefix, num_threads, inflation_rate))
        print("Done")
    else:
        print("icl file exists, skip")

    print("Convert to cluster result")
    if not os.path.exists(file_prefix + ".cluster.txt"):
        System("mcxdump -icl {0}.icl -o {0}.cluster.txt -tabr {0}.tab".format(file_prefix))
        print("Done")
    else:
        print("cluster result exists")
    return


# goal use hic reads to partition contigs
if __name__ == "__main__":
    if len(sys.argv) < 4 or len(sys.argv) > 6:
        print("{0} <contigs> <hic_forward> <hic_reverse> <mapq: 0> <inflation_rate: 2.0>".format(sys.argv[0]))
        sys.exit(1)

    _, contig_file, fwd_file, rve_file = sys.argv[:4]
    
    mapq = int(sys.argv[4]) if len(sys.argv) >= 5 else 0
    inflation_rate = float(sys.argv[5]) if len(sys.argv) >= 6 else 2.0
    print("getting hic-contig alignment")
    mapped_fwd, mapped_rve = alignment_bwa(contig_file, fwd_file, rve_file)

    utgs = get_utgs(contig_file)
    utg_ids = list(utgs.keys())
    mat, adj_list = build_mat_bam(mapped_fwd, mapped_rve, utg_ids, mapq)
    # mat, adj_list = build_mat_minimap2(mapped_fwd_paf, mapped_rve_paf, utg_ids)

    cluster_file = mcl_clustering(adj_list, utg_ids, inflation_rate, mapq, 16)
    sys.exit(0)

# bash mmp.sh references/chrAll.fna data/assembly.fasta
# # sort and index for fast access
# samtools sort mapped.bam > mapped.sorted.bam
# samtools index mapped.sorted.bam

# ## Extract chromsosome names:
# samtools idxstats mapped.sorted.bam | cut -f1 | grep -v '*' > chr.names

# ## Split bam file with w while loop
# while read p
#   do
#   samtools view -o mapped_${p}.bam mapped.sorted.bam ${p}
#   done < chr.names

# while read p
#   do
#   samtools sort mapped_${p}.bam | samtools fastq -@ 16 -1 hic_${p}.forward.fastq -2 hic_${p}.reverse.fastq -0 /dev/null -s hic_${p}.single.fastq
#   done < chr.names
