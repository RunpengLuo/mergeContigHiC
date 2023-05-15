import sys, os
import pysam

from utils import *


def alignment_bwa(ref_file: str, fwd_file: str, rve_file: str):
    print("bwa mem alignment")
    mapped_fwd_bam = "mapped.forward.bam"
    mapped_rve_bam = "mapped.reverse.bam"

    if not os.path.exists(mapped_fwd_bam) or not os.path.exists(mapped_rve_bam):
        System("bwa index {0}".format(ref_file))
        System("bwa mem -A1 -B4 -E50 -L0 {0} 2>>fwd.log | samtools view -Shb -F 8 > {1}".format(fwd_file, mapped_fwd_bam))
        System("bwa mem -A1 -B4 -E50 -L0 {0} 2>>rve.log | samtools view -Shb -F 8 > {1}".format(rve_file, mapped_rve_bam))
        # "bwa mem -5SP -t8 {0} {1} {2} | samtools view -S -h -b -F 2316 > {3}"
    else:
        print("alignment exists")
    print("done")
    return mapped_fwd_bam, mapped_rve_bam

def build_mat_bam(fwd_bam: str, rve_bam: str, utgs: list):
    utg2idx = {}
    for idx, utg in enumerate(utgs):
        utg2idx[utg] = idx
    
    hic2utg_dict = {}
    count_diff = 0
    count_pair = 0
    bam_fwd_reader = pysam.AlignmentFile(fwd_bam, 'rb')  # BAM file reader.
    for read in bam_fwd_reader:
        hic2utg_dict[read.query_name] = Pair.init_pair_fwd(read.reference_name)

    bam_rve_reader = pysam.AlignmentFile(rve_bam, 'rb')  # BAM file reader.
    for read in bam_rve_reader:
        if read.query_name in hic2utg_dict:
            hic2utg_dict[read.query_name].assign(read.reference_name, False)
            if hic2utg_dict[read.query_name].F != hic2utg_dict[read.query_name].R:
                count_diff += 1
            count_pair += 1

    print("Total aligned pair: ", count_pair)
    print("Total rdiff pair: ", count_diff)

    return c2s_mat(hic2utg_dict, utg2idx)

    # Iterate through reads.
    # read1 = None
    # read2 = None
    # num_pair = 0
    # total_iter = 0
    # for read in bam_reader:

    # hic2utg_dict[seg_no] = Pair.init_pair_fwd(ref_no)
    return


# within contig alignment
def build_mat_bam_single(bam_file: str, utgs: list):
    utg2idx = {}
    for idx, utg in enumerate(utgs):
        utg2idx[utg] = idx
    
    mat = [[0 for _ in range(len(utgs))] for _ in range(len(utgs))]
    # mat[i][j] represents the hic linkage between ith and jth contig

    bam_reader = pysam.AlignmentFile(bam_file, 'rb')  # BAM file reader.
    # Iterate through reads.
    read1 = None
    read2 = None
    num_pair = 0
    total_iter = 0
    for read in bam_reader:
        total_iter += 1
        if not read.is_paired or read.mate_is_unmapped or read.is_duplicate:
            continue
        if read.is_read2:
            read2 = read
            assert read1 != None and read2 != None and read1.query_name == read2.query_name
            num_pair += 1
            
            mat[utg2idx[read.reference_name]][utg2idx[read2.reference_name]] += 1
            mat[utg2idx[read2.reference_name]][utg2idx[read.reference_name]] += 1
                
        else:
            read1 = read
            read2 = None
    print("num pairs: ", num_pair)
    print("total iter: ", total_iter)

    System("echo "" > mat.txt")
    with open("mat.txt", "w") as mat_fd:
        for arr in mat:
            mat_fd.write(",".join([str(i) for i in arr]) + "\n")
        mat_fd.close()
    return mat


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

    return c2s_mat(hic2utg_dict, utg2idx)

def c2s_mat(hic2utg_dict: dict, utg2idx: dict):
    # assign to adjacency matrix
    mat = None
    if len(utgs) < 500:
        mat = [[0 for _ in range(len(utgs))] for _ in range(len(utgs))]
        for pair in hic2utg_dict.values():
            if pair.F != None and pair.R != None:
                idx1 = utg2idx[pair.F]
                idx2 = utg2idx[pair.R]

                mat[idx1][idx2] += 1

        System("echo "" > mat.txt")
        with open("mat.txt", "w") as mat_fd:
            for arr in mat:
                mat_fd.write(",".join([str(i) for i in arr]) + "\n")
            mat_fd.close()
    else:
        print("matrix too large, skip..")

    adj_list = dict()
    for pair in hic2utg_dict.values():
        if pair.F != None and pair.R != None:
            idx1 = utg2idx[pair.F]
            idx2 = utg2idx[pair.R]

            min_pair = (min(idx1, idx2), max(idx1, idx2))
            if min_pair not in adj_list:
                adj_list[min_pair] = 0
            adj_list[min_pair] += 1

    System("echo "" > adj_list.txt")
    with open("adj_list.txt", "w") as adj_fd:
        for (idx1, idx2), val in adj_list.items():
            adj_fd.write(f"{utgs[idx1]}\t{utgs[idx2]}\t{val}\n")
        adj_fd.close()

    return mat, adj_list

def mcl_clustering(mat, utgs: list):
    # params
    inflation_rate=2.0
    num_threads=16

    print("Convert graph to abc format")

    file_prefix="mcl_graph"
    if not os.path.exists(file_prefix + ".abc"):
        os.system("echo "" > {0}.abc".format(file_prefix))
        with open(file_prefix + ".abc", "w") as fd:
            if isinstance(mat, list):
                for i, vs in enumerate(mat):
                    for j, v in enumerate(vs):
                        fd.write(f"{utgs[i]}\t{utgs[j]}\t{v}\n")
            else:
                for (idx1, idx2), val in mat.items():
                    fd.write(f"{utgs[idx1]}\t{utgs[idx2]}\t{val}\n")
            fd.close()
        print("done")
    else:
        print("abc format file exists, skip")

    print("Convert abc format to binary mci & tab format")
    if not os.path.exists(file_prefix + ".mci"):
        os.system("mcxload -abc {0}.abc --stream-mirror --write-binary -o {0}.mci -write-tab {0}.tab".format(file_prefix))
        print("done")
    else:
        print("mci format file exists, skip")

    print("Solve mcl clusterig")
    if not os.path.exists(file_prefix + ".icl"):
        os.system("mcl {0}.mci -te {1} -I {2} -o {0}.icl".format(
            file_prefix, num_threads, inflation_rate))
        print("Done")
    else:
        print("icl file exists, skip")

    print("Convert to cluster result")
    if not os.path.exists(file_prefix + ".cluster.txt"):
        os.system("mcxdump -icl {0}.icl -o {0}.cluster.txt -tabr {0}.tab".format(file_prefix))
        print("Done")
    else:
        print("cluster result exists")
    return


# goal use hic reads to partition contigs
if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("{0} <contigs> <hic_forward> <hic_reverse>".format(sys.argv[0]))
        sys.exit(1)
    
    _, contig_file, fwd_file, rve_file = sys.argv

    mapped_fwd_paf = "mapped.forward.paf"
    mapped_rve_paf = "mapped.reverse.paf"
    # mapped_bam = "mapped.bam"
    # sorted_bam = "mapped.sorted.bam"
    print("getting hic-contig alignment")
    mapped_fwd, mapped_rve = alignment_bwa(contig_file, fwd_file, rve_file)

    if contig_file.endswith("fasta"):
        utgs = get_edges(contig_file)
        utg_ids = list(utgs.keys())
    else:
        # fastq
        if not os.path.exists("reads.txt"):
            utgs = get_reads_fastq(contig_file)
            store_reads(utgs, True)
        utg_ids = get_reads(True)

    # mat = build_mat(sorted_bam, utg_ids)
    mat, adj_list = build_mat_bam(mapped_fwd, mapped_rve, utg_ids)
    # mat, adj_list = build_mat_minimap2(mapped_fwd_paf, mapped_rve_paf, utg_ids)
    if mat == None:
        mcl_clustering(adj_list, utg_ids)
    else:
        mcl_clustering(mat, utg_ids)
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
