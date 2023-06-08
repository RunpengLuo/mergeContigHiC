import os
import sys

def System(command):
    return os.system(command)

def get_utgs(file: str):
    if file.endswith("fasta") or file.endswith("fna"):
        return get_edges(file)
    else:
        return get_reads_fastq(file)

def get_edges(fasta_file):
    unitigs = {}
    prev_id = ""
    with open(fasta_file, "r") as fd:
        for l in fd:
            if l.startswith(">"):
                prev_id = l.strip().split()[0][1:]
            else:
                unitigs[prev_id] = l.strip()
                
        fd.close()
    return unitigs

def get_reads_fastq(fastq_file):
    reads = {}
    counter = 0
    with open(fastq_file, "r") as fd:
        prev_id = ""
        i = 0
        for l in fd:
            if i == 0:
                counter += 1
                # read id, get first non-space id
                prev_id = l.strip().split()[0][1:]
            elif i == 1:
                # read sequence
                reads[prev_id] = l.strip()
            i = (i+1) % 4
        fd.close()
        if prev_id != "":
            reads.append((prev_id, prev_seq))

    print("total reads: ", counter)
    return reads

def store_reads(reads: dict, only_id=True):
    # store basic read informations
    System("echo "" > reads.txt")
    with open("reads.txt", "w") as fd:
        for rid, rseq in reads.items():
            if only_id:
                fd.write(f"{rid}\n")
            else:
                fd.write(f"{rid}:{rseq}\n")
        fd.close()

def get_reads(only_id=True):
    ids = []
    with open("reads.txt", "r") as fd:
        for line in fd:
            if only_id:
                ids.append(line.strip())
            else:
                raise Exception
        fd.close()
    return ids

def c2s_mat(hic2utg_dict: dict, utgs: list, utg2idx: dict, mapq=0):
    # assign to adjacency matrix
    mat = None
    if len(utgs) < 500:
        print("storing matrix")
        mat = [[0 for _ in range(len(utgs))] for _ in range(len(utgs))]
        for pair in hic2utg_dict.values():
            if pair.F != None and pair.R != None:
                idx1 = utg2idx[pair.F]
                idx2 = utg2idx[pair.R]

                mat[idx1][idx2] += 1
        mat_f = f"mat_mapq_{mapq}.txt"
        System("echo "" > {0}".format(mat_f))
        with open(mat_f, "w") as mat_fd:
            for arr in mat:
                mat_fd.write(",".join([str(i) for i in arr]) + "\n")
            mat_fd.close()
        print("matrix is stored in : " + mat_f)
    else:
        print("matrix too large, skip..")

    print("storing adjacent list")
    adj_list = dict()
    for pair in hic2utg_dict.values():
        if pair.F != None and pair.R != None:
            idx1 = utg2idx[pair.F]
            idx2 = utg2idx[pair.R]

            min_pair = (min(idx1, idx2), max(idx1, idx2))
            if min_pair not in adj_list:
                adj_list[min_pair] = 0
            adj_list[min_pair] += 1
    adj_f = "adj_mapq_{0}.txt".format(mapq)
    System("echo "" > {0}".format(adj_f))
    with open(adj_f, "w") as adj_fd:
        for (idx1, idx2), val in adj_list.items():
            adj_fd.write(f"{utgs[idx1]}\t{utgs[idx2]}\t{val}\n")
        adj_fd.close()
    print("adj list is stored in : " + adj_f)
    return mat, adj_list

class Pair:
    def __init__(self) -> None:
        self.F = None
        self.R = None
    
    def assign(self, val, isF):
        if isF:
            self.F = val
        else:
            self.R = val
    
    def init_pair_fwd(val):
        pair = Pair()
        pair.F = val
        return pair



# bin
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