import sys, os
import pysam

def System(command):
    return os.system(command)

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


def build_mat(bam_file: str, utgs: list):
    utg2idx = {}
    for idx, utg in enumerate(utgs):
        utg2idx[utg] = idx
    
    adj_list = dict()
    print("getting bam reader...")
    bam_reader = pysam.AlignmentFile(bam_file, 'rb')  # BAM file reader.
    # Iterate through reads.
    read1 = None
    read2 = None
    num_pair = 0
    total_iter = 0
    print("Start building matrix")
    for read in bam_reader:
        total_iter += 1

        if total_iter % 1000000 == 0:
            print(f"bam file iterated up to {total_iter}")

        if not read.is_paired or read.mate_is_unmapped or read.is_duplicate:
            continue
        if read.is_read2:
            read2 = read
            assert read1 != None and read2 != None and read1.query_name == read2.query_name
            num_pair += 1
            
            idx1 = utg2idx[read.reference_name]
            idx2 = utg2idx[read2.reference_name]
            min_pair = (min(idx1, idx2), max(idx1, idx2))
            if min_pair not in adj_list:
                adj_list[min_pair] = 0
            adj_list[min_pair] += 1

        else:
            read1 = read
            read2 = None
    print("num pairs: ", num_pair)
    print("total iter: ", total_iter)

    System("echo "" > adj_list.txt")
    with open("adj_list.txt", "w") as adj_fd:
        for (idx1, idx2), val in adj_list.items():
            adj_fd.write(f"{utgs[idx1]}\t{utgs[idx2]}\t{val}\n")
        adj_fd.close()
    return adj_list

def mcl_clustering(mat: list, utgs: list):
    # params
    inflation_rate=2.0
    num_threads=16

    print("Convert graph to abc format")

    file_prefix="mcl_graph"
    if not os.path.exists(file_prefix + ".abc"):
        os.system("echo "" > {0}.abc".format(file_prefix))
        with open(file_prefix + ".abc", "w") as fd:
            for i, vs in enumerate(mat):
                for j, v in enumerate(vs):
                    fd.write(f"{utgs[i]}\t{utgs[j]}\t{v}\n")
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

# attempt to link duplex reads using hic
# goal use hic reads to partition contigs
if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("{0} <sorted_bam> <duplex_read>".format(sys.argv[0]))
        sys.exit(1)
    
    _, sorted_bam, duplex_read = sys.argv
    if not os.path.exists("reads.txt"):
        reads = get_reads_fastq(duplex_read)
        store_reads(reads, True)

    read_ids = get_reads(True)
    if not os.path.exists("adj_list.txt"):
        mat = build_mat(sorted_bam, read_ids)

    same_counter = 0
    with open("adj_list.txt", "r") as adj_fd:
        for line in adj_fd:
            uid, vid, links = line.strip().split()
            if uid != vid:
                print(uid, vid)
            else:
                same_counter += 1
        adj_fd.close()
    
    print(same_counter)


    # mcl_clustering(mat, utg_ids)



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
