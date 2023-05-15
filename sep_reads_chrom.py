# separate reads using minimap2 and chromosome reference directly.
import sys
import os
import pysam


def get_reads(reads_file, only_id=True):
    ids = []
    with open(reads_file, "r") as fd:
        for line in fd:
            if only_id:
                ids.append(line.strip())
            else:
                raise Exception
        fd.close()
    return ids

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print(f"{sys.argv[0]} <bam> <read_info>")
        sys.exit(0)
    
    _, bam_file, read_info = sys.argv

    ref_dict = {}
    ref_dict["unknown"] = []
    read_id_set = set()
    bam_reader = pysam.AlignmentFile(bam_file, 'rb')  # BAM file reader.
    # Iterate through reads.
    total_iter = 0
    for read in bam_reader:
        total_iter += 1

        if total_iter % 1000000 == 0:
            print(f"bam file iterated up to {total_iter}")

        if read.mate_is_unmapped or read.is_duplicate:
            continue

        # valid read mapping, map the read to best aligned reference
        if read.query_name in read_id_set:
            print("duplicated read id found, error: ", read.query_name)
            sys.exit(1)

        ref_id = str(read.reference_name)
        if not str.startswith(ref_id, "chr"):
            ref_dict["unknown"].append(read.query_name)
        else:
            if ref_id not in ref_dict:
                ref_dict[ref_id] = []
            ref_dict[ref_id].append(read.query_name)

    print("total iter: ", total_iter)

    read_ids = get_reads(read_info)
    read_dict = dict.fromkeys(read_ids, "")


    for ref_id, read_ids in ref_dict.items():
        with open(f"gtruth_duplex_{ref_id}.txt", "w") as fd:
            for rid in read_ids:
                fd.write(f"{rid}\n")
            fd.close()
        for read_id in read_ids:
            read_dict[read_id] = ref_id
        print(f"reference {ref_id} duplex read processed: ", len(read_ids))
    
    with open("reads_gtruth.txt", "w") as rfd:
        for read_id, ref_id in read_dict.items():
            rfd.write(f"{read_id}\t{ref_id}\n")
        rfd.close()
