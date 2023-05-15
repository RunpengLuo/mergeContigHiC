import os
import sys

def System(command):
    return os.system(command)

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