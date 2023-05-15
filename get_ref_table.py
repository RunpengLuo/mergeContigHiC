import os
import sys


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print(f"{sys.argv[0]} <ref.fasta> <out_table_file>")
        sys.exit(1)
    
    _, ref_file, table_file = sys.argv

    table = {}

    with open(ref_file, "r") as ref:
        name_in_file = ""
        content = ""
        for line in ref.readlines():
            if line.startswith(">"):
                if name_in_file != "":
                    # process current content and move to next header
                    sname = str(name_in_file[1:-1].split()[0])
                    table[sname] = content
                # clear previous content
                name_in_file = line
                content = ""
            else:
                content += line.strip()
        ref.close()
    
    sorted_table = sorted(table.items(), key=lambda tple: tple[0])

    with open(table_file, "w") as tfd:
        for (name, seq) in sorted_table:
            tfd.write(f"{name}\t{len(seq)}\n")
        tfd.close()
