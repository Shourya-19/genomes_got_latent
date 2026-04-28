# filter_no_hits.py
import sys
fasta_file = sys.argv[1] #translated seq file
blast_file = sys.argv[2] #blast output file: "blast_results.tsv"
output_file =sys.argv[3] # output name: "no_hits.fasta"


# Step 1: collect query IDs with hits
hit_ids = set()

with open(blast_file) as f:
    for line in f:
        if line.strip():
            qseqid = line.split("\t")[0]
            hit_ids.add(qseqid)


# Step 2: filter FASTA
total = 0
kept = 0

with open(fasta_file) as fin, open(output_file, "w") as fout:
    header = None
    seq_lines = []

    for line in fin:
        line = line.rstrip()

        if line.startswith(">"):
            if header is not None:
                total += 1
                seq_id = header[1:]

                if seq_id not in hit_ids:
                    fout.write(header + "\n")
                    fout.write("\n".join(seq_lines) + "\n")
                    kept += 1

            header = line
            seq_lines = []
        else:
            seq_lines.append(line)

    # last record
    if header is not None:
        total += 1
        seq_id = header[1:]

        if seq_id not in hit_ids:
            fout.write(header + "\n")
            fout.write("\n".join(seq_lines) + "\n")
            kept += 1


print("Total sequences:", total)
print("No-hit sequences:", kept)
print("Sequences with hits:", total - kept)
