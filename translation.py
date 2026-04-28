"""
DNA to Amino Acid Sequence Translator
Translates all 6 reading frames for one or more FASTA sequences.

Usage:
    python dna_translator.py <input.fasta> <output.fasta> [min_length]

Arguments:
    input.fasta   – FASTA file with one or more sequences
    output.fasta  – output file for translated ORFs
    min_length    – (optional) minimum protein length in amino acids to report
                    default: 60
"""

import sys
import time


# ── Constants ────────────────────────────────────────────────────────────────

DEFAULT_MIN_LENGTH = 60   # minimum ORF length (amino acids) if not specified by user

BASE_PAIR = {
    'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G',
    'U': 'A', 'N': 'N',
}

CODON_TABLE = {
    'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
    'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
    'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
    'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
    'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
    'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
    'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
    'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
    'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
    'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
    'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
    'TAC': 'Y', 'TAT': 'Y', 'TAA': '*', 'TAG': '*',
    'TGC': 'C', 'TGT': 'C', 'TGA': '*', 'TGG': 'W',
}


# ── I/O helpers ──────────────────────────────────────────────────────────────

def read_fasta(filepath):
    """Parse a FASTA file and return a list of (header, sequence) tuples."""
    records = []
    header = None
    seq_parts = []

    with open(filepath, 'r') as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if header is not None:
                    records.append((header, ''.join(seq_parts).upper()))
                header = line[1:]   # strip the '>'
                seq_parts = []
            else:
                seq_parts.append(line)

        # flush last record
        if header is not None:
            records.append((header, ''.join(seq_parts).upper()))

    return records


# ── Sequence operations ───────────────────────────────────────────────────────

def reverse_complement(sequence):
    """Return the reverse complement of a DNA sequence."""
    return ''.join(BASE_PAIR[base] for base in reversed(sequence))


# ── Translation ───────────────────────────────────────────────────────────────

def translate_frame(sequence, frame_offset, frame_number, seq_header, output_fh, min_length):
    """
    Translate one reading frame of *sequence* and write ORFs to *output_fh*.

    Parameters
    ----------
    sequence     : DNA string (forward or reverse-complement strand)
    frame_offset : 0, 1, or 2 – how many bases to skip at the start
    frame_number : 1–6  – used in the FASTA header of each ORF
    seq_header   : original FASTA header of the source sequence
    output_fh    : open file handle for writing results
    min_length   : minimum ORF length (aa) required to write the result
    """
    is_reverse = frame_number > 3
    seq_len = len(sequence)
    count = 0

    aa_seq = ""
    orf_start = 1 + frame_offset   # 1-based start on the *source* strand

    i = frame_offset
    while i + 3 <= seq_len:
        codon = sequence[i:i + 3]
        codon_pos = i + 1           # 1-based position in translated strand

        if codon in CODON_TABLE:
            amino_acid = CODON_TABLE[codon]

            if amino_acid == '*':
                # compute stop position on the original strand
                stop_pos = codon_pos if not is_reverse else seq_len - codon_pos + 1

                orf_len = len(aa_seq)
                if orf_len >= min_length:
                    header = (
                        f">{seq_header} | F{frame_number} | "
                        f"st {orf_start} sp {stop_pos} L {orf_len}"
                    )
                    output_fh.write(f"\n{header}\n{aa_seq}*\n")
                    count += 1

                aa_seq = ""
                next_codon_pos = codon_pos + 3
                orf_start = (
                    next_codon_pos if not is_reverse
                    else seq_len - next_codon_pos + 1
                )
            else:
                aa_seq += amino_acid
        else:
            # non-standard codon (e.g. contains N) → placeholder
            if len(codon) == 3:
                aa_seq += 'X'

        i += 3

    # write any trailing (open) ORF
    if aa_seq and len(aa_seq) >= min_length:
        stop_pos = i if not is_reverse else seq_len - i + 1
        orf_len = len(aa_seq)
        header = (
            f">{seq_header} | F{frame_number} | "
            f"st {orf_start} sp {stop_pos} L {orf_len} (no stop codon)"
        )
        output_fh.write(f"\n{header}\n{aa_seq}\n")
        count += 1

    return count


def translate_all_frames(seq_header, sequence, output_fh, min_length):
    """Translate all 6 reading frames for *sequence* and write to *output_fh*.
    Returns the total number of ORFs written."""
    rev_comp = reverse_complement(sequence)

    total = 0
    for frame_offset in range(3):
        total += translate_frame(sequence, frame_offset, frame_offset + 1, seq_header, output_fh, min_length)

    for frame_offset in range(3):
        total += translate_frame(rev_comp, frame_offset, frame_offset + 4, seq_header, output_fh, min_length)

    return total


# ── Main ──────────────────────────────────────────────────────────────────────

def parse_min_length(args):
    """
    Return the minimum ORF length to report.
    If a third argument is supplied, parse and validate it.
    Otherwise, ask the user interactively (default: DEFAULT_MIN_LENGTH).
    """
    if len(args) >= 4:
        try:
            value = int(args[3])
            if value < 0:
                raise ValueError
            return value
        except ValueError:
            print(f"  Warning: invalid min_length '{args[3]}', using default ({DEFAULT_MIN_LENGTH} aa).")
            return DEFAULT_MIN_LENGTH

    # interactive prompt
    raw = input(
        f"Minimum protein length to report (amino acids) "
        f"[press Enter for default = {DEFAULT_MIN_LENGTH}]: "
    ).strip()

    if raw == "":
        print(f"  Using default minimum length: {DEFAULT_MIN_LENGTH} aa")
        return DEFAULT_MIN_LENGTH

    try:
        value = int(raw)
        if value < 0:
            raise ValueError
        print(f"  Minimum length set to: {value} aa")
        return value
    except ValueError:
        print(f"  Invalid input '{raw}', using default ({DEFAULT_MIN_LENGTH} aa).")
        return DEFAULT_MIN_LENGTH


def main():
    if len(sys.argv) < 3:
        print("Usage: python dna_translator.py <input.fasta> <output.fasta> [min_length]")
        sys.exit(1)

    input_file  = sys.argv[1]
    output_file = sys.argv[2]

    min_length = parse_min_length(sys.argv)

    start_time = time.time()
    print(f"\nTranslating genome in all 6 frames: {input_file}")
    print(f"Reporting ORFs >= {min_length} amino acids\n")

    records = read_fasta(input_file)
    print(f"  Found {len(records)} sequence(s) in {input_file}")

    with open(output_file, 'w') as out_fh:
        # out_fh.write("# DNA → Protein translation – all 6 reading frames\n")
        # out_fh.write(f"# Input      : {input_file}\n")
        # out_fh.write(f"# Min length : {min_length} aa\n")

        grand_total = 0
        for idx, (header, sequence) in enumerate(records, start=1):
            print(f"  [{idx}/{len(records)}] {header[:60]}  ({len(sequence)} bp)")
            count = translate_all_frames(header, sequence, out_fh, min_length)
            print(f"    → {count} ORF(s) found (>= {min_length} aa)")
            grand_total += count

        print (f"\n{'='*70}\n# Grand total ORFs reported: {grand_total}\n{'='*70}\n")

    elapsed = time.time() - start_time
    print(f"\nTotal ORFs reported across all sequences: {grand_total}")
    print(f"Done. Results written to: {output_file}")
    print(f"Execution time: {elapsed:.3f} seconds")


if __name__ == "__main__":
    main()
