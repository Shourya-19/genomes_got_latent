# Written by Shourya Burnwal shouryaburnwal19@gmail.com

'''
Copyright (C) 2026 Shourya Burnwal <shouryaburnwal19@gmail.com>

You can use the contents of this file under the conditions of Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International deed (https://creativecommons.org/licenses/by-nc-sa/4.0/). 
A copy of this deed in provided in the LICENSE.md file.

You are free to copy, distribute, remix, transform and build upon this program as long as 1) it is for non commercial purpose, 2) Attributions are provided 3) Shared under the same license. 

'''
"""
BLAST Results Enricher
======================
Reads a BLAST tabular output file (-outfmt 6), a translated FASTA file,
and a UniRef proteome FASTA file, then produces a CSV with the original
BLAST columns plus:
  - sequence  : full amino acid sequence of the query ORF
  - entropy   : Shannon entropy (–Σ p·log2(p)) over amino acid frequencies
  - acc_id    : UniRef accession (e.g. Q63036)
  - function  : protein description from the UniRef header
  - genus     : genus name parsed from Tax= field
  - species   : species epithet parsed from Tax= field

Rows matching the user-supplied genus are excluded from the output.

Usage:
    python blast_to_csv.py <blast.tsv> <translations.fasta> <uniref.fasta> <output.csv> [exclude_genus]
"""

import sys
import csv
import math
from collections import Counter


# ── Constants ────────────────────────────────────────────────────────────────

BLAST_COLUMNS = [
    "qseqid", "sseqid", "pident", "length",
    "qstart", "qend", "sstart", "send",
    "evalue", "bitscore",
]


# ── Helpers ──────────────────────────────────────────────────────────────────

def read_fasta(filepath):
    """
    Parse a FASTA file and return a dict mapping header → sequence.
    The full header line (minus '>') is used as the key.
    """
    sequences = {}
    header = None
    parts = []

    with open(filepath, 'r') as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            if line.startswith('>'):
                if header is not None:
                    sequences[header] = ''.join(parts)
                header = line[1:]
                parts = []
            else:
                parts.append(line)

        if header is not None:
            sequences[header] = ''.join(parts)

    return sequences


def shannon_entropy(sequence):
    """
    Calculate Shannon entropy of an amino acid sequence.
    H = –Σ p(aa) · log2(p(aa))
    Returns 0.0 for sequences of length < 2.
    """
    seq = sequence.replace('*', '').replace('X', '')   # ignore stop / unknown
    n = len(seq)
    if n < 2:
        return 0.0

    counts = Counter(seq)
    entropy = 0.0
    for count in counts.values():
        p = count / n
        entropy -= p * math.log2(p)

    return round(entropy, 4)


def read_uniref_headers(filepath):
    """
    Scan a UniRef FASTA file and return a dict mapping the short ID
    (e.g. 'UniRef50_Q63036') to the full header line (without '>').

    Only headers are read — sequences are skipped — so this is fast
    even for very large proteome files.
    """
    id_to_header = {}
    print(f"Reading UniRef headers: {filepath}")

    with open(filepath, 'r') as fh:
        for line in fh:
            if not line.startswith('>'):
                continue
            header = line[1:].strip()
            short_id = header.split()[0]   # e.g. 'UniRef50_Q63036'
            id_to_header[short_id] = header

    print(f"  Loaded {len(id_to_header)} UniRef entries")
    return id_to_header



def parse_sseqid(sseqid):
    """
    Parse a UniRef sseqid string into its components.

    Example input:
        UniRef50_Q63036 Rat albumin (Fragment) n=1 Tax=Rattus norvegicus TaxID=10116 RepID=Q63036_RAT

    Returns a dict with keys: acc_id, function, genus, species.
    Fields that cannot be parsed are returned as empty strings.
    """
    result = {"acc_id": "", "function": "", "genus": "", "species": ""}

    # acc_id – everything after the last underscore in the first token
    # e.g. UniRef50_Q63036 → Q63036
    first_token = sseqid.split()[0]             # 'UniRef50_Q63036'
    result["acc_id"] = first_token.rsplit('_', 1)[-1]

    rest = sseqid[len(first_token):].strip()    # 'Rat albumin (Fragment) n=1 Tax=...'

    # function – text before the first metadata tag
    func_end = len(rest)
    for marker in (' n=', ' Tax=', ' TaxID=', ' RepID='):
        idx = rest.find(marker)
        if idx != -1 and idx < func_end:
            func_end = idx
    result["function"] = rest[:func_end].strip()

    # genus + species – from 'Tax=Rattus norvegicus'
    tax_idx = rest.find('Tax=')
    if tax_idx != -1:
        tax_str = rest[tax_idx + len('Tax='):]
        # stop at the next tag (TaxID=, RepID=, etc.)
        for stop in (' TaxID=', ' RepID=', ' n='):
            stop_idx = tax_str.find(stop)
            if stop_idx != -1:
                tax_str = tax_str[:stop_idx]
        tax_parts = tax_str.strip().split()
        if len(tax_parts) >= 1:
            result["genus"] = tax_parts[0]
        if len(tax_parts) >= 2:
            result["species"] = tax_parts[1]

    return result


def match_query_to_sequence(qseqid, fasta_dict):
    """
    Find the FASTA sequence whose header matches qseqid.
    The query ID uses commas as separators (e.g. 'test_chr,F6,st...,sp...,L...')
    which corresponds to the FASTA header written by dna_translator.py:
      'test_chr | F6 | st... sp... L...'
    We try three strategies in order:
      1. Exact match on the full header.
      2. Match after normalising separators (comma ↔ pipe/space).
      3. Match on the first token (chromosome / sequence name) alone.
    """
    # Strategy 1 – exact
    if qseqid in fasta_dict:
        return fasta_dict[qseqid]

    # Strategy 2 – normalise separators
    # BLAST query ID has commas; FASTA headers use ' | ' as separator
    normalised = qseqid.replace(',', ' | ')
    for header, seq in fasta_dict.items():
        if header == normalised:
            return seq
        # also try stripping trailing metadata like ' (no stop codon)'
        if header.split(' (')[0].strip() == normalised.strip():
            return seq

    # Strategy 3 – match on all comma-separated tokens present anywhere in header
    tokens = [t.strip() for t in qseqid.split(',')]
    for header, seq in fasta_dict.items():
        if all(tok in header for tok in tokens):
            return seq

    return None   # not found


def ask_exclude_genus(args):
    """
    Return a set of genus names to exclude (lowercase for case-insensitive comparison).
    Accepts a comma-separated list as optional 6th+ CLI argument or prompts interactively.
    Returns an empty set if the user skips.
    """
    if len(args) >= 6:
        raw = args[5].strip()
    else:
        raw = input(
            "Enter genus name(s) to exclude (comma-separated, e.g. Rattus,Homo) "
            "[press Enter to keep all genera]: "
        ).strip()

    if not raw:
        return set()

    genera = {g.strip().lower() for g in raw.split(',') if g.strip()}
    return genera


# ── Main ─────────────────────────────────────────────────────────────────────

def main():
    if len(sys.argv) < 5:
        print("Usage: python blast_to_csv.py <blast.tsv> <translations.fasta> <uniref.fasta> <output.csv> [exclude_genus]")
        sys.exit(1)

    blast_file        = sys.argv[1]
    translations_file = sys.argv[2]
    uniref_file       = sys.argv[3]
    output_file       = sys.argv[4]

    exclude_genus = ask_exclude_genus(sys.argv)
    if exclude_genus:
        print(f"  Excluding genera           : {', '.join(sorted(exclude_genus))}")
    else:
        print("  No genus exclusion applied.")

    # ── Load translation sequences ──
    print(f"Reading translations: {translations_file}")
    fasta_dict = read_fasta(translations_file)
    print(f"  Loaded {len(fasta_dict)} sequences")

    # ── Load UniRef headers ──
    uniref_headers = read_uniref_headers(uniref_file)

    MIN_ENTROPY = 3.4

    # ── Parse BLAST file into rows ──
    print(f"Reading BLAST results: {blast_file}")

    raw_rows = []
    malformed = 0

    with open(blast_file, 'r') as blast_fh:
        for line in blast_fh:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            fields = line.split('\t')
            if len(fields) != len(BLAST_COLUMNS):
                print(f"  Warning: skipping malformed line ({len(fields)} fields): {line[:80]}")
                malformed += 1
                continue
            raw_rows.append(dict(zip(BLAST_COLUMNS, fields)))

    print(f"  Total BLAST hits read      : {len(raw_rows)}")

    # ── Keep only best e-value hit per qseqid ──
    best_hits = {}
    for row in raw_rows:
        qid = row["qseqid"]
        try:
            evalue = float(row["evalue"])
        except ValueError:
            evalue = float('inf')

        if qid not in best_hits or evalue < best_hits[qid][1]:
            best_hits[qid] = (row, evalue)

    best_rows = [row for row, _ in best_hits.values()]
    print(f"  Unique query IDs (best hit): {len(best_rows)}")

    # ── Attach sequence + entropy + parsed sseqid fields; apply filters ──
    OUTPUT_COLUMNS = BLAST_COLUMNS + ["acc_id", "function", "genus", "species", "sequence", "entropy"]

    rows_written   = 0
    rows_no_seq    = 0
    rows_low_ent   = 0
    rows_genus_hit = 0
    missing_ids    = set()

    with open(output_file, 'w', newline='') as out_fh:
        writer = csv.DictWriter(out_fh, fieldnames=OUTPUT_COLUMNS)
        writer.writeheader()

        for row in best_rows:
            qseqid = row["qseqid"]

            # ── Look up full UniRef header for this sseqid ──
            short_id = row["sseqid"].strip()
            full_header = uniref_headers.get(short_id, short_id)  # fallback to raw id if not found

            # ── Parse sseqid ──
            parsed = parse_sseqid(full_header)
            row["acc_id"]   = parsed["acc_id"]
            row["function"] = parsed["function"]
            row["genus"]    = parsed["genus"]
            row["species"]  = parsed["species"]

            # ── Genus filter ──
            if exclude_genus and parsed["genus"].lower() in exclude_genus:
                rows_genus_hit += 1
                continue

            # ── Sequence lookup ──
            seq = match_query_to_sequence(qseqid, fasta_dict)

            if seq:
                h = shannon_entropy(seq)
                if h < MIN_ENTROPY:
                    rows_low_ent += 1
                    continue
                row["sequence"] = seq
                row["entropy"]  = h
                rows_written += 1
            else:
                row["sequence"] = "NOT_FOUND"
                row["entropy"]  = "N/A"
                rows_no_seq += 1
                missing_ids.add(qseqid)

            writer.writerow(row)

    # ── Summary ──
    print(f"  Removed (genus filter)     : {rows_genus_hit}")
    print(f"  Removed (low entropy < {MIN_ENTROPY}) : {rows_low_ent}")
    print(f"  Sequences not found        : {rows_no_seq}")
    print(f"  Rows written to CSV        : {rows_written + rows_no_seq}")
    if missing_ids:
        print(f"\n  Query IDs with no matching FASTA sequence:")
        for mid in sorted(missing_ids):
            print(f"    {mid}")
    print(f"\nDone. CSV written to: {output_file}")


if __name__ == "__main__":
    main()
