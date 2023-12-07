#!/usr/bin/env python
# Created: 24 Nov 2023
# Author: Gomathinayagam, gomathinayagam.s@vit.ac.in

import sys
import argparse
from Bio import SeqIO
import re

def stop_codons(code):
    stops = {
        1: "TAA|TAG|TGA",
        2: "TAA|TAG|AGA|AGG",
        3: "TAA|TAG",
        # ... (remaining entries)
    }

    if code not in stops:
        raise ValueError(f"Not a valid genetic code: {code}")

    return stops[code]

def cigar_to_frameshift(cigar, qstart, qend, qlen):
    cig_iter = re.finditer(r'(\d+)([MID])', cigar)
    fpos = []
    fsym = []

    rev = qstart > qend
    pos = qlen - qstart + 1 if rev else qstart

    for match in cig_iter:
        num, sym = match.groups()
        num = int(num)

        if sym == 'M' or sym == 'I':
            pos += num * 3

        elif sym == '\\':
            pos += 1
            if cig_iter.__next__().group(2) == "D":
                fpos.append(pos - 2)
                fsym.append(2)
            else:
                fpos.append(pos - 2)
                fsym.append(-1)

        elif sym == '/':
            pos -= 1
            if cig_iter.__next__().group(2) == "I":
                fpos.append(pos - 3)
                fsym.append(-2)
            else:
                fpos.append(pos - 2)
                fsym.append(1)

    if rev:
        fpos = [qlen - x + 1 for x in reversed(fpos)]
        fsym.reverse()

    return fpos, fsym

def mask_internal_stops(seq, qmin, qmax, rev):
    m = 5  # n codons from start&end ignored

    ss = seq[qmin-1:qmax]
    if len(ss) < m * 7:
        return seq

    if rev:
        ss = ss.reverse_complement()

    codons = [ss[i:i+3] for i in range(0, len(ss), 3)]
    for i in range(m, len(codons) - m):
        if re.match(stop_codons(genetic_code), codons[i], re.IGNORECASE):
            codons[i] = 'nnn'
            # Increment the count of masked internal stops
            global total_stop_masked
            total_stop_masked += 1

    seq[qmin-1:qmax] = ''.join(codons)

    if rev:
        seq[qmin-1:qmax] = seq[qmin-1:qmax].reverse_complement()

    return seq

def get_alignments(id):
    global next_id
    alns = []

    if next_id is not None and id != next_id:
        return []

    while True:
        line = mapps.readline().rstrip('\n')
        if not line:
            break

        qid, sid, pid, length, mismatches, gap, qstart, qend, sstart, send, evalue, bitscore, qlen, cigar = line.split('\t')
        if 'I' in cigar or 'D' in cigar:
            if qid == id:
                alns.append([cigar, int(qstart), int(qend), int(qlen), sid])
            else:
                alns = sorted(alns, key=lambda x: max(x[1], x[2]), reverse=True)
                next_id = qid
                return alns

    if not line and (next_id is None or id == next_id):
        alns = sorted(alns, key=lambda x: max(x[1], x[2]), reverse=True)

    return alns

def main():
    global total_stop_masked
    total_aln_len = 0
    total_read_len = 0
    n_fix = 0
    n_reads = 0

    print("Processing seqs ...", file=sys.stderr)

    for seq in SeqIO.parse(open(sys.argv[1]), "fasta"):
        n_reads += 1
        total_read_len += len(seq.seq)
        print(f"\r{n_reads}", end='', file=sys.stderr, flush=True)

        alns = get_alignments(seq.id)
        fpos_prev = 1_000_000_000_000

        for aln in alns:
            cigar, qstart, qend, qlen, sid = aln
            qmin, qmax = sorted([qstart, qend])
            rev = qstart > qend

            total_aln_len += abs(qend - qstart) + 1
            fpos, fsym = cigar_to_frameshift(cigar, qstart, qend, qlen)

            if fpos:
                for j in range(len(fpos) - 1, -1, -1):
                    if fpos[j] < fpos_prev:
                        n_fix += 1
                        if fsym[j] > 0:  # del
                            seq.seq = seq.seq[:fpos[j] + 1] + "n" * fsym[j] + seq.seq[fpos[j] + 1:]
                        else:  # ins
                            seq.seq = seq.seq[:fpos[j]] + seq.seq[fpos[j] - fsym[j]:]
                        fpos_prev = fpos[j]
                        qmax += fsym[j]

            if not no_stop_masking:
                seq.seq = mask_internal_stops(seq.seq, qmin, qmax, rev)

        desc = seq.description
        desc = re.sub(r'SUBSTR:.*', '', desc)
        seq.description = desc

        print(seq.seq)

    print("", file=sys.stderr)

    if n_reads:
        print(f"Total processed seqs: {n_reads}", file=sys.stderr)
        print(f"Average seq coverage: {total_aln_len / total_read_len * 100:.1f}%", file=sys.stderr)
        print(f"Sites modified: {n_fix}", file=sys.stderr)
        print(f"Internal stops masked: {total_stop_masked}", file=sys.stderr)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="proovframe fix")
    parser.add_argument("input_fasta", help="Input FASTA file")
    parser.add_argument("input_tsv", help="Input TSV file")
    parser.add_argument("output_file", nargs='?', type=argparse.FileType('w'), default=sys.stdout, help="Output file (default: stdout)")
    parser.add_argument("-g", "--genetic-code", type=int, default=11, help="Genetic code table, sets stop codons (default: 11)")
    parser.add_argument("-S", "--no-stop-masking", action="store_true", help="Disable internal stop codon masking")
    parser.add_argument("-D", "--debug", action="store_true
