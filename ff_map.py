#!/usr/bin/env python
# Created: 20 Nov 2023
# Author: Gomathinayagam, gomathinayagam.s@vit.ac.in

import sys
import argparse
import os
import subprocess
from os.path import basename
from textwrap import wrap

def create_database(aa, db, threads):
    if not os.path.exists(db):
        if aa is None:
            raise ValueError("--aa is required")

        dmnd_makedb = f"diamond makedb -p {threads} --in {aa} --db {db}"
        print(wrap("  ", "", dmnd_makedb), "\n\n")
        if not dry:
            subprocess.run(dmnd_makedb, shell=True, check=True)
    elif aa is not None:
        print(f"  {db} already exists, will reuse\n")

def map_proteins(reads, db, out, threads, mode):
    if not mode or mode == "fast":
        mode = ""
    else:
        mode = f"--{mode}"

    dmnd_blastx = f"diamond blastx --db {db} --query {reads} --out {out} " \
                  f"-p {threads} {mode} --range-culling -F 1 -k 1 " \
                  "--outfmt 6 qseqid sseqid pident length mismatch gapopen qstart " \
                  "qend sstart send evalue bitscore qlen cigar sseq btop {' '.join(sys.argv[1:])}"

    print(wrap("  ", "", dmnd_blastx), "\n\n")
    if not dry:
        subprocess.run(dmnd_blastx, shell=True, check=True)

def main():
    global dry

    parser = argparse.ArgumentParser(description="proovframe map")
    parser.add_argument("reads", help="Input sequence file (e.g., FASTA)")
    parser.add_argument("-t", "--threads", type=int, default=4, help="Number of CPU threads")
    parser.add_argument("-d", "--db", help="Database file path")
    parser.add_argument("-a", "--aa", help="Protein file, not required if --db provided")
    parser.add_argument("-o", "--out", help="Write alignments to this file")
    parser.add_argument("-m", "--diamond-mode", default="more-sensitive", help="Diamond mode (default: more-sensitive)")
    parser.add_argument("-y", "--dry-run", action="store_true", help="Print the diamond command, but don't run it")
    parser.add_argument("-h", "--help", action="store_true", help="Show this help")
    parser.add_argument("-D", "--debug", action="store_true", help="Show debug messages")

    args = parser.parse_args()

    if args.help or not args.reads:
        print("Usage: proovframe map [-a|-d proteins] -o alignments.o6 seqs.fa -- extra-diamond-params")
        print(f"{'%-19s  %s' % ('-a/--aa', 'Protein file, not required if --db provided')}")
        print(f"{'%-19s  %s' % ('-d/--db', 'Created if not existing and --aa given [basename(aa).dmnd]')}")
        print(f"{'%-19s  %s' % ('-o/--out', 'Write alignments to this file [basename(seqs).o6]')}")
        print(f"{'%-19s  %s' % ('-t/--threads', 'Number of CPU threads')}")
        print(f"{'%-19s  %s' % ('-m/--diamond-mode', 'One of fast, sensitive, {mid,more,very,ultra}-sensitive' ' [more-sensitive]')}")
        print(f"{'%-19s  %s' % ('-y/--dry-run', 'Print the diamond command, but don\'t run it')}")
        print(f"{'%-19s  %s' % ('-h/--help', 'Show this help')}")
        print(f"{'%-19s  %s' % ('-D/--debug', 'Show debug messages')}")
        print("\nFor consensus sequences with rather low expected error rates "
              "and if your reference database has a good representation of similar sequences, "
              "you might want to switch to '-m fast' or '-m sensitive' to speed things up.")
        print("Also note, I've experienced inefficient parallelization if "
              "correcting a small number of Mb-sized genomes (as opposed to thousands "
              "of long-reads) - presumably because diamond threads on a per-sequence basis.")
        exit(0)

    dry = args.dry_run
    create_database(args.aa, args.db, args.threads)
    map_proteins(args.reads, args.db, args.out, args.threads, args.diamond_mode)
    print("\nDone")

if __name__ == "__main__":
    main()
