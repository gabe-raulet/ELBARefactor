#!/usr/bin/env python3

import sys
import getopt

base_for = "ACGT"
base_rev = "TGCA"
comp_tab = str.maketrans(base_for, base_rev)

kmer_size = 31

def usage():
    global kmer_size
    sys.stderr.write("Usage: {} [options] <reads.fa>\n".format(sys.argv[0]))
    sys.stderr.write("    -k INT   k-mer size [{}]\n".format(kmer_size))
    sys.stderr.write("    -H FILE  histogram file\n")
    sys.stderr.write("    -h       help message\n")
    sys.exit(1)

def read_fasta(fasta_fname):
    names = []
    reads = []
    seq = []
    name = ""
    for line in open(fasta_fname, "r"):
        if line.startswith(">"):
            if len(seq) > 0:
                reads.append("".join(seq))
                names.append(name)
                seq = []
            name = line.lstrip(">").split()[0].rstrip()
        else:
            seq.append(line.rstrip())
    if len(seq) > 0:
        reads.append("".join(seq))
        names.append(name)
    return reads, names

def main(argc, argv):
    global kmer_size
    histofname = None
    try: opts, args = getopt.gnu_getopt(argv[1:], "k:H:h")
    except getopt.GetoptError as err:
        sys.stderr.write("error: {}\n".format(err))
        sys.exit(1)
    for o, a in opts:
        if o == "-k": kmer_size = int(a)
        elif o == "-H": histofname = a
        elif o == "-h": usage()
    if len(args) != 1:
        sys.stderr.write("error: missing <reads.fa>\n")
        sys.exit(1)
    reads, names = read_fasta(args[0])
    kmertable = {}
    for seq in reads:
        l = len(seq)
        if l < kmer_size: continue
        for i in range(l-kmer_size+1):
            kmer_for = seq[i:i+kmer_size]
            kmer_rev = kmer_for.translate(comp_tab)[::-1]
            kmer = min(kmer_for, kmer_rev)
            if kmer in kmertable: kmertable[kmer] += 1
            else: kmertable[kmer] = 1
    if not histofname is None:
        histo = [0] * max(kmertable.values())
        for kmer, cnt in kmertable.items():
            histo[cnt-1] += 1
        with open(histofname, "w") as f:
            f.write("#count\tnumkmers\n")
            for i in range(len(histo)):
                if histo[i] > 0:
                    f.write("{}\t{}\n".format(i+1, histo[i]))
    distinctmers = len(kmertable)
    totalmers = sum(kmertable.values())
    sys.stdout.write("#distinct {}-mers: {}\n".format(kmer_size, distinctmers))
    sys.stdout.write("average #{}-mers per distinct {}-mer: {:.3f}\n".format(kmer_size, kmer_size, totalmers/distinctmers))

if __name__ == "__main__":
    sys.exit(main(len(sys.argv), sys.argv))
