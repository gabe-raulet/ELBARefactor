import sys

base_for = "ACGT"
base_rev = "TGCA"
comp_tab = str.maketrans(base_for, base_rev)

def read_fasta(fname):
    seq = []
    seqs = []
    for line in open(fname, "r"):
        if line[0] == ">":
            if len(seq) > 0:
                seqs.append("".join(seq))
                seq = []
        else: seq.append(line.rstrip())
    if len(seq) > 0:
        seqs.append("".join(seq))
    return seqs



if __name__ == "__main__":
    if len(sys.argv) != 5:
        sys.stderr.write("usage: {} [reads.fa] [idx] [start] [len]\n".format(sys.argv[0]))
        sys.exit(-1)
    seqs = read_fasta(sys.argv[1])
    idx = int(sys.argv[2])
    start = int(sys.argv[3])
    length = int(sys.argv[4])
    seq = seqs[idx]
    mer = seq[start:start+length]
    revmer = mer.translate(comp_tab)[::-1]
    sys.stdout.write("{}\n{}\n".format(mer, revmer))
