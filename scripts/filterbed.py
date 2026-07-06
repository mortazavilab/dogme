import gzip
import sys


def open_maybe_gzip(path, mode):
    if path.endswith('.gz'):
        return gzip.open(path, mode)
    return open(path, mode)


minCov = int(sys.argv[1])
perMod = int(sys.argv[2])
infilename = sys.argv[3]
outfilename = sys.argv[4]

with open_maybe_gzip(infilename, 'rt') as infile, open_maybe_gzip(outfilename, 'wt') as outfile:
    incount = 0
    outcount = 0
    for line in infile:
        incount += 1
        fields = line.split('\t')
        fields[4] = fields[10]
        if (int(fields[11]) > minCov) & (float(fields[10]) > perMod):
            outcount += 1
            outfile.write("\t".join(fields))

print('input read ' + str(incount))
print('after filtering ' + str(outcount))
