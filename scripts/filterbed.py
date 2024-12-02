import sys

minCov = int(sys.argv[1])
perMod = int(sys.argv[2])
infilename = sys.argv[3]
outfilename = sys.argv[4]

infile = open(infilename,'r')
outfile = open(outfilename,'w')

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

infile.close()
outfile.close()
