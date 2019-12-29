import re
import sys

in_db = sys.argv[1]
out_db = re.sub(".fasta$","_format.fasta",in_db)
print(str(out_db)+"\n")
of = open(out_db,"w")
f = open(in_db,"r")
for line in f:
    if line.startswith(">"):
        of.write(line.split(" ")[0]+"\n")
    else:
        of.write(line)
f.close()
of.close()
