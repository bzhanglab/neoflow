import pandas as pd
import argparse
import sys
import os
import time
import datetime
import re
import subprocess


def do_annovar(vcf:str, v_type:str, sample_name:str, out_dir:str, paramap:dict):

	out_name = out_dir + "/" + os.path.basename(vcf)
	vcf_format = False

	cmd = "perl " + paramap['table_annovar.pl'] + " " + vcf + " " + paramap["anno_dir"] + \
		" -buildver " + paramap['buildver'] + \
		" -protocol " + paramap['protocol'] + \
		" -operation " + paramap['operation'] + \
		" -nastring . " + \
		" --thread " + str(paramap['thread']) + \
		" --maxgenethread " + str(paramap['thread']) + \
		" -polish "

	if vcf.endswith(".vcf") or vcf.endswith(".VCF"):
		#out_name = re.sub(r".vcf$","",out_name)
		#out_name = re.sub(r".VCF$","",out_name)
		out_name = str(sample_name) + "_" + str(v_type)
		vcf_format = True
		cmd = cmd + " -out " + out_name + " -vcfinput "
	else:
		#out_name = re.sub(r".txt$","",out_name)
		out_name = str(sample_name) + "_" + str(v_type)
		vcf_format = False
		cmd = cmd + " -out " + out_name

	print("Run %s\n" % (cmd))
	result = subprocess.run(cmd, stderr=sys.stderr, stdout=sys.stdout, check=True, shell=True)
	#print(result.stdout)
	#print(result.stderr)

	res_file = out_name + "." + paramap['buildver'] + "_multianno.txt"
	return res_file


def main():
	parser = argparse.ArgumentParser(description='Variant annotation')
	parser.add_argument('-i', '--input', default=None, type=str, required=True,
		help="Mapping file")
	parser.add_argument('-d', '--database', default=None, type=str,
		help="Variant annotation database folder")
	parser.add_argument('-a', '--annovar', default=None, type=str,
		help="Annovar table_annovar.pl path")

	parser.add_argument('-o', '--out_dir', default="./", type=str,
		help="Output directory")
	parser.add_argument('-f', '--outfile', default="anno.txt", type=str,
		help="A file contains annotation results")

	parser.add_argument('-b', '--buildver', default="hg38", type=str)
	parser.add_argument('-p', '--protocol', default="refGene", type=str)
	parser.add_argument('-r', '--operation', default="g", type=str)
	parser.add_argument('-c', '--cpu', default=0, type=int)
	
	args = parser.parse_args(sys.argv[1:len(sys.argv)])

	annovar = args.annovar
	map_file = args.input
	anno_res = args.outfile
	anno_dir = args.database
	out_dir = args.out_dir
	buildver = args.buildver
	protocol = args.protocol
	operation = args.operation
	cpu = args.cpu
	if cpu == 0:
		cpu = os.cpu_count()

	paramap = dict()

	paramap['table_annovar.pl'] = annovar
	paramap["anno_dir"] = anno_dir
	paramap['buildver'] = buildver
	paramap['protocol'] = protocol
	paramap['operation'] = operation
	paramap['thread'] = cpu

	if not os.path.exists(out_dir):
		os.makedirs(out_dir)

	dat = pd.read_csv(map_file, sep="\t", header=0, low_memory=False)
	res = dict()

	vtypes = set()
	for i, row in dat.iterrows():
		sample = row["sample"]
		experiment = row['experiment']
		file_row = row["file"]
		file_type_row = row['file_type']

		files = file_row.split(";")
		file_types = file_type_row.split(";")
		for k in range(len(files)):
			vf = files[k]
			vf_type = file_types[k]
			res_file = do_annovar(vf,v_type=vf_type,sample_name=sample,out_dir=out_dir,paramap=paramap)
			if sample in res:
				res[sample][vf_type] = res_file
			else:
				res[sample] = dict()
				res[sample][vf_type] = res_file
			
			vtypes.add(vf_type)

	vtypes = sorted(vtypes)
	of = open(anno_res,"w")
	of.write("sample\t" + "\t".join(vtypes) + "\n")
	for sample in res.keys():
		rf = list()
		for vt in vtypes:
			if vt in res[sample]:
				rf.append(res[sample][vt])
			else:
				rf.append("-")
		of.write(sample + "\t" + "\t".join(rf) + "\n")	
	of.close()


if __name__ == '__main__':
	main()



