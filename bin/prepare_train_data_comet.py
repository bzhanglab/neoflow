
from sklearn.model_selection import train_test_split

import sys

fraction=sys.argv[1]
global_NM_fraction_file=sys.argv[2]
output_path=sys.argv[3]
data_type=sys.argv[4]
print(global_NM_fraction_file)
print(output_path)
with open(global_NM_fraction_file) as data:
    print(global_NM_fraction_file)
    line = data.readline().strip()
    line = data.readline().strip()
    train_file = open(
            output_path + fraction + "_"+ data_type +".txt",
            "w")
    train_file.write("x\ty\tevalue\n")
    while line:
        peptide = line.split("\t")[1]
        modification = line.split("\t")[2]
        rt = line.split("\t")[5]
        evalue = line.split("\t")[3]
        
        if "V" in modification:
            mod_list = modification.split(";")
            for each_mod in mod_list:
                if "V" in each_mod:
                    print("M oxidation")
                    index = int(each_mod.split("_V")[0])
                    peptide = peptide[:index - 1] + "1" + peptide[index:]

        train_file.write(peptide + "\t" + rt + "\t" + evalue + "\n")
        line = data.readline().strip()
