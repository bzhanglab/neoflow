import os
import sys
import glob
import re
import subprocess

input_data_dir = sys.argv[1]
out_dir = sys.argv[2]

if not os.path.exists(out_dir):
    os.mkdir(out_dir)


training_data_files = glob.glob(str(input_data_dir)+"/train_data/*_train.tsv")

autort_path = "/opt/AutoRT/autort.py"
base_model_dir = "/opt/AutoRT/models/base_models_PXD006109/"

for train_file in training_data_files:
    file_name = os.path.basename(train_file)
    data_dir  = os.path.dirname(train_file)

    sample_name = re.sub("_train.tsv$","",file_name)
    test_file = str(data_dir) + "/" + str(sample_name) + "_test.tsv"
    o_dir = str(out_dir) + "/" + sample_name
    if not os.path.exists(o_dir):
        os.mkdir(o_dir)

    new_model_dir = str(o_dir) + "/models"
    if not os.path.exists(new_model_dir):
        os.mkdir(new_model_dir)

    pred_dir = str(o_dir) + "/prediction"
    if not os.path.exists(pred_dir):
        os.mkdir(pred_dir)

    train_cmd = "python "+ autort_path +" train -i "+ str(train_file) + " -o " + str(new_model_dir) + " -e 40 -b 64 -u m -m "+ str(base_model_dir)+"/model.json -rlr -n 10"
    print(train_cmd)
    try:
        rv = subprocess.run([train_cmd],shell=True,stdout=subprocess.PIPE, stderr=subprocess.STDOUT, check=True)
    except subprocess.CalledProcessError as sub_error:
        print("Training error: %s => %s" % (file_name, sub_error.stdout.decode('utf-8')),file=sys.stderr)
        sys.exit(1)

    ## testing
    pred_cmd = "python " + autort_path + " predict -t "+ str(test_file) + " -s "+ str(new_model_dir) + "/model.json -o "+ str(pred_dir) + " -p test"
    print(pred_cmd)
    try:
        rv = subprocess.run([pred_cmd],shell=True,stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError as sub_error:
        print("Prediction error: %s => %s" % (file_name, sub_error.stdout.decode('utf-8')),file=sys.stderr)
        sys.exit(1)
