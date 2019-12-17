#!/bin/sh
####For tmt

set -e -x

mgf_file=${1}
database=${2}
output=${3}
mod_file=${4}

java -Xmx48g -jar /opt/MSGFPlus.jar \
	-s $mgf_file \
        -d $database \
        -mod $mod_file \
        -o $output \
	-m 3 \
        -maxLength 45 \
	-n 1 \
        -minlength 7 \
	-addFeatures 1 \
	-t 20.0ppm \
	-minCharge 2 \
        -maxCharge 3 \
        -tda 0 \
	-protocol 4 \
	-ntt 2 \
        -inst 3 \
        -ti 0,0 \
        -e 1 \
	-maxMissedCleavages 2
