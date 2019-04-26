# Variables being inserted

Assc_List=$1
#Output_Dir=$2

#!/usr/bin/env bash

readarray a < $1 # Reading the file into an array

for id in ${a[@]}; do
	/usr/local/sratoolkit/fastq-dump -A $id --gzip --split-files
done
