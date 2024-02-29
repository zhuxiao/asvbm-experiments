#!/bin/bash
if [ $# -ne 2 ]
then
    echo "Usage : $0 <input_file> <output_file>"
    echo "    input_file   TXT"
    echo "    output_file  FASTQ"
    #echo "      example A: ./$0 SRR_Acc_List.txt SRR885_whole.fastq"
    #echo "      example B: bash $0 SRR_Acc_List.txt SRR885_whole.fastq"
    exit 1
fi

input_file=$1
output_file=$2
> "$output_file"

if [ ! -f "$input_file" ]; then
    echo "Input file not found."
    exit 1
fi
mkdir SRX5327410
while IFS= read -r line;do
    cd SRX5327410
    echo "$line downloading..." 
    prefetch $line
    echo "$line download successfully"
    cd $line
    echo "fasterq-dump $line"
    fasterq-dump $line.sra
    echo "done."
    echo "cat starting"
    cat $line.fastq >> ../../$output_file
    echo "cat done"
    cd ../../
done < "$input_file"
echo "All file are download and dump successfully"


