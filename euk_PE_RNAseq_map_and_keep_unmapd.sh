#!/bin/bash

#usage: euk_PE_RNAseq_map_and_keep_unmapd.sh <threads> <genome.fa> <readfile_list.txt>

# the script aligns paired RNA-seq fastq files to reference and keeps only the unmapped reads, saved in bam format

# readfile_list.txt example:
# sample_1 s1_forward.fq s1_reverse.fq
# sample_2 s2_forward.fq s2_reverse.fq

if [[ -n $1 ]] && [[ -n $2 ]] && [[ -n $3 ]]
then
    echo ""; echo ""
    echo "number of processing cores requested: "$1
    echo "reference fasta filename: "$2
    echo "list of input fastq filenames: "$3
    echo ""; echo ""


    #load env variables, etc
    module load trimmomatic
    module load fastqc
    module load Java
    module load bowtie2
    module load tophat

    #take prefixes from reference filnames
    refprefix=$(echo "$2" | cut -f 1 -d '.')
    counter=0
    mkdir fastQC_files

    #assign call and options to variable for use in trimmomatic
    call_trim="java -jar /opt/biotools/trimmomatic/trimmomatic-0.35.jar PE -threads $1 -phred33"
    trimoptions="ILLUMINACLIP:/opt/biotools/trimmomatic/adapters/TruSeq3-PE.fa:2:30:10 SLIDINGWINDOW:4:20 MINLEN:20" \

    #for each line of text file, split into 3 strings with whitespace delimiter, and assign each to sample, read1, and read2
    input=$3
    while IFS= read -r line; do
        sample=$(echo $line | cut -f 1 -d ' ')
        read1=$(echo $line | cut -f 2 -d ' ')
        read2=$(echo $line | cut -f 3 -d ' ')
        counter=$(expr $counter + 1)
        echo "Processing sample No."$counter
        echo $sample

        #Trim reads for adapter sequences
        $call_trim \
        $read1 $read2 \
        $read1".trimmed.paired" $read1".trimmed.unpaired" \
        $read2".trimmed.paired" $read2".trimmed.unpaired" \
        $trimoptions; trimx=$?
        echo ""; echo ""

        #generate QC report for paired reads
        if [ $trimx -eq 0 ]
        then
            echo "running fastqc on "$read1".trimmed.paired"; fastqc -t $1 -q -o ./fastQC_files -q $read1".trimmed.paired"
            echo "running fastqc on "$read2".trimmed.paired"; fastqc -t $1 -q -o ./fastQC_files -q $read2".trimmed.paired"
            echo "running fastqc on "$read1".trimmed.unpaired"; fastqc -t $1 -q -o ./fastQC_files -q $read1".trimmed.unpaired"
            echo "running fastqc on "$read2".trimmed.unpaired"; fastqc -t $1 -q -o ./fastQC_files -q $read2".trimmed.unpaired"
            fqcx=$?; echo ""; echo ""
        else
            echo "trimmomatic failed"; echo ""; echo ""
            break 2
        fi

        #map reads to host, collect unmapped reads and lof files
        if [[ $fqcx -eq 0 ]]
        then
            tophat2 -o tophat_temp -p $1 --b2-fast --max-intron-length 5000 --max-segment-intron 5000 $refprefix \
            $read1".trimmed.paired",$read1".trimmed.unpaired" $read2".trimmed.paired",$read2".trimmed.unpaired"; thx=$?
            echo ""; echo ""

        else
            echo "fastqc failed"; echo ""; echo ""
            break 2
        fi


        if [[ $thx -eq 0 ]]
        then
            mv ./tophat_temp/unmapped.bam ./tophat_temp/$sample"_unm.bam"; mv ./tophat_temp/$sample"_unm.bam" ./
            mv ./tophat_temp/logs ./tophat_temp/$sample"_tophat_logs"; mv ./tophat_temp/$sample"_tophat_logs" ./
            echo ""; echo ""

        else
            echo "tophat2 failed"; echo ""; echo ""
            break 2
        fi

        echo "Finished processing sample No."$counter

    done < "$input"

else
  echo "An argument is missing"
fi
