
#!/bin/bash

# usage: map_small_RNAs_to_ref_and_denovo_assemble_unmppd.sh <threads(int)> <memory(int)> <reference.fa> <reads_suffix>

# For a set of fastq files, the script performs adapter trimming and QC, maps to a prokaryotic reference, collects unmapped reads, and denovo assembles these.

if [[ -n $1 ]] && [[ -n $2 ]] && [[ -n $3 ]] && [[ -n $4 ]]
then
  echo "number of processing cores requested: "$1
  echo "amount of memory requested: "$2
  echo "reference genome file: "$3
  echo "using input read files with "$4" suffix"
  echo ""; echo ""

  #load env variables, etc
  module load Java
  module load trimmomatic
  module load fastqc
  module load bowtie2
  module load samtools
  module load bedtools
  module load trinity

  #make a bowtie2 prefix
  idxprefix=$(echo "$3" | cut -f 1 -d '.')

  #create an array of fastq filenames
  shopt -s nullglob
  inreads=(*.$4)

  #create memory allocation argument for trinity
  mem=$2"G"

  #print number of read files detected, and list them.
  echo "found ${#inreads[@]} read files in working directory"
  for i in "${inreads[@]}"; do echo "$i"; done

  #iterate over inreads array
  for i in "${inreads[@]}"
  do
	  #make prefix for readfile
	  readspref=$(echo $i | cut -f 1 -d ".")

	  #make filenames for trimmed fastqs
	  trimmed=$readspref".trimmed.fq"

	  #make prefix for sam, bam and txt files
	  aligned=$readspref"_to_miRs"
	  #run the pipeline

	  #Trim reads for adapter sequences
	  java -jar /opt/biotools/trimmomatic/trimmomatic-0.35.jar SE -threads $1 -quiet\
	  $i \
	  $trimmed\
	  ILLUMINACLIP:/opt/biotools/trimmomatic/adapters/TruSeq3-SE.fa:2:30:9
	  trimx=$?
	  echo ""; echo ""

	  #generate QC report
	  if [ $trimx -eq 0 ]
	  then
	    echo "running fastqc..."
	    fastqc $trimmed
	    fqcx=$?
	    echo ""; echo ""
	  else
	    echo "trimmomatic failed"
	  fi

	  #map reads to reference, convert to bam, keep only unmapped reads, and sort by query name
	  if [[ $trimx -eq 0 ]] && [[ $fqc -eq 0 ]] #&& [[ $filx -eq 0 ]]
	  then
	    echo "generating sorted bam file..."
	    btcores=$(expr $1 - 1)
	    bowtie2 -p $btcores --fast-local -x $idxprefix -U $trimmed | samtools view -h -@ 1 -uT $3 > $aligned.bam
		samtools view -h -@ $1 -f 4 $aligned.bam > $aligned"_unm.bam"
		samtools sort -n -@ $1 $aligned"_unm.bam" > $aligned"_unm_sorted.bam"; rm $aligned"_unm.bam" $aligned.bam
	    btstx=$?
	    echo ""; echo ""
	  else
	    echo "fastqc failed"
	  fi

	  #convert bam to fastq
	  if [[ $trimx -eq 0 ]] && [[ $fqc -eq 0 ]] && [[ $btstx -eq 0 ]] && [[ $trinx -eq 0 ]]
	  then
	    echo "converting bam to fastq"
	  	bedtools bamtofastq -i $aligned"_unm_sorted.bam" -fq $aligned"_unm.fq"
	  	bedx=$?
	  	echo ""; echo ""
	  else
	  	echo "bowtie2 or samtools failed"
	  fi

	  #make denovo assemblies
	  if [[ $trimx -eq 0 ]] && [[ $fqc -eq 0 ]] && [[ $btstx -eq 0 ]] && [[ $bedx -eq 0 ]]
	  then
	    echo "creating denovo assemblies"
	  	Trinity --seqType fq --CPU $1 --max_memory $mem --single $aligned"_unm.fq" --output $readspref"_trinity_out"
		trinx=$?
		echo ""; echo ""
	  else
		echo "bedtools failed"
	  fi

	  if [[ $trimx -eq 0 ]] && [[ $fqc -eq 0 ]] && [[ $btstx -eq 0 ]] && [[ $bedx -eq 0 ]] && [[ $trinx -eq 0 ]]
	  then
		  echo "denovo assembly completed"
	  else
		  echo "trinity failed"
  done
	#move temp files
	mkdir temp_files
	mv $trimmed temp_files
	mv $aligned"_unm.fq" temp_files
	mv $aligned"_unm_sorted.bam" temp_files
	echo "job completed"
fi

else
  echo "An argument is missing"
fi
