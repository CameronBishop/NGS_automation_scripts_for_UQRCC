#!/bin/bash/

# USAGE: run the script with no arguments. Script will use as input any file in the working dir that has a ".fa" extension.
# It will create a two-column file with each row representing a read and a NCBI taxonomic identifier for use in Kronatools.

#create an array of fasta filenames
shopt -s nullglob
inreads=(*.fa)
echo "FILES TO BE USED AS INPUT:"
for i in "${inreads[@]}"; do echo "$i"; done


for i in "${inreads[@]}"
do
######################### RUN KRAKEN, CREATE KRONA INPUT, AND RUN KRONA ##################################
readspref=$(echo $i | cut -f 1 -d ".")
report=$readspref".report.txt"
krakenout=$readspref".krakenout.txt"
kronain=$readspref"_krona_input.txt"
kronaout=$readspref"_krona.html"
kraken2 --threads 6 --db /home/cam/kraken/minikraken2_v2_8GB_201904_UPDATE --report $report $i > $krakenout
grep "^C" $krakenout | cut -f 2,3 > $kronain
ktImportTaxonomy $kronain -o $kronaout
echo $i" completed"
done
echo "all jobs completed"


######################### REMOVE WOLBACHIA READS, AND RUN KRONA AGAIN ####################################
shopt -s nullglob
kronin=(*krona_input.txt)
echo "FILES TO BE USED AS INPUT:"
for i in "${kronin[@]}"; do echo "$i"; done

#iterate over inreads array
for i in "${kronin[@]}"
do
#
noWolbprefix=$(echo $i | cut -f 1 -d "_")
noWolbkronain=$noWolbprefix"no_Wolb_krona_input.txt"
noWolbkronaout=$noWolbprefix"no_Wolb_krona.html"
awk -F "\t" '{ if ( $2 != 953 && $2 != 570417 && $2 != 1236908 && $2 != 1236909 && $2 != 225364 && $2 != 77038 && $2 != 125593 && $2 != 169402 && $2 != 66084 && $2 != 292805 && $2 != 246273 && $2 != 1633785 ) {print} }' $i > $noWolbkronain
ktImportTaxonomy $noWolbkronain -o $noWolbkronaout
echo $i" Wolbachia reads removed, and krona plot created"
done
echo "all jobs completed"
