Raw_Data=$1
Processed_Data=$2
nthreads=$3

#####################
###Quality Control###
#####################

echo "FastQC"

mkdir -p Quality_Check

parallel -j $nthreads "fastqc --noextract -o Quality_Check {}" ::: $Raw_Data/*.fastq.gz

########################
###Trimming the Files###
########################

echo "Trimming"

mkdir -p Trimmed_out

#Performing the loop
for f in $Raw_Data/*R1*; do
	Raw_Out=$(echo Trimmed_out/$(basename $f))
	Log_Out=$(echo Trimmed_out/$(basename $f ".fastq.gz").trimlog)
	java -jar /usr/local/trimmomatic/trimmomatic-0.36.jar PE -threads $nthreads -phred33 -trimlog $Log_Out -basein $f -baseout $Raw_Out ILLUMINACLIP:/usr/local/trimmomatic/adapters/TruSeq3-PE-2.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
done

########################
###Paired End Merging###
########################
echo "Now merging our Paired End Reads"

mkdir -p merge_int
mkdir -p merge_out
parallel -j $nthreads "/usr/local/seqprep/SeqPrep -f {.}_1P.fastq.gz -r {.}_2P.fastq.gz -1 /dev/null -2 /dev/null -s merge_int/{/.}.fastq.gz;
cat merge_int/{/.}.fastq.gz {.}_1U.fastq.gz > merge_out/{/.}.fastq.gz" ::: Trimmed_out/*.trimlog
rm -rf merge_int

#########################
###Removing Duplicates###
#########################
echo "Deduplicating"

mkdir -p Dedup_out
parallel -j $nthreads "zcat {} | perl ~/Applications/prinseq-lite-0.20.4/prinseq-lite.pl -fastq stdin -derep 14 -out_good stdout -out_bad null 2>/dev/null | gzip > Dedup_out/{/}" ::: merge_out/*.fastq.gz

############################
###Metagenomics or 16sRNA###
############################
mkdir -p kraken_uniq_out

for f in Dedup_out/*fastq.gz; do
	/home/keaton/myapps/krakenuniq-0.5.7/krakenuniq --db /2/scratch/keaton/Kraken --threads $nthreads --fastq-input --gzip-compressed --output kraken_uniq_out/$(basename $f ".fastq.gz").out --report-file kraken_uniq_out/$(basename $f ".fastq.gz").tab $f #Running Kraken
done

mkdir -p kraken_MPA
for f in kraken_uniq_out/*.out; do
	/home/keaton/myapps/krakenuniq-0.5.7/krakenuniq-mpa-report --db /2/scratch/keaton/Kraken $f > kraken_MPA/$(basename $f .out).mpa
done
