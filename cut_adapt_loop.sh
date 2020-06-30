### cutadapt loop run with cutadapt 1.18 
### cutadapt options
### -a ADAPTER, --adapter=ADAPTER. Sequence of an adapter that was ligated to the 3' end.
### -g ADAPTER, --front=ADAPTER. Sequence of an adapter that was ligated to the 5' end.
### -A ADAPTER  3' adapter to be removed from second read in a pair.
### -G ADAPTER  5' adapter to be removed from second read in a pair.
### --quality-cutoff 5',3' . Quality threshold used for trimming. It actually applies the algorithm implemented in BWA. It is applied to both reads
### -u number of bases to trimm from 5' in READ1
### -U number of bases to trimm from 5' in READ2
### path to directory containing reads
READS=/home/eugenia/raw_data/orchid_Rawdata
CUT_READS=/home/eugenia/qiime2-ITSxpress-full-dataset/cutreads
SAMPLE_NAMES=/home/eugenia/qiime2-ITSxpress/all_samples_list.txt
CUT_READ1_suf="_1_cut.fastq.gz"
CUT_READ2_suf="_2_cut.fastq.gz"
IN_READ1_suf="_1.fastq.gz"
IN_READ2_suf="_2.fastq.gz"

ls $READS
while read line
do
	echo $line
	sample=$line
	cutadapt --quality-cutoff 20,15 --cores 11 \
	-u 9 \
	-U 9 \
	--adapter GCATATCAATAAGCGGAGGA -A GCTGCGTTCTTCATCGATGC \
	--front GCATCGATGAAGAACGCAGC -G TCCTCCGCTTATTGATATG \
	-o $CUT_READS/$sample$CUT_READ1_suf \
	-p $CUT_READS/$sample$CUT_READ2_suf \
	$READS/$sample$IN_READ1_suf $READS/$sample$IN_READ2_suf > $CUT_READS/$sample"_report.txt"
done < $SAMPLE_NAMES


