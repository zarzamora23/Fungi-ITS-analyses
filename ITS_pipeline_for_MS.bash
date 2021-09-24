### Created by E.Zarza in 2020 to analyse ITS data for manuscript: Fungal diversity in shade-coffee plantations in Soconusco, Mexico
### Eugenia Zarza, Alejandra López-Pastrana, Anne Damon, Griselda Karina Guillén-Navarro, Luz Verónica García-Fajardo. 
### El Colegio de la Frontera Sur - Tapachula - Mexico
### Raw sequences are available form Genbank bioproject: PRJNA610266 (At publication time)
### Analyses were executed on the Linux Subsytem for Windows. It requires QIIME2, R, vegan and other libraries indicated bellow.
### This is not intended to run as script. Commands should be executed one by one. Main reason for this is that it is necessary to check outcome of applied filters.
### These notes start indicating the database for taxonomy assignment, followed by importing data to QIIME2, taxonomic assignment, diversity calculations considering all taxa included in UNITE, heatmap creation, GUILDFun analysis and CLAM test 
### Define some variables if needed
my_online_folder=path_to_folder_available_online
### Define path to directory where you may store raw data or results under Windows OS, usually has a path starting /mnt/c/Users/ ...etc. 
my_windows_folder=path_to_directory_in_windows_OS 

### Follow ITSxpress pipeline: https://forum.qiime2.org/t/q2-itsxpress-a-tutorial-on-a-qiime-2-plugin-to-trim-its-sequences/5780
### Rivers AR, Weber KC, Gardner TG, Liu S, Armstrong SD. ITSxpress: Software to rapidly trim internally transcribed spacer sequences with quality scores for marker gene analysis. 
### F1000Res. 2018;7:1418. Published 2018 Sep 6. doi:10.12688/f1000research.15704.1
source activate qiime2-2019

### Install ITSxpress

pip install q2-itsxpress

### copy data to linux
cp $my_windows_folder/Raw_Data/ $HOME/raw_data/orchid_Rawdata

### do fastqc analysis
#http://lab.loman.net/high-throughput%20sequencing/genomics/bioinformatics/2013/04/17/adaptor-trim-or-die-experiences-with-nextera-libraries/

fastqc $HOME/raw_data/orchid_Rawdata/*.gz 

### run cutadapt loop to remove adapters

bash cut_adapt_loop.sh

### do fastqc analysis again to verify adapter removal. 
cd $HOME/qiime2-ITSxpress-full-dataset/cutreads

fastqc *.gz 

### summarize with muliqc

multiqc .

### change to working directory
cd $HOME/qiime2-ITSxpress-full-dataset
### create manifest

### extract sample name

ls $HOME/raw_data/orchid_Rawdata | grep '.gz' | awk -F $'_' '{print $1",$HOME/raw_data/orchid_Rawdata/" $0 ","}' | sort -d |  sed -e 's/_1.fastq.gz,/_1.fastq.gz,forward/g' -e 's/_2.fastq.gz,/_2.fastq.gz,reverse/g' -e '1 i\sample-id,absolute-filepath,direction' > orchid_manifest.txt

### import samples produced with cutadapt, path indicated in orchid_manifest.txt 

qiime tools import \
  --type SampleData[PairedEndSequencesWithQuality] \
  --input-format PairedEndFastqManifestPhred33 \
  --input-path orchid_manifest.txt \
  --output-path sequences.qza
  
  qiime demux summarize \
  --i-data sequences.qza \
  --o-visualization sequences.qzv
  
cp sequences.qzv $my_windows_folder/results_qiime2

### run itsxpress, specify sequenced region. According to macrogen provided document OTUanalysis.html, a fragment of 300bp was obtained, ITS3-4, 
### ITSxpress article (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6206612/) mentions  " ITS2 reads were amplified using the ITS3/ITS4 primer set ( White et al., 1990)".
### Thus I assume ITS2 was amplified. To confirm, I performed a BLAST search of sequence @M01406:32:000000000-BT455:1:1102:11492:1276 1:N:0:GGACTCCT+ACTGCATA, from file BOP_1_1
### result confirms it belongs to ITS2: Uncultured fungus clone OTU_F131_R184 small subunit ribosomal RNA gene, partial sequence; internal transcribed spacer 1, 5.8S ribosomal RNA gene, and internal transcribed spacer 2, complete sequence; and large subunit ribosomal RNA gene, partial sequence 

### As previous macrogen analyses suggested that the sequences contained plants and animals, I selected ALL taxa

qiime itsxpress trim-pair-output-unmerged \
  --i-per-sample-sequences sequences.qza \
  --p-region ITS2 \
  --p-taxa ALL \
  --p-threads 11 \
  --verbose \
  --o-trimmed trimmed.qza
  
### According to https://bmcmicrobiol.biomedcentral.com/articles/10.1186/1471-2180-10-189 the position of the primers and thus the expected amplicon length are:
### ITS3 2024-2045; ITS4 2390-2409; Expected length 385

####TRY ITSEXPRESS, with 6 ERRORS, 
 ### In this thread Benjamin, a dada2 developer discusses 
 ### "The biggest benefit of a maxEE filter is it reduces computation time, by removing a lot of reads with unique sequences because of their many errors. 
 ### This is recommended, because reads with errors don't improve the inference of what is really there. 
 ### Additionally, in some unusual cases there are pathological artefacts that show up in very low quality reads which are best to get rid of this way (eg. non-target sequences).
 ### However, because DADA2 incorporates the quality scores in its error model, it is quite robust to lower quality sequences (for example). 
 ### I've never seen it "break down", so if you need to raise the maxEE cutoff to get more reads through, you can do so. https://github.com/benjjneb/dada2/issues/248 """
 ### check this other thread too https://github.com/benjjneb/dada2/issues/232  
### analyze ITSxpress trimmed sequences, try --p-max-ee 6 [default 2.0] . For this samples, the prefic full_dataset prefix was used because I had done analyses on a subsampled dataset

qiime dada2 denoise-paired \
  --i-demultiplexed-seqs trimmed.qza \
  --p-trunc-len-r 0 \
  --p-trunc-len-f 0 \
  --p-max-ee 6 \
  --p-n-threads 11 \
  --verbose \
  --output-dir full_dataset_ITSx_dada2out
  
  qiime feature-table summarize \
  --i-table full_dataset_ITSx_dada2out/table.qza \
  --o-visualization full_dataset_ITSx_tableviz.qzv
  
  ### create summaries for the output of previous analyses
 
 ### to view stats 
  qiime metadata tabulate \
  --m-input-file full_dataset_ITSx_dada2out/denoising_stats.qza \
  --o-visualization full_dataset_ITSx_stats-dada2.qzv

### to view representative sequences
  qiime feature-table tabulate-seqs \
  --i-data full_dataset_ITSx_dada2out/representative_sequences.qza \
  --o-visualization full_dataset_ITSx_rep-seqs.qzv
 
 cp full_dataset_ITSx*.qzv $my_windows_folder/results_qiime2
 

################################################## TAXONOMIC ASSIGNMENT ################################################################ 
### Get database for TAXONOMIC ASSIGNMENT with UNITE DATABASE
### mkdir for database
mkdir ITS_database
cd ITS_database
### Download ITS UNITE database
### Following Kõljalg et al. (2013), each terminal fungal taxon for which two or more ITS sequences are available is referred to as a species hypothesis (SH).
### One sequence is chosen to represent each SH;
### these sequences are called representative sequences (RepS) when chosen automatically by the computer and reference sequences (RefS) when those choices are overridden (or confirmed) by users with expert knowledge of the taxon at hand.
### Three sets of QIIME files are released, corresponding to the SHs resulting from clustering at the 97% and 99% threshold levels.
### The third set of files is the result of a dynamic use of clustering thresholds, such that some SHs are delimited at the 97% level, some at the 97.5% level, some at the 98% level, and so on; these choices were made manually by experts of those particular lineages of fungi. The syntax is the same throughout the three sets of files. 
### Includes singletons set as RefS (in dynamic files).

### I downloaded several databases, but used only one 'All eukaryotes	9 409 (RefS)	54 013 (RepS)'

### 8.0	2018-11-18	All eukaryotes	9 409 (RefS)	54 013 (RepS)	Current	https://doi.org/10.15156/BIO/786335	When using this resource, please cite it as follows:
### UNITE Community (2019): UNITE QIIME release for eukaryotes. Version 18.11.2018. UNITE Community. https://doi.org/10.15156/BIO/786335 
### Includes singletons set as RefS (in dynamic files).
wget https://files.plutof.ut.ee/public/orig/6D/A3/6DA3F3DDFBACB7D3031FE12EF14C4F5464235C73230C217D14E8C0EAB1042E86.zip

### Unzip database for taxonomy analyses in database directory

unzip 6DA3F3DDFBACB7D3031FE12EF14C4F5464235C73230C217D14E8C0EAB1042E86.zip
  
### import into qiime2. The trimmed database is recommended in the tutorial (untrimmed database is included in the develope folder)
qiime tools import \
  --type 'FeatureData[Sequence]' \
  --input-path $HOME/ITS_database/sh_refs_qiime_ver8_dynamic_all_02.02.2019.fasta \
  --output-path unite.qza

### import taxonomyImport the associated UNITE taxonomy file.
qiime tools import \
  --type 'FeatureData[Taxonomy]' \
  --input-format HeaderlessTSVTaxonomyFormat \
  --input-path $HOME/ITS_database/sh_taxonomy_qiime_ver8_dynamic_all_02.02.2019.txt \
  --output-path unite-taxonomy.qza

### Train the classifier as: classifier.qza
qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads unite.qza \
  --i-reference-taxonomy unite-taxonomy.qza \
  --o-classifier classifier.qza

### Once the classifier is trained sequences can be classified.
qiime feature-classifier classify-sklearn \
  --i-classifier classifier.qza \
  --i-reads full_dataset_ITSx_dada2out/representative_sequences.qza \
  --o-classification full_dataset_ITSx_taxonomy.qza

### Summarize the results
### Summarize the results for visualization in the QIIME 2 viewer.

qiime metadata tabulate \
  --m-input-file full_dataset_ITSx_taxonomy.qza \
  --o-visualization full_dataset_ITSx_taxonomy.qzv

### Create an interactive bar plot figure
qiime taxa barplot \
  --i-table full_dataset_ITSx_dada2out/table.qza  \
  --i-taxonomy full_dataset_ITSx_taxonomy.qza \
  --m-metadata-file orchid_metadata.tsv \
  --o-visualization full_dataset_ITSx_taxa-bar-plots.qzv
cp full_dataset_ITSx_taxa-bar-plots.qzv $my_windows_folder/results_qiime2
cp full_dataset_ITSx_taxonomy.qzv $my_windows_folder/results_qiime2
 
 ### generate a tree for phylogenetic analyses, but as it contains not only fungi, would not be very informative
 
 #qiime phylogeny align-to-tree-mafft-fasttree \
  #--i-sequences full_dataset_ITSx_dada2out/representative_sequences.qza \
  #--o-alignment full_dataset_ITSx_aligned-rep-seqs.qza \
  #--o-masked-alignment full_dataset_ITSx_masked-aligned-rep-seqs.qza \
  #--o-tree full_dataset_ITSx_unrooted-tree.qza \
  #--o-rooted-tree full_dataset_ITSx_rooted-tree.qza
 
### this threshold was set because it represents a drop to 0.85 from previous value 
  qiime diversity core-metrics-phylogenetic \
  --i-phylogeny full_dataset_ITSx_rooted-tree.qza \
  --i-table full_dataset_ITSx_dada2out/table.qza \
  --p-sampling-depth 77000 \
  --m-metadata-file orchid_metadata.tsv \
  --output-dir full_dataset_ITSx_core-metrics-results
  
  cp full_dataset_ITSx_core-metrics-results/*.qzv $my_windows_folder/results_qiime2
  
  cp full_dataset_ITSx_core-metrics-results $my_windows_folder/results_qiime2


### from moving-pictures tutorial: --p-max-depth should be determined by reviewing the “Frequency per sample” information presented in the table.qzv file that was created above. 
  ### I choosed a value around the median 122616
 
  qiime diversity alpha-rarefaction \
  --i-table full_dataset_ITSx_dada2out/table.qza \
  --i-phylogeny full_dataset_ITSx_rooted-tree.qza \
  --p-max-depth 123000 \
  --m-metadata-file orchid_metadata.tsv \
  --o-visualization full_dataset_alpha-rarefaction.qzv

### alpha group significance
  qiime diversity alpha-group-significance \
  --i-alpha-diversity full_dataset_ITSx_core-metrics-results/shannon_vector.qza \
  --m-metadata-file orchid_metadata.tsv \
  --o-visualization full_dataset_ITSx_core-metrics-results/shannon-group-significance.qzv
 
  cp full_dataset_ITSx_core-metrics-results/*.qzv $my_windows_folder/results_qiime2
  
### create emperor plot 
qiime emperor plot \
  --i-pcoa full_dataset_ITSx_core-metrics-results/bray_curtis_pcoa_results.qza \
  --m-metadata-file orchid_metadata.tsv \
  --o-visualization full_dataset_ITSx_core-metrics-results/bray-curtis-emperor_pcoa.qzv
  
  cp full_dataset_ITSx_core-metrics-results/bray-curtis-emperor_pcoa.qzv $my_windows_folder/results_qiime2
  
###perform bray curtis group association because it is more appropiate for this kind of data
### beta group significance bray_curtis
qiime diversity beta-group-significance \
  --i-distance-matrix full_dataset_ITSx_core-metrics-results/bray_curtis_distance_matrix.qza \
  --m-metadata-file orchid_metadata.tsv \
  --m-metadata-column Type \
  --o-visualization full_dataset_ITSx_core-metrics-results/bray_curtis-Type-group-significance.qzv \
  --p-pairwise

qiime diversity beta-group-significance \
  --i-distance-matrix full_dataset_ITSx_core-metrics-results/bray_curtis_distance_matrix.qza \
  --m-metadata-file orchid_metadata.tsv \
  --m-metadata-column Site \
  --o-visualization full_dataset_ITSx_core-metrics-results/bray_curtis-site-group-significance.qzv \
  --p-pairwise
  
 qiime diversity beta-group-significance \
  --i-distance-matrix full_dataset_ITSx_core-metrics-results/bray_curtis_distance_matrix.qza \
  --m-metadata-file orchid_metadata.tsv \
  --m-metadata-column Group1 \
  --o-visualization full_dataset_ITSx_core-metrics-results/bray_curtis-Group1-group-significance.qzv \
  --p-pairwise
  
  qiime diversity beta-group-significance \
  --i-distance-matrix full_dataset_ITSx_core-metrics-results/bray_curtis_distance_matrix.qza \
  --m-metadata-file orchid_metadata.tsv \
  --m-metadata-column Group2 \
  --o-visualization full_dataset_ITSx_core-metrics-results/bray_curtis-Group2-group-significance.qzv \
  --p-pairwise
  
######################################################### ANALYZE ONLY FUNGI ###################################################################

qiime itsxpress trim-pair-output-unmerged \
  --i-per-sample-sequences sequences.qza \
  --p-region ITS2 \
  --p-taxa F \
  --p-threads 11 \
  --verbose \
  --o-trimmed full_Fungi_sequences_trimmed.qza
  
### analyze ITSxpress trimmed sequences, try --p-max-ee 6 [default 2.0]

qiime dada2 denoise-paired \
  --i-demultiplexed-seqs full_Fungi_sequences_trimmed.qza \
  --p-trunc-len-r 0 \
  --p-trunc-len-f 0 \
  --p-max-ee 6 \
  --p-n-threads 11 \
  --verbose \
  --output-dir full_Fungi_ITSx_dada2out
  
  qiime feature-table summarize \
  --i-table full_Fungi_ITSx_dada2out/table.qza \
  --o-visualization full_Fungi_ITSx_tableviz.qzv
  
### create summaries for the output of previous analyses
 
 ### view stats 
  qiime metadata tabulate \
  --m-input-file full_Fungi_ITSx_dada2out/denoising_stats.qza \
  --o-visualization full_Fungi_ITSx_stats-dada2.qzv

### view representative sequences
  qiime feature-table tabulate-seqs \
  --i-data full_Fungi_ITSx_dada2out/representative_sequences.qza \
  --o-visualization full_Fungi_ITSx_rep-seqs.qzv
 
 cp full_Fungi_ITSx*.qzv $my_windows_folder/results_qiime2
 
 
### Database was unzipped and imported above
### Classified was already trained above
### Once the classifier is trained sequences can be classified.

qiime feature-classifier classify-sklearn \
  --i-classifier classifier.qza \
  --i-reads full_Fungi_ITSx_dada2out/representative_sequences.qza \
  --o-classification full_Fungi_ITSx_taxonomy.qza

### Summarize the results
### Summarize the results for visualization in the QIIME 2 viewer.

qiime metadata tabulate \
  --m-input-file full_Fungi_ITSx_taxonomy.qza \
  --o-visualization full_Fungi_ITSx_taxonomy.qzv

### Create an interactive bar plot figure
qiime taxa barplot \
  --i-table full_Fungi_ITSx_dada2out/table.qza  \
  --i-taxonomy full_Fungi_ITSx_taxonomy.qza \
  --m-metadata-file orchid_metadata.tsv \
  --o-visualization full_Fungi_ITSx_taxa-bar-plots.qzv
cp full_Fungi_ITSx_taxa-bar-plots.qzv $my_windows_folder/results_qiime2
cp full_Fungi_ITSx_taxonomy.qzv $my_windows_folder/results_qiime2

----
qiime taxa filter-table \
  --i-table full_Fungi_ITSx_dada2out/table.qza \
  --i-taxonomy full_Fungi_ITSx_taxonomy.qza \
  --p-include k__Fungi \
  --o-filtered-table filtered-Fungi-table.qza

qiime metadata tabulate \
  --m-input-file filtered-Fungi-table.qza \
  --o-visualization filtered-Fungi-table.qzv


### barplot visualization of taxonomy
qiime taxa barplot \
  --i-table filtered-Fungi-table.qza \
  --i-taxonomy full_Fungi_ITSx_taxonomy.qza \
  --m-metadata-file orchid_metadata.tsv \
  --o-visualization filtered_Fungi_ITSx_taxa-bar-plots.qzv
cp filtered_Fungi_ITSx_taxa-bar-plots.qzv $my_windows_folder/results_qiime2

###################################TO RUN AFTER TAXONOMIC ANALYSIS####################################################
### some hints to apply the filter and execute the analyses came from https://forum.qiime2.org/t/compute-alpha-beta-diversity-at-different-taxonomic-levels/2615

### this threshold was set because it represents a drop to 0.84 from previous value (under 0.85), but this value was estimated prior to removing unassigned taxa.
  qiime diversity core-metrics \
  --i-table filtered-Fungi-table.qza \
  --p-sampling-depth 20000 \
  --m-metadata-file orchid_metadata.tsv \
  --output-dir filtered_Fungi_ITSx_core-metrics-results
cp -r filtered_Fungi_ITSx_core-metrics-results $my_windows_folder/results_qiime2 

### calculated and copied on sept 04-09-2019 
qiime diversity alpha-rarefaction \
  --i-table filtered-Fungi-table.qza \
  --p-max-depth 20000 \
  --m-metadata-file orchid_metadata.tsv \
  --o-visualization filtered_Fungi_ITSx_core-metrics-results/alpha-rarefaction.qzv
cp filtered_Fungi_ITSx_core-metrics-results/alpha-rarefaction.qzv $my_windows_folder/results_qiime2/filtered_Fungi_ITSx_core-metrics-results

### alpha group significance
 qiime diversity alpha-group-significance \
  --i-alpha-diversity filtered_Fungi_ITSx_core-metrics-results/shannon_vector.qza \
  --m-metadata-file orchid_metadata.tsv \
  --o-visualization filtered_Fungi_ITSx_core-metrics-results/shannon-group-significance.qzv
   
cp filtered_Fungi_ITSx_core-metrics-results/*.qzv $my_windows_folder/results_qiime2
   
###perform bray curtis group association because it is more appropiate for this kind of data 
### beta group significance Bray Curtis
qiime diversity beta-group-significance \
  --i-distance-matrix filtered_Fungi_ITSx_core-metrics-results/bray_curtis_distance_matrix.qza \
  --m-metadata-file orchid_metadata.tsv \
  --m-metadata-column Type \
  --o-visualization filtered_Fungi_ITSx_core-metrics-results/bray_curtis-Type-group-significance.qzv \
  --p-pairwise

qiime diversity beta-group-significance \
  --i-distance-matrix filtered_Fungi_ITSx_core-metrics-results/bray_curtis_distance_matrix.qza \
  --m-metadata-file orchid_metadata.tsv \
  --m-metadata-column Site \
  --o-visualization filtered_Fungi_ITSx_core-metrics-results/bray_curtis-site-group-significance.qzv \
  --p-pairwise
  
 qiime diversity beta-group-significance \
  --i-distance-matrix filtered_Fungi_ITSx_core-metrics-results/bray_curtis_distance_matrix.qza \
  --m-metadata-file orchid_metadata.tsv \
  --m-metadata-column Group1 \
  --o-visualization filtered_Fungi_ITSx_core-metrics-results/bray_curtis-Group1-group-significance.qzv \
  --p-pairwise
  
  qiime diversity beta-group-significance \
  --i-distance-matrix filtered_Fungi_ITSx_core-metrics-results/bray_curtis_distance_matrix.qza \
  --m-metadata-file orchid_metadata.tsv \
  --m-metadata-column Group2 \
  --o-visualization filtered_Fungi_ITSx_core-metrics-results/bray_curtis-Group2-group-significance.qzv \
  --p-pairwise
 
 cp filtered_Fungi_ITSx_core-metrics-results/shannon-group-significance.qzv $my_windows_folder/results_qiime2/filtered_Fungi_ITSx_core-metrics-results
 cp filtered_Fungi_ITSx_core-metrics-results/bray_curtis-*-group-significance.qzv $my_windows_folder/results_qiime2/filtered_Fungi_ITSx_core-metrics-results

#############################ANALYSIS WITH FUNGUILD########################

###Create a biom table with taxonomy annotations
qiime tools export --input-path full_Fungi_ITSx_dada2out/table.qza --output-path exported_tvs
biom convert -i exported_tvs/feature-table.biom -o exported_tvs/feature-table.tsv --to-tsv
biom head -i exported_tvs/feature-table.tsv

qiime tools export --input-path full_Fungi_ITSx_taxonomy.qza --output-path exported_tvs

cp exported_tvs/taxonomy.tsv biom-taxonomy.tsv

###Change the first line of biom-taxonomy.tsv (i.e. the header) to this:
#OTUID	taxonomy	confidence

### From:
#		Feature ID	Taxon	Confidence

sed -i -e 's/Feature ID/#OTUID/g' -e 's/Taxon/taxonomy/g' -e 's/Confidence/confidence/g' exported_tvs/biom-taxonomy.tsv
head exported_tvs/biom-taxonomy.tsv | cat -A


### try with join, but first need to remove top header
 grep -v 'biom' feature-table.tsv > feature_table_wo_biom_header.tsv
 
 ### change column orders in biom-taxonomy.tsv

awk -F $'\t' 'BEGIN {OFS = FS} {print $1,$3"|"$2}' biom-taxonomy.tsv | sed -e 's/confidence//g' -e 's/|taxonomy/taxonomy/g' > biom-taxonomy-reordered.tsv
##output of wc -l: 4415 biom-taxonomy-reordered.tsv

join -t $'\t' --header feature_table_wo_biom_header.tsv biom-taxonomy-reordered.tsv | sed -e 's/#OTU ID/OTU ID/g' > Fungi_OTU_table.txt
cp Fungi_OTU_table.txt $my_windows_folder/results_qiime2

### submitted to http://www.stbates.org/guilds/app.php
### resulting file Fungi_OTU_table.guilds.txt

##Get Trophic Mode from FunGUILD output file
awk -F $'\t' 'BEGIN {OFS = FS}{print $53}' Fungi_OTU_table.guilds.txt | sort | uniq -c > TrophicMode_count.txt

##Get Guild from FunGUILD output file
awk -F $'\t' 'BEGIN {OFS = FS}{print $54}' Fungi_OTU_table.guilds.txt | sort | uniq -c > Guild_count.txt
cp Guild_count.txt /mnt/c/Users/Biotecnologia/Documents/Projects/ITS_Orchids
cp TrophicMode_count.txt /mnt/c/Users/Biotecnologia/Documents/Projects/ITS_Orchids

###########################CREATE HEATMAPS################################################

###  To generate a heatmap containing taxonomic annotations, use `qiime taxa collapse` to collapse the feature table at the desired taxonomic level.

### collapse taxa at taxon level2
qiime taxa collapse \
  --i-table filtered-Fungi-table.qza \
  --i-taxonomy full_Fungi_ITSx_taxonomy.qza \
  --p-level 2 \
  --o-collapsed-table filtered-Fungi-table-level2.qza \
  --verbose

### export table - creates a biom artifact
qiime tools export \
  --input-path filtered-Fungi-table-level2.qza \
  --output-path exported-filtered-tables
### change name 
mv exported-filtered-tables/feature-table.biom exported-filtered-tables/filtered-Fungi-table-level2.biom
### export biom to tsv
biom convert -i exported-filtered-tables/filtered-Fungi-table-level2.biom -o filtered-Fungi-table-level2.tsv --to-tsv


### get first line with header
grep 'OTU ID' filtered-Fungi-table-level2.tsv | sed 's/\#OTU ID/Phylum/g'> filtered-Fungi-table-phylum.tsv

###get only phylum name and frequencies per sample
awk 'NR>2' filtered-Fungi-table-level2.tsv | awk -F ";" '{print $2}' | grep -v 'p__'$'\t' | grep -v 'p__unidentified' | grep -v '^__' | sed -e 's/p__//g' -e 's/\[//g' -e 's/\]/_b/g' >> filtered-Fungi-table-phylum.tsv

### this table needs to be transposed (https://www.thelinuxrain.com/articles/transposing-rows-and-columns-3-methods) to follow R heatmap tutorial. The program datamash was installed and executed.
cat filtered-Fungi-table-phylum.tsv | datamash transpose > filtered-Fungi-table-phylum.transposed.tsv

### To create a heatmap with gplots, it is necessary to install fortran library 
### conda install -c anaconda libgfortran
### This is needed to install a dependency for gplots - but still couldn't install KernSmooth a dependency of gplots
### Tried with ggplots2
### follow tutorial https://www.molecularecologist.com/2013/08/making-heatmaps-with-r-for-microbiome-analysis/
### load R libraries
R
library(ggplot2)
library(Heatplus)
library(vegan)
library(RColorBrewer)

## If needed Remove/clear environment in R 
### rm(list=ls())

## create vector to color samples according to tratment
## COFFEE_BRANCH	black
## COFFEE_TRUNK HIGH	orange
## COFFEE_TRUNK LOW	yellow
## COFFEE_TWIG	tomato
## TREES_BRANCH	mediumseagreen
## TREES_TRUNK HIGH	orchid
## TREES_TRUNK LOW	skyblue
## TREES_TWIG	blue


treatment <- c("tomato","yellow","blue","orange","mediumseagreen","orchid","skyblue","yellow","tomato","black","orange","yellow","blue","mediumseagreen","orchid","skyblue","black","tomato","black","orange","yellow","tomato","black","orange","tomato","yellow","blue","orange","mediumseagreen","orchid","skyblue","yellow","tomato","black","orange","yellow","blue","mediumseagreen","orchid","skyblue","black","tomato","black","orange","yellow","tomato","black","orange")

### colorRampPalette is in the RColorBrewer package.
### This creates a colour palette that shades from light yellow to red in RGB space with 100 unique colours
scaleyellowred <- colorRampPalette(c("lightyellow", "red"), space = "rgb")(100) 

###################### heatmaps for Phyla
### import data
phylum.data <- read.csv("filtered-Fungi-table-phylum.transposed.tsv", sep="\t", header=TRUE)
### strip off the sample ids and convert them to row names so that the data matrix contains only sequence count data.
### the header for sample names was "Phylum" due to previous table transposing
row.names(phylum.data) <- phylum.data$Phylum
phylum.data <- phylum.data[, -1]
### transform the raw counts of reads to proportions within a sample:
phylum.prop <- phylum.data/rowSums(phylum.data)
phylum.prop[1:3, 1:3]
# calculate the Bray-Curtis dissimilarity matrix on the full dataset:
phylum.dist <- vegdist(phylum.prop, method = "bray")
### Do average linkage hierarchical clustering. Other options are 'complete' or 'single'. 
### You'll need to choose the one that best fits the needs of your situation and your data.
row.clus <- hclust(phylum.dist, "aver")
###  add a column dendrogram to cluster the phyla that occur more often together.
### you have to transpose the dataset to get the genera as rows
phylum.dist.t <- vegdist(t(phylum.prop), method = "bray")
col.clus <- hclust(phylum.dist.t, "aver")
### make the heatmap with Rowv = as.dendrogram(row.clus)
### add colors to side rows for plant type

pdf("phylum_cluster_bray_treatment2.pdf")
heatmap(as.matrix(phylum.prop), Rowv = as.dendrogram(row.clus), Colv = as.dendrogram(col.clus), col = scaleyellowred, RowSideColors = treatment, margins = c(12, 3))
dev.off()

### copy to windows using another terminal
 cp phylum_cluster_bray_treatment2.pdf $my_windows_folder/results_qiime2



##################################### CLAM TEST #######################################################

### export sheet from excel file  $my_online_folder/Fungi_OTU_table.guilds_reorganizedv2
### It is one sheet with counts of sequences identified as mycorrhizal. Copy to current directory

cp $my_online_folder/mycorrhizal_counts_for_R.txt .

### Transpose it to import into R
cat mycorrhizal_counts_for_R.txt | datamash transpose > mycorrhizal_counts_for_R.transposed.tsv

### run clam test with R, library vegan
R
library(vegan)

mycorriza.data <- read.csv("mycorrhizal_counts_for_R.transposed.tsv", sep="\t", header=TRUE)
### strip off the sample ids and convert them to row names so that the data matrix contains only sequence count data.
### the header for sample names was "Symbiont" due to previous table transposing
row.names(mycorriza.data) <- mycorriza.data$Symbiont
mycorriza.data <- mycorriza.data[, -1]

## test for mycorrhiza

mycorriza_test <- clamtest(mycorriza.data, host)
summary(mycorriza_test)

pdf("mycorrhiza_specialists.pdf")
plot(mycorriza_test)
dev.off()

png("mycorrhiza_specialists.png")
plot(mycorriza_test)
dev.off()

cp mycorrhiza_specialists.p* $my_online_folder/

write.table(mycorriza_test, "mycorriza_specialists.txt", sep="\t")

### CLAM test for all ASV identified as Fungi

### export table - creates a biom artifact
qiime tools export \
  --input-path filtered-Fungi-table.qza \
  --output-path exported-filtered-tables
### change name 
mv exported-filtered-tables/feature-table.biom exported-filtered-tables/filtered-Fungi-all-levels.biom
### export biom to tsv
biom convert -i exported-filtered-tables/filtered-Fungi-all-levels.biom -o filtered-Fungi-all-levels.tsv --to-tsv

## prepare to import to R
### get first line with header
grep 'OTU ID' filtered-Fungi-all-levels.tsv | sed 's/\#OTU ID/All_Levels/g'> filtered-Fungi-table-All_Levels.tsv

###get only frequencies per sample
awk 'NR>2' filtered-Fungi-all-levels.tsv >> filtered-Fungi-table-All_Levels.tsv

### this table needs to be transposed (https://www.thelinuxrain.com/articles/transposing-rows-and-columns-3-methods) 
cat filtered-Fungi-table-All_Levels.tsv | datamash transpose > filtered-Fungi-table-All_Levels.transposed.tsv

## now run R
R
library(vegan)

all_levels_Fungi.data <- read.csv("filtered-Fungi-table-All_Levels.transposed.tsv", sep="\t", header=TRUE)
### strip off the sample ids and convert them to row names so that the data matrix contains only sequence count data.
### the header for sample names was "All_Levels" due to previous table transposing
row.names(all_levels_Fungi.data) <- all_levels_Fungi.data$All_Levels 
all_levels_Fungi.data <- all_levels_Fungi.data[, -1]

#Apply clamtest using default values:  coverage.limit = 10, specialization = 2/3, npoints = 20, alpha = 0.05/20, and using host -plant type - as grouping criterium
all_levels_Fungi_test <- clamtest(all_levels_Fungi.data, host)
summary(all_levels_Fungi_test)

pdf("all_levels_Fungi_specialists.pdf")
plot(all_levels_Fungi_test)
dev.off()

png("all_levels_Fungi_specialists.png")
plot(all_levels_Fungi_test)
dev.off()

cp all_levels_Fungi_specialists.p* $my_online_folder/

write.table(all_levels_Fungi_test, "all_levels_Fungi_specialists.txt", sep="\t")





  
