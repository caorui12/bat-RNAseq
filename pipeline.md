## Assembly Rhinolophus_cornutus and Vespertilio_sinensis

```
nohup ~/trinityrnaseq-v2.12.0/Trinity --seqType fq \
--left ../rawdata/QC/S8_1.clean.fq.gz,../rawdata/QC/S10_1.clean.fq.gz,../rawdata/QC/S12_1.clean.fq.gz,../rawdata/QC/S13_1.clean.fq.gz \
--right ../rawdata/QC/S8_2.clean.fq.gz,../rawdata/QC/S10_2.clean.fq.gz,../rawdata/QC/S12_2.clean.fq.gz,../rawdata/QC/S13_2.clean.fq.gz \
--CPU 30 --min_contig_length 150 --output Rhinolophus_cornutus_trinity --max_memory 100G
```

## Check assembly quality
(1) N50 for Rhinolophus_cornutus
```
~/trinityrnaseq-v2.12.0/util/TrinityStats.pl Trinity.fasta
################################
## Counts of transcripts, etc.
################################
Total trinity 'genes':	333798
Total trinity transcripts:	421903
Percent GC: 46.56

########################################
Stats based on ALL transcript contigs:
########################################

	Contig N10: 7095
	Contig N20: 5345
	Contig N30: 4229
	Contig N40: 3360
	Contig N50: 2555

	Median contig length: 327
	Average contig: 895.42
	Total assembled bases: 377778420


#####################################################
## Stats based on ONLY LONGEST ISOFORM per 'GENE':
#####################################################

	Contig N10: 6121
	Contig N20: 4179
	Contig N30: 2715
	Contig N40: 1504
	Contig N50: 838

	Median contig length: 277
	Average contig: 546.37
	Total assembled bases: 182377547
```
(2) N50 for Vespertilio_sinensis

```
Total trinity 'genes':	197255
Total trinity transcripts:	286637
Percent GC: 49.55

########################################
Stats based on ALL transcript contigs:
########################################

	Contig N10: 5277
	Contig N20: 3992
	Contig N30: 3130
	Contig N40: 2476
	Contig N50: 1915

	Median contig length: 351
	Average contig: 851.99
	Total assembled bases: 244211681


#####################################################
## Stats based on ONLY LONGEST ISOFORM per 'GENE':
#####################################################

	Contig N10: 4689
	Contig N20: 3250
	Contig N30: 2264
	Contig N40: 1457
	Contig N50: 847

	Median contig length: 273
	Average contig: 534.70
	Total assembled bases: 105471994
```

(3) remove redundant 
```
~/trinityrnaseq-v2.12.0/util/misc/get_longest_isoform_seq_per_trinity_gene.pl Trinity.fasta >VS_unigene.fasta
~/trinityrnaseq-v2.12.0/util/misc/get_longest_isoform_seq_per_trinity_gene.pl Trinity.fasta >RS_unigene.fasta
```
(4) index and alignment by bowtiew
```
bowtie2-build RS_unigene.fasta RS_unigene.fasta
bowtie2-build VS_unigene.fasta VS_unigene.fasta
```
## Annotation
### (1) Identification of likely protein-coding regions in transcripts
```
mv ~/s3_d4/caorui/batRNAseq/reference
mkdir && cd Trinotate
ransDecoder.LongOrfs -t ../Rhinolophus_cornutus_trinity/RS_unigene.fasta
TransDecoder.Predict -t ../Rhinolophus_cornutus_trinity/RS_unigene.fasta
```
### (2) Sequence homology searches
```
~/diamond blastx -p 10 -q ../Rhinolophus_cornutus_trinity/RS_unigene.fasta -d /s3_d3/xiuwan/reference_genome/nr.dmnd.tax.v2.0.dmnd  --max-target-seqs 1 --outfmt 6 --evalue 1e-5 -o RS_unigene.fasta.transdecoder.outfmt6
```
### (3) HMMER search against the Pfam database, and identify conserved domains 

```
nohup hmmscan --cpu 20 --domtblout TrinotateRS_PFAM.out /s3_d4/caorui/database/pfam/Pfam-A.hmm RS_unigene.fasta.transdecoder.pep
```
