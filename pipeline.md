## Assembly Rhinolophus_cornutus and Vespertilio_sinensis

```
### Rhinolophus_cornutus 
~/trinityrnaseq-v2.12.0/Trinity --seqType fq \
--left ../rawdata/QC/S8_1.clean.fq.gz,../rawdata/QC/S10_1.clean.fq.gz,../rawdata/QC/S12_1.clean.fq.gz,../rawdata/QC/S13_1.clean.fq.gz \
--right ../rawdata/QC/S8_2.clean.fq.gz,../rawdata/QC/S10_2.clean.fq.gz,../rawdata/QC/S12_2.clean.fq.gz,../rawdata/QC/S13_2.clean.fq.gz \
--CPU 30 --min_contig_length 150 --output Rhinolophus_cornutus_trinity --max_memory 100G
### Vespertilio_sinensis
nohup ~/trinityrnaseq-v2.12.0/Trinity --seqType fq \
--left /s3_d4/caorui/batRNAseq/rawdata/QC/S7_1.clean.fq.gz,/s3_d4/caorui/batRNAseq/rawdata/QC/S9_1.clean.fq.gz,/s3_d4/caorui/batRNAseq/rawdata/QC/S11_1.clean.fq.gz,/s3_d4/caorui/batRNAseq/rawdata/QC/JP21_041aL_1.clean.fq.gz,/s3_d4/caorui/batRNAseq/rawdata/QC/JP21_059aL_1.clean.fq.gz,/s3_d4/caorui/batRNAseq/rawdata/QC/JP21_049aR_1.clean.fq.gz \
--right /s3_d4/caorui/batRNAseq/rawdata/QC/S7_2.clean.fq.gz,/s3_d4/caorui/batRNAseq/rawdata/QC/S9_2.clean.fq.gz,/s3_d4/caorui/batRNAseq/rawdata/QC/S11_2.clean.fq.gz,/s3_d4/caorui/batRNAseq/rawdata/QC/JP21_041aL_2.clean.fq.gz,/s3_d4/caorui/batRNAseq/rawdata/QC/JP21_059aL_2.clean.fq.gz,/s3_d4/caorui/batRNAseq/rawdata/QC/JP21_049aR_2.clean.fq.gz  \
--CPU 200 --min_contig_length 150 --output Vespertilio_sinensis_trinity --max_memory 100G
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
################################
## Counts of transcripts, etc.
################################
Total trinity 'genes':	278882
Total trinity transcripts:	395383
Percent GC: 50.14

########################################
Stats based on ALL transcript contigs:
########################################

	Contig N10: 5731
	Contig N20: 4213
	Contig N30: 3230
	Contig N40: 2505
	Contig N50: 1881

	Median contig length: 324
	Average contig: 803.67
	Total assembled bases: 317756635


#####################################################
## Stats based on ONLY LONGEST ISOFORM per 'GENE':
#####################################################

	Contig N10: 4850
	Contig N20: 3158
	Contig N30: 2031
	Contig N40: 1177
	Contig N50: 683

	Median contig length: 266
	Average contig: 495.38
	Total assembled bases: 138153436

```

(3) use cd-hit to remove redundance, set  identity threshold of 95%
```
cd-hit-est -i Trinity.fasta -o CD-hit-RC_Trinity.fasta -c 0.95 
cd-hit-est -i Trinity.fasta -o CD-hit-VS_Trinity.fasta -c 0.95 
```
(4) BUSCO to evaulate assembly quality 
```
nohup busco -i Trinity.fasta -m tran -f --offline -l ~/database/laurasiatheria_odb10 -o RC_busco -c 15
```
(5) Bowtie2 mapping 
```
bowtie2 --local --no-unal -x Trinity.fasta -q -1 /s3_d4/caorui/batRNAseq/rawdata/QC/S20_1.clean.fq.gz -2 /s3_d4/caorui/batRNAseq/rawdata/QC/S20_2.clean.fq.gz | samtools view -Sb - | samtools sort -o bowtie2.coordSorted.bam
------------------------------------------------------------------------------------------------------------------------------
Result
22015425 (100.00%) were paired; of these:
    715016 (3.25%) aligned concordantly 0 times
    1995997 (9.07%) aligned concordantly exactly 1 time
    19304412 (87.69%) aligned concordantly >1 times
    ----
    715016 pairs aligned concordantly 0 times; of these:
      14608 (2.04%) aligned discordantly 1 time
    ----
    700408 pairs aligned 0 times concordantly or discordantly; of these:
      1400816 mates make up the pairs; of these:
        394376 (28.15%) aligned 0 times
        104113 (7.43%) aligned exactly 1 time
        902327 (64.41%) aligned >1 times
99.10% overall alignment rate
```

## Annotation
### 1. homology search
#### (1) Identification of likely protein-coding regions in transcripts
```
mv ~/s3_d4/caorui/batRNAseq/reference
mkdir && cd Trinotate
TransDecoder.LongOrfs -t ../../Rhinolophus_cornutus_trinity/c-0.95_RC-trintiy.fasta
TransDecoder.Predict -t ../../Rhinolophus_cornutus_trinity/c-0.95_RC-trintiy.fasta
```
#### (2) Sequence homology searches
```
diamond blastx -q ../../Rhinolophus_cornutus_trinity/c-0.95_RC-trintiy.fasta -d uniprot_sprot.pep -p 100 --max-target-seqs 1 --outfmt 6 --evalue 1e-5 -o blastx.outfmrt
diamond blastp -q c-0.95_RC-trintiy.fasta.transdecoder.pep -d uniprot_sprot.pep -p 100 --max-target-seqs 1 --outfmt 6 --evalue 1e-5 -o transdecoder_blastp.outfmrt
```
#### (3) HMMER search against the Pfam database, and identify conserved domains 

```
hmmpress Pfam-A.hmm 
hmmscan --cpu 100 --domtblout TrinotatePFAM.out Pfam-A.hmm c-0.95_RC-trintiy.fasta.transdecoder.pep >pfam.log
```
#### (4) input to the sqlite
```
~/trinityrnaseq-v2.12.0/util/support_scripts/get_Trinity_gene_to_trans_map.pl c-0.95_RC-trintiy.fasta >c-0.95_RC-trintiy.fasta.gene_trans_map ### prepare transcript map
## inital the sqlite 
~/Trinotate-Trinotate-v3.2.2/Trinotate Trinotate.sqlite init \
--gene_trans_map ../../Rhinolophus_cornutus_trinity/c-0.95_RC-trintiy.fasta.gene_trans_map \
--transcript_fasta ../../Rhinolophus_cornutus_trinity/c-0.95_RC-trintiy.fasta \
--transdecoder_pep c-0.95_RC-trintiy.fasta.transdecoder.pep

~/Trinotate-Trinotate-v3.2.2/Trinotate Trinotate.sqlite LOAD_pfam TrinotatePFAM.out ### load PFam
~/Trinotate-Trinotate-v3.2.2/Trinotate Trinotate.sqlite  LOAD_swissprot_blastx blastx.outfmrt ### load blastx
~/Trinotate-Trinotate-v3.2.2/Trinotate Trinotate.sqlite  LOAD_swissprot_blastp transdecoder_blastp.outfmrt ### load blastp
```
#### (5) annotated by interprotein
```
nohup interproscan.sh -i ../../Rhinolophus_cornutus_trinity/c-0.95_RC-trintiy.fasta -appl Pfam -f GFF3 -goterms -cpu 100 -dp
```
#### (6) annotated by eggmapper
```
nohup ~/eggnog-mapper/emapper.py -i CD-hit-VS_Trinity.fasta.transdecoder.pep --output VS_eggnog -m diamond --cpu 30 --override
```

### 2. orthology assignment (use mouse_pep.faa)
```
formatdb -i Mus_musculus.GRCm39.pep.all.fa -p T
# -i specifies the input file
# -p lets the program know whether the input is protein (T) or nucleotide (F)
blastp -query /s3_d4/caorui/batRNAseq/reference/annotation/RC/c-0.95_RC-trintiy.fasta.transdecoder.pep -db Mus_musculus.GRCm39.pep.all.fa -num_threads 80 -max_target_seqs 1 -outfmt 6 -evalue 1e-5 > blastp_RC_Mus.outfmt6 
formatdb -i /s3_d4/caorui/batRNAseq/reference/annotation/RC/c-0.95_RC-trintiy.fasta.transdecoder.pep -p T ## second blastp
blastp -query Mus_musculus.GRCm39.pep.all.fa -db /s3_d4/caorui/batRNAseq/reference/annotation/RC/c-0.95_RC-trintiy.fasta.transdecoder.pep -num_threads 80 -max_target_seqs 1 -outfmt 6 -evalue 1e-5 > blastp_Mus_RC.outfmt6
```
## Quantify by RSEM and Bowtie2
```
for sample in `awk '{print $1}' ./sample.txt`; do echo $sample\
~/trinityrnaseq-v2.12.0/util/align_and_estimate_abundance.pl --seqType fq \
--left /s3_d4/caorui/batRNAseq/rawdata/QC/"$sample"_1.clean.fq.gz \
--right /s3_d4/caorui/batRNAseq/rawdata/QC/"$sample"_2.clean.fq.gz \
--transcripts /s3_d4/caorui/batRNAseq/reference/Rhinolophus_cornutus_trinity/RS_unigene.fasta \
--est_method RSEM  --aln_method bowtie2 --trinity_mode --prep_reference --output_dir "$sample".RSEM; done
```
## combine the expression matrix
```
~/trinityrnaseq-v2.12.0/util/abundance_estimates_to_matrix.pl --est_method RSEM --out_prefix RS_Trinity_trans *genes*  --gene_trans_map ../../reference/Rhinolophus_cornutus_trinity/Trinity.fasta.gene_trans_map
```
