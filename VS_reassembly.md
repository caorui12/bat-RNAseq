##  Trinity assembly of Vespertilio_sinensis, each stages select one sample
```
nohup ~/trinityrnaseq-v2.12.0/Trinity --seqType fq \
--left /s3_d4/caorui/batRNAseq/rawdata/QC/JP21_041aL_1.clean.fq.gz,/s3_d4/caorui/batRNAseq/rawdata/QC/JP21_037aL_1.clean.fq.gz,/s3_d4/caorui/batRNAseq/rawdata/QC/S27_1.clean.fq.gz,/s3_d4/caorui/batRNAseq/rawdata/QC/S35_1.clean.fq.gz,/s3_d4/caorui/batRNAseq/rawdata/QC/S15_1.clean.fq.gz,/s3_d4/caorui/batRNAseq/rawdata/QC/S28_1.clean.fq.gz,/s3_d4/caorui/batRNAseq/rawdata/QC/S14_1.clean.fq.gz,/s3_d4/caorui/batRNAseq/rawdata/QC/S30_1.clean.fq.gz \
--right /s3_d4/caorui/batRNAseq/rawdata/QC/JP21_041aL_2.clean.fq.gz,/s3_d4/caorui/batRNAseq/rawdata/QC/JP21_037aL_2.clean.fq.gz,/s3_d4/caorui/batRNAseq/rawdata/QC/S27_2.clean.fq.gz,/s3_d4/caorui/batRNAseq/rawdata/QC/S35_2.clean.fq.gz,/s3_d4/caorui/batRNAseq/rawdata/QC/S15_2.clean.fq.gz,/s3_d4/caorui/batRNAseq/rawdata/QC/S28_2.clean.fq.gz,/s3_d4/caorui/batRNAseq/rawdata/QC/S14_2.clean.fq.gz,/s3_d4/caorui/batRNAseq/rawdata/QC/S30_2.clean.fq.gz \
--CPU 30 --min_contig_length 150 --output Vespertilio_sinensis_trinity --max_memory 100G &
```
## filter transcript based on length and cd-hit
```
 ~/trinityrnaseq-v2.12.0/util/misc/get_longest_isoform_seq_per_trinity_gene.pl Trinity.fasta >longest_isoform_Trinity.fasta
~/cd-hit-v4.8.1-2019-0228/cd-hit-est -i longest_isoform_Trinity.fasta -c 0.95 -o cd-hit-c_0.95_longest_isofrom_Trinity.fasta
```
final retain 402385 transcript
## Identify ORF
```
TransDecoder.LongOrfs -t cd-hit-c_0.95_longest_isofrom_Trinity.fasta
TransDecoder.Predict -t cd-hit-c_0.95_longest_isofrom_Trinity.fasta
```
finally, 47278 transcript has potential ORF
## mulitple species annotation (MSA)
### We obtain six bat genome from 2020 nature study, identify ortholog genes for mapping. Genomes downloaded from here https://bds.mpi-cbg.de/hillerlab/Bat1KPilotProject/ 
(1) using gffread to extract cds peptide. 
(2) use in-house script to add species name and revise gene ID to symbol based on gff files
(3) use in-house script to extract longest isoform
Here show one example
#### (1) prepare files
```
gffread HLmolMol2.Bat1Kannotation.gff3 -g HLmolMol2.fa -y HLmolMol2_protein.fa
python extract_gene_ID.py HLmolMol2.Bat1Kannotation.gff3 HLmolMol2_protein.fa HLmolMol2_protein_revised.fa 
python extract_longest_isoform.py HLmolMol2_revised_protein.fa
```
#### (2) Index the file
```
formatdb -i HLmyoMyo6longest.txt -p T
formatdb -i cd-hit-c_0.95_longest_isofrom_Trinity.fasta.transdecoder.pep -p T
```
#### (3) reciporal blast
```
blastp -query cd-hit-c_0.95_longest_isofrom_Trinity.fasta.transdecoder.pep -db HLmyoMyo6longest.txt -num_threads 80 -max_target_seqs 1 -outfmt 6 -evalue 1e-5 > blastp_VS_Myo.outfmt6  
blastp -db cd-hit-c_0.95_longest_isofrom_Trinity.fasta.transdecoder.pep -query HLmyoMyo6longest.txt -num_threads 80 -max_target_seqs 1 -outfmt 6 -evalue 1e-5 > blastp_Myo_VS.outfmt6  
```
#### (4) extract one-to-one ortholg gene using in house script
```
python /s3_d4/caorui/batRNAseq/reference/ortholog/homo_gene_Extract.py blastp_VS_Myo.outfmt6 blastp_Myo_VS.outfmt6 > VS_Myo
```
#### (5) finished the rest of blastp in bash script
```
#!/usr/bin/bash

VS=cd-hit-c_0.95_longest_isofrom_Trinity.fasta.transdecoder.pep
pip=HLpipKuh2longest.txt
rhi=HLrhiFer5longest.txt
rou=HLrouAeg4longest.txt

blastp -query $VS -db $pip -num_threads 80 -max_target_seqs 1 -outfmt 6 -evalue 1e-5 > blastp_vs_pip.outfmt6
blastp -db $VS -query $pip -num_threads 80 -max_target_seqs 1 -outfmt 6 -evalue 1e-5 > blastp_pip_vs.outfmt6
python /s3_d4/caorui/batRNAseq/reference/ortholog/homo_gene_Extract.py blastp_vs_pip.outfmt6 blastp_pip_vs.outfmt6 > vs_pip_orthlog.txt
blastp -query $VS -db $rhi -num_threads 80 -max_target_seqs 1 -outfmt 6 -evalue 1e-5 > blastp_vs_rhi.outfmt6
blastp -db $VS -query $rhi -num_threads 80 -max_target_seqs 1 -outfmt 6 -evalue 1e-5 > blastp_rhi_vs.outfmt6
python /s3_d4/caorui/batRNAseq/reference/ortholog/homo_gene_Extract.py blastp_vs_rhi.outfmt6 blastp_rhi_vs.outfmt6 > vs_rhi_orthlog.txt
blastp -query $VS -db $rou -num_threads 80 -max_target_seqs 1 -outfmt 6 -evalue 1e-5 > blastp_vs_rou.outfmt6
blastp -db $VS -query $rou -num_threads 80 -max_target_seqs 1 -outfmt 6 -evalue 1e-5 > blastp_rou_vs.outfmt6
python /s3_d4/caorui/batRNAseq/reference/ortholog/homo_gene_Extract.py blastp_vs_rou.outfmt6 blastp_rou_vs.outfmt6 > vs_rou_orthlog.txt
```
#### (6) merge result
## eggnog-mapper annotation 
### eggNOG-mapper  is a tool for functional annotation of large sets of sequences based on fast orthology assignments
```
/home/caorui/eggnog-mapper/emapper.py -i cd-hit-c_0.95_longest_isofrom_Trinity.fasta.transdecoder.pep --output VS_eggnog -m diamond --cpu 30 --override
```
### Transcript factor annotation 
```
formatdb -i Myotis_lucifugus_TF_protein.fasta -p T
blastp -query ../annotation/Trinity.fasta.transdecoder.pep -db Myotis_lucifugus_TF_protein.fasta -num_threads 80 -max_target_seqs 1 -outfmt 6 -evalue 1e-5 > blastp_VS_TF.outfmt6 
formatdb -i ../annotation/Trinity.fasta.transdecoder.pep -p T
blastp -query Myotis_lucifugus_TF_protein.fasta -db ../annotation/Trinity.fasta.transdecoder.pep -num_threads 80 -max_target_seqs 1 -outfmt 6 -evalue 1e-5  > blastp_TF_VS.outfmt6 
```
### Then use in-house script to extract reciporal blast result

```
python ../../ortholog/homo_gene_Extract.py blastp_TF_VS.outfmt6 blastp_VS_TF.outfmt6 >TF_VS.RBH
```
