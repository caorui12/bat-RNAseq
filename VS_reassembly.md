##  Trinity assembly of Vespertilio_sinensis, each stages select one sample
```
nohup ~/trinityrnaseq-v2.12.0/Trinity --seqType fq \
--left /s3_d4/caorui/batRNAseq/rawdata/QC/JP21_041aL_1.clean.fq.gz,/s3_d4/caorui/batRNAseq/rawdata/QC/JP21_037aL_1.clean.fq.gz,/s3_d4/caorui/batRNAseq/rawdata/QC/S27_1.clean.fq.gz,/s3_d4/caorui/batRNAseq/rawdata/QC/S35_1.clean.fq.gz,/s3_d4/caorui/batRNAseq/rawdata/QC/S15_1.clean.fq.gz,/s3_d4/caorui/batRNAseq/rawdata/QC/S28_1.clean.fq.gz,/s3_d4/caorui/batRNAseq/rawdata/QC/S14_1.clean.fq.gz,/s3_d4/caorui/batRNAseq/rawdata/QC/S30_1.clean.fq.gz \
--right /s3_d4/caorui/batRNAseq/rawdata/QC/JP21_041aL_2.clean.fq.gz,/s3_d4/caorui/batRNAseq/rawdata/QC/JP21_037aL_2.clean.fq.gz,/s3_d4/caorui/batRNAseq/rawdata/QC/S27_2.clean.fq.gz,/s3_d4/caorui/batRNAseq/rawdata/QC/S35_2.clean.fq.gz,/s3_d4/caorui/batRNAseq/rawdata/QC/S15_2.clean.fq.gz,/s3_d4/caorui/batRNAseq/rawdata/QC/S28_2.clean.fq.gz,/s3_d4/caorui/batRNAseq/rawdata/QC/S14_2.clean.fq.gz,/s3_d4/caorui/batRNAseq/rawdata/QC/S30_2.clean.fq.gz \
--CPU 30 --min_contig_length 150 --output Vespertilio_sinensis_trinity --max_memory 100G &
```
## Identify ORF
```
TransDecoder.LongOrfs -t trintiy.fasta
TransDecoder.Predict -t trintiy.fasta
```
## map ansembly transcript to Myotis lucifugus (little brown bat)
(1) download protein sequence 
```
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/147/115/GCF_000147115.1_Myoluc2.0/GCF_000147115.1_Myoluc2.0_protein.faa.gz 
```
(2) reciporal blast hit
```
formatdb -i GCF_000147115.1_Myoluc2.0_protein.faa -p T
```
## eggnog-mapper annotation 
### eggNOG-mapper  is a tool for functional annotation of large sets of sequences based on fast orthology assignments
```
/home/caorui/eggnog-mapper/emapper.py -i Trinity.fasta.transdecoder.pep --output VS_eggnog -m diamond --cpu 30 --override
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
