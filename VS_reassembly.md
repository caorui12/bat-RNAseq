## use Trinity to assembly the Vespertilio_sinensis, each stages select one sample
nohup ~/trinityrnaseq-v2.12.0/Trinity --seqType fq \
--left /s3_d4/caorui/batRNAseq/rawdata/QC/JP21_041aL_1.clean.fq.gz,/s3_d4/caorui/batRNAseq/rawdata/QC/JP21_037aL_1.clean.fq.gz,/s3_d4/caorui/batRNAseq/rawdata/QC/S27_1.clean.fq.gz,/s3_d4/caorui/batRNAseq/rawdata/QC/S35_1.clean.fq.gz,/s3_d4/caorui/batRNAseq/rawdata/QC/S15_1.clean.fq.gz,/s3_d4/caorui/batRNAseq/rawdata/QC/S28_1.clean.fq.gz,/s3_d4/caorui/batRNAseq/rawdata/QC/S14_1.clean.fq.gz,/s3_d4/caorui/batRNAseq/rawdata/QC/S30_1.clean.fq.gz \
--right /s3_d4/caorui/batRNAseq/rawdata/QC/JP21_041aL_2.clean.fq.gz,/s3_d4/caorui/batRNAseq/rawdata/QC/JP21_037aL_2.clean.fq.gz,/s3_d4/caorui/batRNAseq/rawdata/QC/S27_2.clean.fq.gz,/s3_d4/caorui/batRNAseq/rawdata/QC/S35_2.clean.fq.gz,/s3_d4/caorui/batRNAseq/rawdata/QC/S15_2.clean.fq.gz,/s3_d4/caorui/batRNAseq/rawdata/QC/S28_2.clean.fq.gz,/s3_d4/caorui/batRNAseq/rawdata/QC/S14_2.clean.fq.gz,/s3_d4/caorui/batRNAseq/rawdata/QC/S30_2.clean.fq.gz \
--CPU 30 --min_contig_length 150 --output Vespertilio_sinensis_trinity --max_memory 100G &