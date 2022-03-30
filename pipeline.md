## Assembly Rhinolophus_cornutus and Vespertilio_sinensis

```
nohup ~/trinityrnaseq-v2.12.0/Trinity --seqType fq \
--left ../rawdata/QC/S8_1.clean.fq.gz,../rawdata/QC/S10_1.clean.fq.gz,../rawdata/QC/S12_1.clean.fq.gz,../rawdata/QC/S13_1.clean.fq.gz \
--right ../rawdata/QC/S8_2.clean.fq.gz,../rawdata/QC/S10_2.clean.fq.gz,../rawdata/QC/S12_2.clean.fq.gz,../rawdata/QC/S13_2.clean.fq.gz \
--CPU 30 --min_contig_length 150 --output Rhinolophus_cornutus_trinity --max_memory 100G
```
