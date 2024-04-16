# Translational-regulation-of-TDMD

## 1 TDMD identification analysis
### 1.1 Cutadapt （https://cutadapt.readthedocs.io/en/stable/） in FASTQ
cutadapt -a TGGAATTCTCGGGTGCCAAG -A GATCGTCGGACTGTAGAACT -o `test_R1_cut.fastq` -p `test_R2_cut.fastq` `test_R1.fastq test_R2.fastq` --minimum-length 18 -j 10 -m 26

### 1.2 Pear (https://cme.h-its.org/exelixis/web/software/pear/doc.html) in FASTQ
pear -f `test_R1_cut.fastq` -r `test_R2_cut.fastq` -j 10 -o `test` -n 26

### 1.3 Collapse (http://hannonlab.cshl.edu/fastx_toolkit/) PCR duplicated reads
fastx_collapser  -i `test.assembled.fastq` -o `test_collapsed.fasta`

### 1.4 Remove UMI sequences from 5' and 3' of FASTA
cutadapt -u 4 -u -4 -m 18 `test_collapsed.fasta` -o `test_cutN.fasta` -j 10              

### 1.5 Run hyb(https://github.com/gkudla/hyb) software
module load hyb

module load unafold

hyb analyse in=`test_cutN.fasta` db=`martquery_0228172604_181_human_unique1_add18S28S` type=mim pref=mim format=fasta

### 1.6 Calculate hybrids reads
python3 CLASH.py hybrid_number_statistic -d `martquery_0228172604_181_human_unique1_add18S28S.fasta` -i `23090FL-01-01-01_S2.cutUMI_comp_martquery_0228172604_181_human_unique1_add18S28S_hybrids_ua.hyb`c

### 1.7 Potential TDMD miRNA-target RNA hybrids identification
python3 CLASH.py Viennad_to_Table  -c `human_transcripts_CS_20220422.txt` -t `martquery_0228172604_181_human_unique1_add18S28S.fasta` -n `martquery_0228172604_181_name.txt` -i `23090FL-01-01-01_S2.cutUMI_comp_martquery_0228172604_181_human_unique1_add18S28S_hybrids_ua`

```ruby
martquery_0228172604_181_human_unique1_add18S28S.fasta file size is ???MB. The dataset folder has a small size file, called ???.
  
human_transcripts_CS_20220422.txt file size is ???G. The dataset folder has a small size file, called ???.
```


## 2 miRNA analysis
### 2.1 Deduplicate clean-miRNA-seq reads from BGI
  python3 CLASH.py deduplicate_BGI -i `test_miRNA_BGI.fq`
### 2.2 miRNA abundance calculation
python3 CLASH.py miRNA_abundance -i `test_miRNA_BGI_deduplicated.fa` -d `20220221_dm6_miRNA_database.fasta`

### 2.3 Differential expression level analysis
run `Deseq.R` code

### 2.4 miRNA length distribution (isoform) count
python3 CLASH.py miRNA_length_distribution -i `test_miRNA_BGI_deduplicated.fa` -d `20220221_dm6_miRNA_database.fasta` -m `[--miRNA_sequence]`

e.g. miR-999 length distribution count

*python3 CLASH.py miRNA_length_distribution -i `test_miRNA_BGI_deduplicated.fa` -d `20220221_dm6_miRNA_database.fasta` -m TGTTAACTGTAAGACTGTGTCT*

