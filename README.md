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

hyb analyse in=`test_cutN.fasta` db=`human_transcript_database` type=mim pref=mim format=fasta

### 1.6 Calculate hybrids reads
python3 CLASH.py hybrid_number_statistic -d `human_transcript_database.fasta` -i `output_ua.hyb`

### 1.7 Potential TDMD miRNA-target RNA hybrids identification
python3 CLASH.py Viennad_to_Table  -c `human_transcripts_CS_20220422.txt` -t `human_transcript.fasta` -n `human_transcript_name.txt` -i `output_hybrids_ua`

## 2 miRNA analysis

### 2.1 Cutadapt （https://cutadapt.readthedocs.io/en/stable/） in FASTQ
cutadapt -a TGGAATTCTCGGGTGCCAAG -A GATCGTCGGACTGTAGAACT -o `test_R1_cut.fastq` -p `test_R2_cut.fastq` `test_R1.fastq test_R2.fastq` --minimum-length 18 -j 10 -m 26

### 2.2 Pear (https://cme.h-its.org/exelixis/web/software/pear/doc.html) in FASTQ
pear -f `test_R1_cut.fastq` -r `test_R2_cut.fastq` -j 10 -o `test` -n 26

### 2.3 Collapse (http://hannonlab.cshl.edu/fastx_toolkit/) PCR duplicated reads
fastx_collapser  -i `test.assembled.fastq` -o `test_collapsed.fasta`

### 2.4 Remove UMI sequences from 5' and 3' of FASTA
cutadapt -u 4 -u -4 -m 18 `test_collapsed.fasta` -o `test_cutN.fasta` -j 10 
  
### 2.5 miRNA abundance calculation
python3 CLASH.py all_miRNA_isoform_table_2nd_18th -d `Exp326_Homo_mature_HSUR4_SpikeIn_miRNA.fa` -i `test.UMI.fasta`

### 2.6 Differential expression level analysis
run `EdgeR.R` code

### 2.4 miRNA length distribution (isoform) count
python3 isoform_length.py
