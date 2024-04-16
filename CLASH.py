import csv, getopt, collections, os, re, sys, glob,subprocess
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
from Bio.Seq import Seq
import plotly.express as px
from scipy.stats import mannwhitneyu

class BedGraph():
    def making_simple_version(self, file, chromosome_name):
        with open(file, 'r+') as f1:
            for line1 in f1:
                line2 = line1.strip().split('\t')
                chromosome = line2[0]
                if chromosome_name == chromosome:
                    chr_start = int(line2[1])
                    chr_end = int(line2[2])
                    with open(f"{file.replace('.bedGraph', '')}_{chromosome}.bedGraph", 'a+') as f2:
                        if chr_end - chr_start == 1:
                            # f2.write(chromosome.replace('chrM', 'chrMT') + '\t' + str(chr_end) + '\t' + line2[3] + '\n')

                            f2.write(chromosome.replace('chrM', 'chrMtDNA') + '\t' + str(chr_end) + '\t' + line2[3] + '\n')
                        else:
                            for position in list(range(chr_start + 1, chr_end + 1)):
                                f2.write(chromosome + '\t' + str(position) + '\t' + line2[3] + '\n')

    def no_exist_conservation_socre_in_bedGraph(self, chr, position,
                                                f_transcript_conservation_score):  # I labeled '0' for individuals who do not have a conservation score in their genome.
        mux = pd.MultiIndex.from_arrays([[chr], [position]], names=[0, 1])
        df_if_conservation_score_not_exist = pd.DataFrame(0, index=mux, columns=[2])
        f_transcript_conservation_score = pd.concat(
            [f_transcript_conservation_score, df_if_conservation_score_not_exist])
        return f_transcript_conservation_score


class Database():
    def microRNA_database(self, input):  ## For analyzing the length of each miRNA\n",
        with open(input, 'r+') as f1:
            dict_miRNA, name1, seq1 = {}, '', ''
            for line1 in f1:
                if line1[0] == '>':
                    name1 = line1.strip()[1:]
                else:
                    seq1 = line1.strip()
                    if ('microRNA' in name1) and (seq1 != ''):
                        dict_miRNA[name1.split('_')[2]] = seq1
        return dict_miRNA

    def microRNA_sequence_to_name_database_1st_18nt(self,
                                                    input):  ## dictionary miRNA seq to name, miRNA length at least 18nt
        with open(input, 'r+') as f1:
            dict_miRNA, name1, seq1 = {}, '', ''
            for line1 in f1:
                if line1[0] == '>':
                    name1 = line1.strip().split(' ')[0][1:]
                else:
                    seq1_short = line1.strip()[:18]
                    seq1_raw = line1.strip()
                    name1_seq1raw = f"{name1}&{seq1_raw}"
                    if (('miR-' in name1) or ('let-' in name1)) and (seq1_short != ''):
                        dict_miRNA[seq1_short] = dict_miRNA.get(seq1_short, '')
                        dict_miRNA[seq1_short] = (dict_miRNA[seq1_short] + '_' + name1_seq1raw).lstrip('_')
        return dict_miRNA

    def microRNA_sequence_to_name_database_whole_nt(self,
                                                    input):  ## dictionary miRNA seq to name, miRNA length at least 18nt
        with open(input, 'r+') as f1:
            dict_miRNA, name1, seq1 = {}, '', ''
            for line1 in f1:
                if line1[0] == '>':
                    name1 = line1.strip().split(' ')[0].replace('>hsa-','')
                else:
                    seq1_raw = line1.strip()
                    if (('miR-' in name1) or ('let-' in name1)) and (seq1_raw != ''):
                        dict_miRNA[seq1_raw] = name1
        return dict_miRNA

    def microRNA_sequence_to_name_database_16nt(self,
                                                input):  ## dictionary miRNA seq to name, miRNA length at least 16nt
        with open(input, 'r+') as f1:
            dict_miRNA, name1, seq1 = {}, '', ''
            for line1 in f1:
                if line1[0] == '>':
                    name1 = line1.strip().split(' ')[0][1:]
                else:
                    seq1 = line1.strip()[:16]
                    if (('miR-' in name1) or ('let-' in name1)) and (seq1 != ''):
                        dict_miRNA[seq1] = dict_miRNA.get(seq1, '')
                        dict_miRNA[seq1] = (dict_miRNA[seq1] + '_' + name1).lstrip('_')
        return dict_miRNA

    def making_unique_redundant_database_hg38(self, input):
        normal_chromosome_list = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16',
                                  '17', '18', '19', '20', '21', '22', 'X', 'Y', 'MT']
        unchr_name = {'RNA5-8SN4', 'KIR3DS1', 'KIR2DS1', 'FAM8A5P', 'RNU1-79P', 'KIR2DS5', 'LILRA3', 'HLA-DRB8',
                      'OR8U9', 'KIR2DS3', 'MAFIP', 'CCL3L1', 'TAS2R45', 'PRSS3P2', 'CACNA1C-IT2', 'KIR2DL2', 'KIR2DS2',
                      'HLA-DRB7', 'KIR2DL5A', 'KIR2DL5B', 'RNA5-8SN5', 'C4B_2', 'HLA-DRB3', 'OR8U8', 'GTF2H2C_2',
                      'HLA-DRB2', 'HLA-DRB4', 'OR9G9', 'PRAMEF22', 'GSTT1', 'RNU1-116P'}
        dict_set = set()
        dict_gene_type = set()
        f2 = open(input.replace('.txt', '_human_unique.fasta'), 'w+')
        f3 = open(input.replace('.txt', '_human_redundant.fasta'), 'w+')
        # integrated_name = f'>{Gene_ID}_{Transcript_ID}_{Gene_name}#{Chromosome}#{exons_left_list}#{exons_right_list}#{strand}#{CDS_left}#{CDS_right}#{Gene_type}_mRNA'
        with open(input, 'r+') as f1:
            integrated_name, sequence = '', ''
            for line1 in f1:
                if ('>' in line1):
                    if ('miRNA' not in integrated_name) and (integrated_name != ''):
                        tmp_chr = integrated_name.split('_')[2].split('#')[1]
                        if tmp_chr in normal_chromosome_list:
                            if sequence not in dict_set:
                                f2.write(integrated_name + '\n')
                                f2.write(sequence + '\n')
                                dict_set.add(sequence)
                            elif sequence in dict_set:
                                f3.write(integrated_name + '\n')
                                f3.write(sequence + '\n')
                        elif tmp_chr not in normal_chromosome_list:
                            if sequence not in dict_set:
                                f2.write(integrated_name + '\n')
                                f2.write(sequence + '\n')
                                dict_set.add(sequence)
                            elif sequence in dict_set:
                                f3.write(integrated_name + '\n')
                                f3.write(sequence + '\n')
                    line2 = line1.replace('_', '&').strip().split('|')
                    sequence = ''
                    Gene_ID = line2[0][1:]
                    Transcript_ID = line2[1]
                    Gene_type = line2[2].replace('protein&coding', 'mRNA')
                    exons_left_list = line2[3]
                    exons_right_list = line2[4]
                    try:
                        CDS_left = min(int(x) for x in re.findall('\d+', line2[5]))
                        CDS_right = max(int(x) for x in re.findall('\d+', line2[6]))
                    except:
                        CDS_left, CDS_right = '0', '0'
                    strand = line2[7]
                    Chromosome = line2[8]
                    if len(line2) != 10:
                        Gene_name = 'NoName'
                    else:
                        Gene_name = line2[9]
                    dict_gene_type.add(Gene_type)
                    integrated_name = f'>{Gene_ID}_{Transcript_ID}_{Gene_name}#{Chromosome}#{exons_left_list}#{exons_right_list}#{strand}#{CDS_left}#{CDS_right}#{Gene_type}_mRNA'
                else:
                    sequence += line1.strip()
            if 'miRNA' not in integrated_name:
                tmp_chr = integrated_name.split('_')[2].split('#')[1]
                if tmp_chr in normal_chromosome_list:
                    if sequence not in dict_set:
                        f2.write(integrated_name + '\n')
                        f2.write(sequence + '\n')
                        dict_set.add(sequence)
                    elif sequence in dict_set:
                        f3.write(integrated_name + '\n')
                        f3.write(sequence + '\n')
        f2.close()
        f3.close()
        print(dict_gene_type)

    def making_unique_redundant_database_mm39(self, input):
        normal_chromosome_list = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16',
                                  '17', '18', '19', 'X', 'Y', 'MT']
        dict_set = set()
        dict_gene_type = set()
        f2 = open(input.replace('.txt', '_mm39_unique.fasta'), 'w+')
        f3 = open(input.replace('.txt', '_mm39_redundant.fasta'), 'w+')
        # integrated_name = f'>{Gene_ID}_{Transcript_ID}_{Gene_name}#{Chromosome}#{exons_left_list}#{exons_right_list}#{strand}#{CDS_left}#{CDS_right}#{Gene_type}_mRNA'
        with open(input, 'r+') as f1:
            integrated_name, sequence = '', ''
            for line1 in f1:
                if ('>' in line1):
                    if ('miRNA' not in integrated_name) and (integrated_name != ''):
                        tmp_chr = integrated_name.split('_')[2].split('#')[1]
                        if tmp_chr in normal_chromosome_list:  # The newly generated mouse transcript database will not support genes that are not located on non-normal chromosomes.
                            if sequence not in dict_set:
                                f2.write(integrated_name + '\n')
                                f2.write(sequence + '\n')
                                dict_set.add(sequence)
                            elif sequence in dict_set:
                                f3.write(integrated_name + '\n')
                                f3.write(sequence + '\n')
                    line2 = line1.replace('_', '&').strip().split('|')
                    sequence = ''
                    Gene_ID = line2[0][1:]
                    Transcript_ID = line2[1]
                    Gene_type = line2[2].replace('protein&coding', 'mRNA')
                    exons_left_list = line2[3]
                    exons_right_list = line2[4]
                    try:
                        CDS_left = min(int(x) for x in re.findall('\d+', line2[5]))
                        CDS_right = max(int(x) for x in re.findall('\d+', line2[6]))
                    except:
                        CDS_left, CDS_right = '0', '0'
                    strand = line2[7]
                    Chromosome = line2[8]
                    if len(line2) != 10:
                        Gene_name = 'NoName'
                    else:
                        Gene_name = line2[9]
                    dict_gene_type.add(Gene_type)
                    integrated_name = f'>{Gene_ID}_{Transcript_ID}_{Gene_name}#{Chromosome}#{exons_left_list}#{exons_right_list}#{strand}#{CDS_left}#{CDS_right}#{Gene_type}_mRNA'
                else:
                    sequence += line1.strip()
            if 'miRNA' not in integrated_name:
                tmp_chr = integrated_name.split('_')[2].split('#')[1]
                if tmp_chr in normal_chromosome_list:
                    if sequence not in dict_set:
                        f2.write(integrated_name + '\n')
                        f2.write(sequence + '\n')
                        dict_set.add(sequence)
                    elif sequence in dict_set:
                        f3.write(integrated_name + '\n')
                        f3.write(sequence + '\n')
        f2.close()
        f3.close()
        print(dict_gene_type)

    def making_unique_redundant_database_WBcel235(self, input):
        normal_chromosome_list = ['I', 'II', 'III', 'IV', 'V', 'X', 'MtDNA']
        dict_set = set()
        dict_gene_type = set()
        f2 = open(input.replace('.txt', '_WBcel235_unique.fasta'), 'w+')
        f3 = open(input.replace('.txt', '_WBcel235_redundant.fasta'), 'w+')
        # integrated_name = f'>{Gene_ID}_{Transcript_ID}_{Gene_name}#{Chromosome}#{exons_left_list}#{exons_right_list}#{strand}#{CDS_left}#{CDS_right}#{Gene_type}_mRNA'
        with open(input, 'r+') as f1:
            integrated_name, sequence = '', ''
            for line1 in f1:
                if ('>' in line1):
                    if ('miRNA' not in integrated_name) and (integrated_name != ''):
                        tmp_chr = integrated_name.split('_')[2].split('#')[1]
                        if tmp_chr in normal_chromosome_list:  # The newly generated worm transcript database will not support genes that are not located on non-normal chromosomes.
                            if sequence not in dict_set:
                                f2.write(integrated_name + '\n')
                                f2.write(sequence + '\n')
                                dict_set.add(sequence)
                            elif sequence in dict_set:
                                f3.write(integrated_name + '\n')
                                f3.write(sequence + '\n')
                    line2 = line1.replace('_', '&').strip().split('|')
                    sequence = ''
                    Gene_ID = line2[0][1:]
                    Transcript_ID = line2[1]
                    Gene_type = line2[2].replace('protein&coding', 'mRNA')
                    exons_left_list = line2[3]
                    exons_right_list = line2[4]
                    try:
                        CDS_left = min(int(x) for x in re.findall('\d+', line2[5]))
                        CDS_right = max(int(x) for x in re.findall('\d+', line2[6]))
                    except:
                        CDS_left, CDS_right = '0', '0'
                    strand = line2[7]
                    Chromosome = line2[8]
                    if len(line2) != 10:
                        Gene_name = 'NoName'
                    else:
                        Gene_name = line2[9]
                    dict_gene_type.add(Gene_type)
                    integrated_name = f'>{Gene_ID}_{Transcript_ID}_{Gene_name}#{Chromosome}#{exons_left_list}#{exons_right_list}#{strand}#{CDS_left}#{CDS_right}#{Gene_type}_mRNA'
                else:
                    sequence += line1.strip()
            if 'miRNA' not in integrated_name:
                tmp_chr = integrated_name.split('_')[2].split('#')[1]
                if tmp_chr in normal_chromosome_list:
                    if sequence not in dict_set:
                        f2.write(integrated_name + '\n')
                        f2.write(sequence + '\n')
                        dict_set.add(sequence)
                    elif sequence in dict_set:
                        f3.write(integrated_name + '\n')
                        f3.write(sequence + '\n')
        f2.close()
        f3.close()
        print(dict_gene_type)

    def making_transcript_sequence_genomeposition_conservation_database(self, genome_file, transcript_file, chr_name):
        dict_genome = self.loading_genome_database_as_dict(genome_file)
        dict_transcript = self.loading_transcript_databse(transcript_file)
        # f_bedGraph = pd.read_table(f'phyloP100way.bedGrph_chr{chr_name}.bedGraph', index_col=[0, 1], header=None)#human
        # f_bedGraph = pd.read_table(f'mm39.phyloP35way_chr{chr_name}.bedGraph', index_col=[0, 1], header=None)#mouse
        f_bedGraph = pd.read_table(f'ce11.phyloP135way_chr{chr_name}.bedGraph', index_col=[0, 1], header=None)# Worm

        output_file = f"{transcript_file.replace('.fasta', '_')}{chr_name}_conservation_score.txt"
        if os.path.exists(output_file):
            os.remove(output_file)

        f2 = open(output_file, 'a+')
        num1 = 0
        for transcript_name, sequence in dict_transcript.items():
            num1 += 1
            print(num1)
            if transcript_name.count('#') == 7:
                #target_chromosome = transcript_name.split('_')[2].split('#')[1].replace('MT', 'M')#human/mouse
                target_chromosome = transcript_name.split('_')[2].split('#')[1] # Worm
                target_exons_genome_start_sites = sorted(
                    int(x) for x in transcript_name.split('_')[2].split('#')[2].split(';'))
                target_exons_genome_end_sites = sorted(
                    int(x) for x in transcript_name.split('_')[2].split('#')[3].split(';'))
                target_all_sites = zip(target_exons_genome_start_sites,
                                       target_exons_genome_end_sites)  ## the number show the left and right location of each exon in the genome, 1st, this function is iterator
            else:
                target_chromosome = 'not exist'  ## i didn't label miRNA position in the genome
            target_position_in_genome = []
            target_sequence_in_genome = ''

            if target_chromosome == chr_name:
                for exon_sites in target_all_sites:
                    target_position_in_genome += list(range(exon_sites[0], exon_sites[
                        1] + 1))  ## target whole sequence location or position or site in the genome
                    target_sequence_in_genome += dict_genome[target_chromosome][exon_sites[0] - 1:exon_sites[
                        1]]  ## this sequence could same or reverse complementary of transcript (in the negative strand)
                target_position_in_genome_str = ','.join(
                    str(x) for x in target_position_in_genome)  ## switch target site from number to string in list
                transcripts_site_in_genome = list(map(lambda x: ('chr' + target_chromosome, int(x)),
                                                      target_position_in_genome))  ## e.g. ('2L', 42), in case some position not exist in the bedGraph
                f_transcript_conservation_score = pd.DataFrame()
                for chromosome_site in transcripts_site_in_genome:  ## chromosome site: e.g. ('2L', 42)
                    try:  # if conservation score exist in bedGraph
                        f_transcript_conservation_score = pd.concat(
                            [f_transcript_conservation_score, f_bedGraph.loc[[chromosome_site]]])
                    except:  # if conservation score not exist in bedGraph
                        f_transcript_conservation_score = BedGraph().no_exist_conservation_socre_in_bedGraph(
                            target_chromosome, chromosome_site[1], f_transcript_conservation_score)
                each_transcript_conservation_score = ','.join([str(x) for x in f_transcript_conservation_score[2]])
                f2.write('>' + transcript_name + '\n')  ## write transcript name
                f2.write(target_sequence_in_genome + '\n')
                f2.write(target_position_in_genome_str + '\n')  ## write transcript position from genome
                f2.write(each_transcript_conservation_score + '\n')
        f2.close()

    def loading_genome_database_as_dict(self, file):
        dict_genome = {}
        with open(file, 'r+') as f1:
            chromosome, sequence = '', ''
            for line1 in f1:
                if '>' == line1[0]:
                    if chromosome != '':
                        dict_genome[chromosome] = sequence
                    # chromosome = line1.split(' ')[0].strip('>').replace('MT', 'M') # Human/Mouse
                    chromosome = line1.split(' ')[0].strip('>') # Worm
                    sequence = ''
                else:
                    sequence += line1.strip()
            dict_genome[chromosome] = sequence
        return dict_genome

    def transcript_sequence_genomeposition_conservation_database(self, file):
        dict_CS = {}  ## the database including name, sequence, position and conservation score, the sequence from genome, could be reverse complementary from transcript
        with open(file, 'r+') as f1:
            for line1 in f1:
                if line1[0] == '>':
                    name1 = line1.strip()[1:]
                else:
                    dict_CS[name1] = dict_CS.get(name1, [])
                    dict_CS[name1].append(line1.strip())
        return dict_CS

    def loading_transcript_databse(self, file):
        dict_transcipts = {}  ## the databse including transcript, sequence
        with open(file, 'r+') as f1:
            for line1 in f1:
                if line1[0] == '>':
                    name1 = line1.strip()[1:]
                else:
                    dict_transcipts[name1] = (line1.strip())
        return dict_transcipts

    def sequence_name_database(self, file):  ## sequence as  key, name as item
        dict_seq_name = {}
        with open(file, 'r+') as f1:
            for line1 in f1:
                if line1[0] == '>':
                    name1 = line1.strip()[1:]
                else:
                    sequence = line1.strip()
                    dict_seq_name[sequence] = dict_seq_name.get(sequence, [])
                    dict_seq_name[sequence].append(name1)
        return dict_seq_name

    def unique_redundant_geneID_dict(self, unique_database, redundant_database):
        unique_redundant_name_dict1 = {}
        unique_transcripts = self.sequence_name_database(unique_database)
        redundant_transcripts = self.sequence_name_database(redundant_database)
        for seq, redundant_name_list in redundant_transcripts.items():
            unique_geneID = unique_transcripts[seq][0].split('_')[0]  # remove gene name, only keep gene_ID
            for redundant_name in redundant_name_list:
                redundant_geneID = redundant_name.split('_')[0]
                redundant_transcriptID = redundant_name.split('_')[1]
                redundant_strand = redundant_name.split('_')[2].split('#')[4]
                redundant_exon_min = min(int(x) for x in redundant_name.split('_')[2].split('#')[2].split(';'))
                redundant_exon_max = max(int(x) for x in redundant_name.split('_')[2].split('#')[3].split(';'))
                redundant_CDS_min = min(int(x) for x in redundant_name.split('_')[2].split('#')[5].split(';'))
                redundant_CDS_max = min(int(x) for x in redundant_name.split('_')[2].split('#')[6].split(';'))
                redundant_exon_range_str = f"{redundant_exon_min}-{redundant_exon_max}"
                redundant_CDS_range_str = f"{redundant_CDS_min}-{redundant_CDS_max}"
                if redundant_geneID != unique_geneID:
                    unique_redundant_name_dict1[unique_geneID] = unique_redundant_name_dict1.get(unique_geneID, {})
                    unique_redundant_name_dict1[unique_geneID][redundant_transcriptID] = unique_redundant_name_dict1[
                        unique_geneID].get(redundant_transcriptID, {'strand': '', 'exon_range': '', 'CDS_range': ''})
                    unique_redundant_name_dict1[unique_geneID][redundant_transcriptID]['strand'] = redundant_strand
                    unique_redundant_name_dict1[unique_geneID][redundant_transcriptID][
                        'exon_range'] = redundant_exon_range_str
                    unique_redundant_name_dict1[unique_geneID][redundant_transcriptID][
                        'CDS_range'] = redundant_CDS_range_str
        return unique_redundant_name_dict1


class Viennad_to_table():
    def __init__(self, transcript_ConservationScore_database, transcript_only_database):
        self.dict_CS = Database().transcript_sequence_genomeposition_conservation_database(
            file=transcript_ConservationScore_database)
        self.dict_transcript = Database().loading_transcript_databse(file=transcript_only_database)

    def input_viennad(self, file):
        viennad_file = file + '.viennad'
        print(viennad_file)
        with open(viennad_file, 'r+') as f1:
            list_hyb = []
            dict_hyb = {}
            no_exist_transcript_in_dict_cs = set()
            for line1 in f1:
                line2 = line1.strip().split('_')
                if len(line2) == 15:  # microRNA on the left of hybrid, type=mim pref=mim
                    if (list_hyb != []) and (len(list_hyb[4].split('\t')) == 2) and ('microRNA' in list_hyb[0]) and (
                            'mRNA' in list_hyb[0]):
                        hyb_each_num = int(list_hyb[0].split('_')[1])
                        miRNA_name = list_hyb[0].split('_')[4]
                        miRNA_seq = list_hyb[2].split('\t')[0].strip('-')  ##### miRNA sequence
                        miRNA_length = len(miRNA_seq)
                        miRNA_pattern = list_hyb[4][:miRNA_length]  ##### miRNA pattern
                        miRNA_pattern_short = miRNA_pattern.strip('.')  ## remove '.' in the  outside of '('
                        miRNA_pattern_short_pattern = ''.join('\\' + x for x in miRNA_pattern_short)
                        miRNA_unpaired_5prime_length = re.search(miRNA_pattern_short_pattern, miRNA_pattern).span()[0]
                        miRNA_unpaired_3prime_length = miRNA_length - \
                                                       re.search(miRNA_pattern_short_pattern, miRNA_pattern).span()[1]
                        gene_id = list_hyb[0].split('_')[9]
                        transcript_id = list_hyb[0].split('_')[10]
                        target_name = list_hyb[0].split('_')[11]
                        target_pattern_element_only = list_hyb[4].split('\t')[0][miRNA_length:].strip(
                            '.')  ##### target pattern
                        target_pattern_element_only_format = ''.join('\\' + x for x in target_pattern_element_only)
                        target_seq_in_Viennad = list_hyb[3].split('\t')[0].strip(
                            '-')  ## the sequence from viennad, i will get longer sequence from transcript
                        target_pattern_in_Viennad = list_hyb[4].split('\t')[0][miRNA_length:]
                        transcript_name = '_'.join([gene_id, transcript_id, target_name, 'mRNA'])
                        target_seq_completed = '-' * 20 + self.dict_transcript[
                            transcript_name] + '-' * 20  ## add 20 '-' in the transcript sequence
                        target_seq_in_Viennad_start = (re.search(target_seq_in_Viennad, target_seq_completed).start())
                        target_seq_in_Viennad_end = (re.search(target_seq_in_Viennad, target_seq_completed).end())
                        target_seq_in_Viennad_extend20nt = target_seq_completed[
                                                           target_seq_in_Viennad_start - 20: target_seq_in_Viennad_end + 20]
                        target_pattern_in_Viennad_extended = '.' * 20 + target_pattern_in_Viennad + '.' * 20
                        target_pattern_element_only_span = (
                            re.search(target_pattern_element_only_format, target_pattern_in_Viennad_extended)).span()
                        target_start_num = target_pattern_element_only_span[0] - miRNA_unpaired_3prime_length
                        target_end_num = target_pattern_element_only_span[1] + miRNA_unpaired_5prime_length
                        target_seq_element = target_seq_in_Viennad_extend20nt[
                                             target_start_num:target_end_num]  #### target sequence with unpaired nucleotides
                        target_pattern_element = target_pattern_in_Viennad_extended[target_start_num:target_end_num]
                        strand = list_hyb[0].split('_')[11].split('#')[4]
                        energy = list_hyb[4].split('\t')[1][1:-1]

                        hyb_name = '_'.join(
                            [miRNA_name, miRNA_seq, miRNA_pattern, gene_id, transcript_id, target_name,
                             target_seq_element, target_pattern_element, energy])
                        dict_hyb[hyb_name] = dict_hyb.get(hyb_name,
                                                          {'abundance': 0, 'conservation_score': set(),
                                                           'genome_position': set()})
                        if strand == '-1':
                            target_seq_in_Viennad = str(Seq(target_seq_in_Viennad).reverse_complement())
                            target_pattern_in_Viennad = target_pattern_in_Viennad[
                                                        ::-1]  ## for identify the paired nucleotides position in the genome
                        if transcript_name in self.dict_CS:
                            target_position_in_genome_span = re.search(target_seq_in_Viennad,
                                                                       self.dict_CS[transcript_name][0]).span()
                            target_element_position = self.dict_CS[transcript_name][1].split(',')[
                                                      target_position_in_genome_span[0]:target_position_in_genome_span[
                                                          1]]
                            target_element_conservation_score = [float(x) for x in
                                                                 self.dict_CS[transcript_name][2].split(',')[
                                                                 target_position_in_genome_span[0]:
                                                                 target_position_in_genome_span[1]]]
                            target_element_data = {'pattern': list(target_pattern_in_Viennad),
                                                   'genome_seq': list(target_seq_in_Viennad),
                                                   'position': target_element_position,
                                                   'CS': target_element_conservation_score}
                            target_element_dataframe = pd.DataFrame(data=target_element_data)
                            target_element_dataframe = target_element_dataframe[
                                target_element_dataframe['pattern'] == ')']  ## remove unpaired position
                            target_element_only_conservation_score = {np.around(target_element_dataframe['CS'].mean(),
                                                                                decimals=2)}  ## conservation score without unpaired position
                            dict_hyb[hyb_name]['conservation_score'].update(target_element_only_conservation_score)
                            dict_hyb[hyb_name]['genome_position'].update(set(target_element_dataframe['position']))
                        if (transcript_name not in self.dict_CS) and (
                                transcript_name not in no_exist_transcript_in_dict_cs):
                            no_exist_transcript_in_dict_cs.add(transcript_name)
                            print(f"{transcript_name} is not in dict_CS")  ## no conservation socre out
                        dict_hyb[hyb_name]['abundance'] += hyb_each_num
                    list_hyb = []
                list_hyb.append(line1.strip())
        f2 = open(viennad_file.replace('viennad', 'txt'), 'w+')
        f2.write(
            f'miRNA_name\tmiRNA_sequence\tmiRNA_pattern\tGene_ID\ttranscript_ID\tGene_information\telement_sequence\telement_pattern\tdG\tabundance\tgenome_position\tConservation_score\n')
        for name in dict_hyb:
            if dict_hyb[name][
                'conservation_score'] == set():  #
                f2.write('\t'.join(name.split('_')) + '\t' + str(dict_hyb[name]['abundance']) + '\n')
            elif dict_hyb[name]['conservation_score'] != set():
                genome_position_str = str(sorted([int(x) for x in dict_hyb[name]['genome_position']]))
                if len(dict_hyb[name]['conservation_score']) == 1:
                    f2.write('\t'.join(name.split('_')) + '\t' + str(
                        dict_hyb[name]['abundance']) + '\t' + genome_position_str + '\t' + str(
                        max(dict_hyb[name]['conservation_score'])) + '\n')
                else:  # there are multiple conservation scores, 0.1% type of  hybs have multiple CS, because their repetitive sequence in the genome. Ming and Nick suggest me to keep all those, but delete conservation score. Lable as multiple repetitive element.
                    f2.write(
                        '\t'.join(name.split('_')) + '\t' + str(
                            dict_hyb[name]['abundance']) + '\t' + 'multiple_elements' + '\t' + str(
                            max(dict_hyb[name][
                                    'conservation_score'])) + '\n')  # there are multiple element site in one transcript, i will put the largest conservation score, but i don't show the position
        f2.close()


class Combined_table():
    def __init__(self, replicates=2, outfile='out.csv'):
        self.replicates = replicates
        self.outfile = outfile

    def clash_table(self, file):
        f = pd.read_table(file,
                          index_col=['miRNA_name', 'miRNA_sequence', 'miRNA_pattern', 'Gene_ID', 'element_sequence',
                                     'element_pattern', 'dG', 'genome_position', 'Conservation_score', 'Gene_type',
                                     'Gene_name', 'Chromosome', 'element_region'])
        f['abundance'] = (f['abundance'] / f['abundance'].sum()) * 1000000
        return f

    def input_multiple_raw_clash(self, viennad_files):
        files = [re.findall(".+txt", x)[0] for x in str(viennad_files).split(',') if re.findall(".+txt", x) != []]
        print(files)
        f_all = pd.DataFrame()
        index_col1 = ['miRNA_name', 'miRNA_sequence', 'miRNA_pattern', 'Gene_ID', 'element_sequence',
                      'element_pattern', 'dG', 'genome_position', 'Conservation_score', 'Gene_type',
                      'Gene_name', 'Chromosome', 'element_region']
        for file_name1 in files:
            print(file_name1)
            f = pd.read_table(file_name1, index_col=index_col1)
            f.rename(columns={"abundance": file_name1}, inplace=True)
            f_all = pd.concat([f_all, f], axis=1)
        f_all.to_csv(self.outfile, index_label=index_col1, sep='\t')

    def input_multiple_clash(self, files):  ##the first number is replicates, the others are files name
        files = [re.findall("\w.+\w+.txt", x)[0] for x in str(files).split(',') if re.findall("\w.+\w+.txt", x) != []]
        print(files)
        f = pd.concat([self.clash_table(x) for x in files], axis=1)
        old_columns = set(f.columns)
        print(old_columns)
        f = f.query('dG<-11')

        try:
            count_col_abundance = f['abundance'].shape[1]  # how many columns
            f = f[f['abundance'].isna().sum(axis=1) <= count_col_abundance - int(
                self.replicates)]  ## remove low replicates
            f['mean_abundance'] = f['abundance'].mean(axis=1)
            f['mean_miRNA_ex_ratio'] = f['miRNA_ex_ratio'].mean(axis=1)
            f['mean_hyb_ex_ratio'] = f['hyb_ex_ratio'].mean(axis=1)
        except:  # only one abundance
            f.rename(columns={"abundance": "mean_abundance", "miRNA_ex_ratio": "mean_miRNA_ex_ratio",
                              "hyb_ex_ratio": "mean_hyb_ex_ratio"}, inplace=True)
            for columns_name1 in ['abundance', 'hyb_ex_ratio', 'miRNA_ex_ratio']:
                old_columns.remove(columns_name1)
        f.drop(columns=list(old_columns), inplace=True)
        f = f.round(decimals=3)
        f.to_csv(self.outfile, sep='\t')
        return f

    def add_redundant_hyb_in_table(self, unique_file, redundant_file, CLASH_table):
        dict_unique_redundant_name = Database().unique_redundant_geneID_dict(unique_database=unique_file,
                                                                             redundant_database=redundant_file)
        ## this dictionary contain the redundant geneID_gene name related to unique geneID
        outfile = CLASH_table.replace('.txt', '_redundant.txt')
        if os.path.exists(outfile):
            os.remove(outfile)
        f2 = open(outfile, 'a+')
        with open(CLASH_table, 'r+') as f1:
            lines = f1.readlines()
            line1st = lines[:1][0]
            f2.write(line1st)
            lines_others = lines[1:]
            for line1 in lines_others:
                f2.write(line1)
                line2 = line1.split('\t')
                geneID = line2[3]
                geneName = line2[10]
                if geneID in dict_unique_redundant_name:
                    for redundant_geneID_geneName in dict_unique_redundant_name[geneID]:
                        redundant_geneID = redundant_geneID_geneName.split('_')[0]
                        redundant_geneName = redundant_geneID_geneName.split('_')[1]
                        line_out = line1.replace(geneID, redundant_geneID).replace(geneName, redundant_geneName)
                        f2.write(line_out)
        f2.close()

        line1st_list = line1st.rstrip('\n').split('\t')
        line1st_list.remove('mean_abundance')
        f = pd.read_table(outfile, index_col=line1st_list)
        f.sort_values(by='mean_abundance', ascending=False, inplace=True)
        print(f[f.index.duplicated(keep="first")])
        f = f[~f.index.duplicated(keep="first")]
        os.remove(outfile)
        f.to_csv(outfile, sep='\t')

    def TDMD_candidate(self, input_file):
        output1 = input_file.replace('.txt', '_TDMD.txt')
        if os.path.exists(output1):
            os.remove(output1)

        f1 = pd.read_table(input_file, header=[0])
        f1 = f1.query('dG<-16')  ## TDMD base pattern, dG <= -16
        first_line1 = '\t'.join(list(f1.columns))

        f2 = open(output1, 'a+')
        f2.write(first_line1 + '\n')
        print(
            'TDMD analysis criteria\n1. miRNA seed region 2-8 pair with target\n2. abs|bulge| <= 6nt\n3. P2: 7 continuous pairing in last 8nt or 9 continuous pairing in any position.\n4. dG <= 16kcal/mol\n')
        for index1 in f1.index:
            miRNA_pattern = f1.loc[index1]['miRNA_pattern']
            target_pattern = f1.loc[index1]['element_pattern']
            if (miRNA_pattern[:8] == '((((((((' and target_pattern[-8:] == '))))))))') or (
                    miRNA_pattern[:8] == '.(((((((' and target_pattern[-8:] == '))))))).'):  # allow P1 pairing
                if (('(((((((' in miRNA_pattern[-8:]) and (')))))))' in target_pattern[:8])) or (
                        '(((((((((' in miRNA_pattern[9:] and ')))))))))' in target_pattern[:-9]):  # allow p2 pairing
                    if abs(len(target_pattern) - len(miRNA_pattern)) <= 6:  # allow bulge less than 6nt
                        each_line1 = '\t'.join([str(x) for x in f1.loc[index1]])
                        f2.write(f"{each_line1}\n")
        f2.close()


class Compressed_table():  # current, this code only surpport 1 abundance, because i want to sum all abundance of similiar element

    def miRNA_element_range(self, list1):
        list_cluster = []
        lowest_num, maximum_num = 0, 0
        for number1 in sorted(list(list1)):
            if (lowest_num == 0) and (maximum_num == 0):
                lowest_num, maximum_num = number1, number1
            else:
                if int(number1) - 1 == maximum_num:
                    maximum_num = int(number1)
                elif int(number1) - 1 != maximum_num:
                    cluster1 = f'{lowest_num}-{maximum_num}'
                    list_cluster.append(cluster1)
                    lowest_num, maximum_num = number1, number1
        cluster1 = f'{lowest_num}-{maximum_num}'
        list_cluster.append(cluster1)
        return list_cluster

    def compressed_each_element(self, dict1):
        highest_abundance_each_hyb = 0
        highest_abundance_hyb_name = ''
        total_abundance_each_hyb = 0
        for each_hyb, abundance in dict1.items():
            if highest_abundance_each_hyb < abundance:
                highest_abundance_each_hyb = abundance
                highest_abundance_hyb_name = each_hyb
            total_abundance_each_hyb += abundance
        return highest_abundance_hyb_name, total_abundance_each_hyb

    def compressed_same_index_table_to_dict(self, table, compressed_index):
        dict_compressed_table = {}
        for index1 in compressed_index:
            dict_compressed_table[index1] = int(table.loc[index1].sum())
        return dict_compressed_table

    def files(self, input):
        input_table = input + '.txt'
        output = input_table.replace('.txt', '_compressed.txt')
        if os.path.exists(output):
            os.remove(output)

        f1 = pd.read_table(input_table,
                           index_col=["miRNA_name", "miRNA_sequence", "miRNA_pattern", "Gene_ID",
                                      "element_sequence", "element_pattern", "dG",
                                      "genome_position", "Conservation_score"])

        f1 = (f1[['Gene_information', 'abundance']])
        f1['Gene_name'] = f1['Gene_information'].str.split('#', expand=True)[0]
        f1['Gene_type'] = f1['Gene_information'].str.split('#', expand=True)[7]
        f1['Chromosome'] = f1['Gene_information'].str.split('#', expand=True)[1]
        f1.drop(columns=['Gene_information'], inplace=True)
        f1.set_index('Gene_type', append=True, inplace=True)
        f1.set_index('Gene_name', append=True, inplace=True)
        f1.set_index('Chromosome', append=True, inplace=True)

        f_genome_position_NaN_index_compressed = set()  # make a compressed list
        f_genome_position_NaN_index_compressed_str = set()  # the index is set as string, if not set as string, i found there can be same index1 in the set()
        f_normal_index_compressed = set()  # make a compressed list
        f_normal_index_compressed_str = set()  # the index is set as string, if not set as string, i found there can be same index1 in the set()
        for index1 in f1.index:
            if (np.nan in index1) or ('multiple_elements' in index1):
                if str(index1) not in f_genome_position_NaN_index_compressed_str:
                    f_genome_position_NaN_index_compressed_str.update({
                        str(index1)})  ## this set only for checking the index, if index name is same, do not put real index into f_genome_position_NaN_index_compressed
                    f_genome_position_NaN_index_compressed.update({index1})
            if (np.nan not in index1) and ('multiple_elements' not in index1):
                if str(index1) not in f_normal_index_compressed_str:
                    f_normal_index_compressed_str.update({str(index1)})
                    f_normal_index_compressed.update({index1})

        f_genome_position_NaN_index = list(f_genome_position_NaN_index_compressed)
        f_genome_position_NaN = (f1.loc[f_genome_position_NaN_index])
        dict_f_genome_position_NaN = self.compressed_same_index_table_to_dict(table=f_genome_position_NaN,
                                                                              compressed_index=f_genome_position_NaN_index)  ## compressecd multiple element and no genome position table

        f_output = open(output, 'a+')  ##output Nan conservation table
        f_output.write(
            f'miRNA_name\tmiRNA_sequence\tmiRNA_pattern\tGene_ID\telement_sequence\telement_pattern\tdG\tgenome_position\tConservation_score\tGene_type\tGene_name\tChromosome\tabundance\n')
        for index1, abundance1 in dict_f_genome_position_NaN.items():  # miRNA_site_ranges_dict: {'137-148':{each_hyb: abundance}}
            each_hyb_list = [str(x) for x in index1]
            each_hyb = '\t'.join((each_hyb_list))
            f_output.write(f'{each_hyb}\t{abundance1}\n')
        f_output.close()

        f_normal_index_compressed = list(f_normal_index_compressed)
        f_normal = f1.loc[f_normal_index_compressed]  ##for most of hybrids, including genome position
        dict_f_normal = self.compressed_same_index_table_to_dict(table=f_normal,
                                                                 compressed_index=f_normal_index_compressed)

        dict_miRNA_target_chromosome_position = {}  ## combine all miRNA_target postion in this dictionary, for remove the part element hybrids from clash table, e.g. {'miRNA_target_chromosome': {3,4,5,6,7,8}'}
        for index1, abundance in dict_f_normal.items():
            # for loop 1st, collect all miRNA binding position in the genome
            miRNA_name1 = index1[0]
            gene_id1 = index1[3]
            chromosome1 = index1[11]
            miRNA_geneID_chr = '_'.join([miRNA_name1, gene_id1, chromosome1])
            miRNA_site_genome_set = set(int(x) for x in re.findall(r'\d+', index1[7]))
            miRNA_site_genome_min = int(min(miRNA_site_genome_set))
            miRNA_site_genome_max = int(max(miRNA_site_genome_set))
            dict_miRNA_target_chromosome_position[miRNA_geneID_chr] = dict_miRNA_target_chromosome_position.get(
                miRNA_geneID_chr, set())
            dict_miRNA_target_chromosome_position[miRNA_geneID_chr].update(
                set(range(miRNA_site_genome_min,
                          miRNA_site_genome_max + 1)))  ## this miRNA element region including bulge position

        dict_hyb_abundance = {}  ## combine all hybrids in this dictionary, determine which is the highest abundance one element position in the genome, e.g. {'miRNA_target_chromosome': {'3452-3476': {hyb_name: abundance}}'}
        for index1, abundance in dict_f_normal.items():
            # for loop 2nd, collect each_hyb into element range e.g. {'3452-3476': {hyb_name: abundance}
            chromosome1 = index1[11]
            each_hyb_list = [str(x) for x in index1]
            each_hyb = '\t'.join(each_hyb_list)
            each_hyb_abundance = abundance
            element_position_in_genome = set([int(x) for x in re.findall(r'\d+', index1[7])])
            miRNA_name1 = index1[0]
            gene_id1 = index1[3]
            miRNA_geneID_chr = '_'.join([miRNA_name1, gene_id1, chromosome1])
            list_hyb_position = self.miRNA_element_range(dict_miRNA_target_chromosome_position[
                                                             miRNA_geneID_chr])  ## def, calculate the lowest and highest position in the genome, ['134-156', '167-178']
            for element_range in list_hyb_position:  ## e.g. element_range : '123-145'
                # inner of 2nd for loop
                miRNA_site_range_set = set(range(int(element_range.split('-')[0]), int(
                    element_range.split('-')[1]) + 1))  ##In the each_hyb, there could be multiple range e.g.(167-178),
                if miRNA_site_range_set.intersection(
                        element_position_in_genome) != set():  # Does the element binding site in each row & one of multiple element_range from genome
                    dict_hyb_abundance[miRNA_geneID_chr] = dict_hyb_abundance.get(miRNA_geneID_chr, {})
                    dict_hyb_abundance[miRNA_geneID_chr][element_range] = dict_hyb_abundance[miRNA_geneID_chr].get(
                        element_range, {})
                    dict_hyb_abundance[miRNA_geneID_chr][element_range][each_hyb] = \
                        dict_hyb_abundance[miRNA_geneID_chr][element_range].get(each_hyb, 0)
                    dict_hyb_abundance[miRNA_geneID_chr][element_range][each_hyb] += each_hyb_abundance

        f_output = open(output, 'a+')
        for miRNA_geneID_chr, element_ranges_dict in dict_hyb_abundance.items():
            for each_element_range, each_hyb_dict in element_ranges_dict.items():  ## element_ranges_dict: {'137-148':{each_hyb1: abundance,each1_hyb2: abundance2},' 152-167':{each_hyb: abundance} }
                highest_abundance_hyb_name = self.compressed_each_element(each_hyb_dict)[0]
                total_abundance_each_hyb = self.compressed_each_element(each_hyb_dict)[1]
                f_output.write(f'{highest_abundance_hyb_name}\t{total_abundance_each_hyb}\n')
        f_output.close()


class Gff:
    def __init__(self, hyb_ratio=0.01):
        self.hyb_ratio = hyb_ratio

    def extract_top_abundance_table_dict(self, file):  ## extract top abudance hyb into dict
        f1 = pd.read_table(file, index_col=["miRNA_name", "miRNA_sequence", "miRNA_pattern", "Gene_ID",
                                            "Gene_name", "Gene_type", "Chromosome", "element_sequence",
                                            "element_pattern", "dG", "genome_position", "Conservation_score"])
        dict_geneID_abundance = {}
        for index1 in f1.index:
            # 1st loop, sum each target abundance from viennad file
            geneID = index1[3]
            each_geneID_abundance = f1.loc[index1]
            dict_geneID_abundance[geneID] = dict_geneID_abundance.get(geneID, 0)
            dict_geneID_abundance[geneID] += each_geneID_abundance

        dict_top_abundance_hyb = {}
        for index1 in f1.index:
            each_hyb = '\t'.join([str(x) for x in index1])
            geneID = index1[3]
            each_geneID_abundance = float(f1.loc[index1])
            each_geneID_abundance_ratio = np.around(float(f1.loc[index1] / dict_geneID_abundance[geneID]), decimals=3)
            if each_geneID_abundance_ratio >= self.hyb_ratio:
                dict_top_abundance_hyb[each_hyb] = each_geneID_abundance
        return dict_top_abundance_hyb

    def element_gff_each_row(self,
                             list1):  ## the shortest of intron >= 30bp  <<Piovesan,  et al. DNA Research 22.6 (2015): 495-503.>>
        list_miRNA_sites = []
        lowest_num, maximum_num = 0, 0
        for index, number1 in enumerate(list1):
            if index == 0:
                lowest_num = int(number1)
                maximum_num = int(number1)
            else:
                if int(number1) - 30 <= maximum_num:
                    maximum_num = int(number1)
                else:  ## there is a gap
                    list_miRNA_sites.append([lowest_num, maximum_num])
                    lowest_num = int(number1)
                    maximum_num = int(number1)
        list_miRNA_sites.append([lowest_num, maximum_num])
        return list_miRNA_sites

    def make_gff_file(self, infile):
        outfile = infile.replace('.txt', '.gff')
        if os.path.exists(outfile):
            os.remove(outfile)
        f2 = open(outfile, 'a+')
        dict_top_hyb = self.extract_top_abundance_table_dict(infile)
        for each_hyb, abundance in dict_top_hyb.items():
            element_position = each_hyb.split('\t')[10]
            gene_ID = each_hyb.split('\t')[3]
            gene_seq = each_hyb.split('\t')[7]
            chromosome = each_hyb.split('\t')[6]
            miRNA_name = each_hyb.split('\t')[0]
            miRNA_seq = each_hyb.split('\t')[1]
            if (str(element_position) != str(np.nan)) and (str(element_position) != 'multiple_elements'):
                element_position_list = re.findall('\d+', element_position)
                for miRNA_exon in self.element_gff_each_row(element_position_list):
                    f2.write(
                        f'chr{chromosome}\tFlyBase\texon\t{miRNA_exon[0]}\t{miRNA_exon[1]}\t.\t.\t.\tParent={miRNA_name}_{miRNA_seq}_{gene_ID}_{gene_seq}\n')
        f2.close()


class Gene_region():
    def __init__(self, ensemble_database):
        self.ensemble_database = ensemble_database

    def geneID_transcriptID_CDS_exon(self, file):  # including CDS and exon position
        dict_geneID_transcriptID_CDS = {}
        with open(file, 'r+') as f1:
            for line1 in f1:
                if 'protein_coding' in line1:
                    line2 = line1.strip().split('|')
                    Gene_ID = line2[0][1:]
                    transcript_ID = line2[1]
                    strand = line2[7]
                    exon_left_list = line2[3]
                    exon_right_list = line2[4]
                    CDS_left_list = line2[5]
                    CDS_right_list = line2[6]
                    if CDS_left_list != '':
                        target_exons_genome_start_sites = sorted(
                            int(x) for x in exon_left_list.split(';'))
                        target_exons_genome_end_sites = sorted(
                            int(x) for x in exon_right_list.split(';'))
                        target_CDS_genome_start_sites = sorted(
                            int(x) for x in CDS_left_list.split(';'))
                        target_CDS_genome_end_sites = sorted(
                            int(x) for x in CDS_right_list.split(';'))

                        dict_geneID_transcriptID_CDS[Gene_ID] = dict_geneID_transcriptID_CDS.get(Gene_ID, {})
                        dict_geneID_transcriptID_CDS[Gene_ID][transcript_ID + '_' + strand] = \
                            dict_geneID_transcriptID_CDS[Gene_ID].get(transcript_ID + '_' + strand,
                                                                      {'CDSs': set(), 'Exons': set()})

                        for exon_each in zip(target_exons_genome_start_sites, target_exons_genome_end_sites):
                            dict_geneID_transcriptID_CDS[Gene_ID][transcript_ID + '_' + strand]['Exons'].update(
                                set(range(exon_each[0], exon_each[1] + 1)))
                        for CDS_each in zip(target_CDS_genome_start_sites, target_CDS_genome_end_sites):
                            dict_geneID_transcriptID_CDS[Gene_ID][transcript_ID + '_' + strand]['CDSs'].update(
                                set(range(CDS_each[0], CDS_each[1] + 1)))
        return dict_geneID_transcriptID_CDS

    def table(self, input):
        input_table = input + '_compressed_AUex.txt'
        output = input_table.replace('.txt', '_region.txt')
        if os.path.exists(output):
            os.remove(output)
        dict_cds = self.geneID_transcriptID_CDS_exon(file=self.ensemble_database)
        with open(input_table, 'r+') as f1:
            line1 = f1.readlines()
            line1st = line1[0]
            line_other = line1[1:]

        f2 = open(output, 'a+')
        f2.write(line1st.rstrip('\n') + '\t' + 'element_region' + '\n')
        for each_row in line_other:
            index1 = each_row.rstrip('\n').split('\t')
            element_position = set(int(x) for x in re.findall('\d+', index1[7]))
            list_region = []
            geneID = index1[3]
            if (element_position != set()) and ('mRNA' == index1[9]):
                if geneID in dict_cds:
                    transcriptName_list = dict_cds[geneID]
                    for transcriptName_strand in transcriptName_list:
                        strand = transcriptName_strand.split('_')[1]
                        CDS_region_set = dict_cds[geneID][transcriptName_strand]['CDSs']
                        exon_region_set = dict_cds[geneID][transcriptName_strand]['Exons']
                        if exon_region_set.intersection(
                                element_position) != set():  ## the element in the exon region, cds or utr,
                            if CDS_region_set.intersection(
                                    element_position) != set():  # the transcript exon in genome should overlap with target element
                                list_region.append('CDS')
                            if strand == '-1':
                                if max(element_position) < min(CDS_region_set):
                                    list_region.append('3UTR')
                                if min(element_position) > max(CDS_region_set):
                                    list_region.append('5UTR')
                            if strand == '1':
                                if min(element_position) > max(CDS_region_set):
                                    list_region.append('3UTR')
                                if max(element_position) < min(CDS_region_set):
                                    list_region.append('5UTR')
                    if list_region == []:
                        list_region.append(
                            'intron')  ## this is what i guess, some element in the exon, but not in all transcripts CDS
                    list_region_str = ';'.join(sorted(set(list_region)))
                    f2.write(each_row.rstrip('\n') + '\t' + str(list_region_str) + '\n')
                elif geneID not in dict_cds:  ## not mRNA
                    f2.write(f'{each_row}')
            else:
                f2.write(f'{each_row}')
        f2.close()


class MiRNA_CLASH_extension():
    def __init__(self, infile=''):
        self.infile = infile + '_compressed.txt'

    def viennad_list(self, list1, dict1):
        miRNA_name1 = list1[0].split('_')[4]
        gene_ID1 = list1[0].split('_')[9]
        seriers_number = list1[0].split('_')[0] + '_' + list1[0].split('_')[1]
        miRNA_end = int(list1[2].split('\t')[3])
        miRNA_basepattern = list1[4][:miRNA_end]

        target_pattern_element_only = list1[4].split('\t')[0][miRNA_end:].strip('.')  ##### target pattern
        target_pattern_element_only_format = ''.join('\\' + x for x in target_pattern_element_only)
        target_pattern_extension = list1[4].split('\t')[0][miRNA_end:]
        target_pattern_element_only_span = (
            re.search(target_pattern_element_only_format, target_pattern_extension)).span()
        target_seq_extension = list1[3].split('\t')[0].strip('-')
        target_seq_element_only = target_seq_extension[
                                  target_pattern_element_only_span[0]:target_pattern_element_only_span[
                                      1]]  #### target  sequence
        combine_name1 = f'{miRNA_name1}${gene_ID1}${miRNA_basepattern}${target_pattern_element_only}${target_seq_element_only}'
        dict1[seriers_number] = combine_name1
        return dict1

    def each_pattern_in_viennad(self,
                                viennad_file):  ## Extract the series number and the corresponding base pattern\n",
        with open(viennad_file, 'r+') as f1:
            dict1 = {}  ## find each base pattern in viennad file
            list1 = []
            for line1 in f1:
                line2 = line1.strip().split('_')
                if len(line2) >= 15:
                    if list1 != []:
                        try:
                            dict1 = self.viennad_list(list1, dict1)
                        except:
                            print(list1, 'wrong hyb in viennad')
                    list1 = []
                list1.append(line1.strip())
            try:
                dict1 = self.viennad_list(list1, dict1)
            except:
                print(list1, 'wrong hyb in viennad')
        return dict1

    def miRNA_AU_ex_number(self, hyb_viennad_names):
        hyb_viennad_names = [re.findall("\w.+\w+", x)[0] for x in str(hyb_viennad_names).split(',') if
                             re.findall("\w.+\w+", x) != []]
        dict_miRNAex_in_hyb = {}  ## get complete and short miRNA of each base pattern, equal 'combine_name1'
        dict_eachmiRNA_ex = {}  ## miRNA total number and total extension
        for file1 in hyb_viennad_names:
            viennad_file = file1 + '.viennad'
            hyb_file = file1 + '.hyb'
            dict_series_base = self.each_pattern_in_viennad(viennad_file)
            with open(hyb_file) as f1:
                reader = csv.reader(f1, delimiter='\t')
                for line1 in reader:
                    if line1[
                        0] in dict_series_base:  ## series number from hyb in viennad, very little hyb in veinnad is wrong, they are not in dict_series_base
                        base_pattern_name1 = dict_series_base[
                            line1[0]]  ## base pattern correspond to Series repeat number
                        dict_miRNAex_in_hyb[base_pattern_name1] = dict_miRNAex_in_hyb.get(base_pattern_name1,
                                                                                          {'total': 0, 'extension': 0})
                        dict_miRNAex_in_hyb[base_pattern_name1]['total'] += int(line1[0].split('_')[-1])
                        each_miRNA = base_pattern_name1.split('$')[0]
                        dict_eachmiRNA_ex[each_miRNA] = dict_eachmiRNA_ex.get(each_miRNA, {'total': 0,
                                                                                           'extension': 0})  ## each miRNA extension
                        dict_eachmiRNA_ex[each_miRNA]['total'] += int(line1[0].split('_')[-1])

                        if ('microRNA' in line1[3]):  ## miRNA on the left
                            pos_microRNA_end = int(line1[5])
                            pos_mRNA_start = int(line1[10])
                            seq_insert = line1[1][
                                         pos_microRNA_end:pos_mRNA_start - 1]  ### insert between microRNA and target
                            if (2 <= len(seq_insert) <= 8) and (('A' in seq_insert) or ('T' in seq_insert)):
                                dict_miRNAex_in_hyb[base_pattern_name1]['extension'] += int(line1[0].split('_')[-1])
                                dict_eachmiRNA_ex[each_miRNA]['extension'] += int(line1[0].split('_')[-1])
                        elif ('microRNA' in line1[9]):  ## microRNA on the right
                            seq_insert = line1[1][int(line1[11]):]  ### insert at the end of target
                            if (2 <= len(seq_insert) <= 8) and (('A' in seq_insert) or ('T' in seq_insert)):
                                dict_miRNAex_in_hyb[base_pattern_name1]['extension'] += int(line1[0].split('_')[-1])
                                dict_eachmiRNA_ex[each_miRNA]['extension'] += int(line1[0].split('_')[-1])
                    else:
                        print(f"{line1[0]} is not in viennad dict")
        return dict_eachmiRNA_ex, dict_miRNAex_in_hyb

    def miRNA_AU_ex_ratio(self, hyb_viennad_names):
        outfile = self.infile.replace('.txt', '_AUex.txt')
        print(self.infile)
        print(outfile)
        dict_eachmiRNA_ex = self.miRNA_AU_ex_number(hyb_viennad_names)[0]
        dict_each_miRNA_extension_ratio = {}
        for miRNA_name in dict_eachmiRNA_ex:
            miRNA_extension_ratio = np.around(
                dict_eachmiRNA_ex[miRNA_name]['extension'] / dict_eachmiRNA_ex[miRNA_name]['total'], decimals=3)
            dict_each_miRNA_extension_ratio[miRNA_name] = miRNA_extension_ratio

        dict_miRNAex_in_hyb_ratio = {}
        dict_miRNAex_in_hyb = self.miRNA_AU_ex_number(hyb_viennad_names)[1]
        for each_hyb in dict_miRNAex_in_hyb:
            miRNAex_in_hyb_ratio = np.around(
                dict_miRNAex_in_hyb[each_hyb]['extension'] / dict_miRNAex_in_hyb[each_hyb]['total'], decimals=3)
            dict_miRNAex_in_hyb_ratio[each_hyb] = miRNAex_in_hyb_ratio

        if os.path.exists(outfile):
            os.remove(outfile)

        f_out = open(outfile, 'a+')
        with open(self.infile, 'r+') as f1:
            lines = f1.readlines()
            line1st = lines[:1][0]
            f_out.write(line1st.rstrip(
                '\n') + f'\tmiRNA_tot\tmiRNA_ex_tot\tmiRNA_ex_ratio\thyb_tot\thyb_ex_tot\thyb_ex_ratio\n')
            lines_others = lines[1:]
            for line1 in lines_others:
                line2 = line1.rstrip('\n').split('\t')
                miRNA_name = line2[0]
                miRNA_pattern = line2[2]
                Gene_ID = line2[3]
                target_pattern = line2[5]
                striped_target_pattern = ''.join('\\' + x for x in target_pattern.strip('.'))
                left_dot_target_num = re.search(striped_target_pattern, target_pattern).start()
                right_dot_target_num = re.search(striped_target_pattern, target_pattern).end()
                target_seq = line2[4]
                target_seq_without_unpaired = target_seq[left_dot_target_num:right_dot_target_num]
                each_row_hyb_name = f"{miRNA_name}${Gene_ID}${miRNA_pattern}${target_pattern.strip('.')}${target_seq_without_unpaired}"
                if each_row_hyb_name in dict_miRNAex_in_hyb:
                    f_out.write(line1.rstrip(
                        '\n') + f"\t{dict_eachmiRNA_ex[miRNA_name]['total']}\t{dict_eachmiRNA_ex[miRNA_name]['extension']}\t{dict_each_miRNA_extension_ratio[miRNA_name]}\t{dict_miRNAex_in_hyb[each_row_hyb_name]['total']}\t{dict_miRNAex_in_hyb[each_row_hyb_name]['extension']}\t{dict_miRNAex_in_hyb_ratio[each_row_hyb_name]}\n")
                else:
                    print(f"{each_row_hyb_name}  not in dict_miRNAex_in_hyb, check Lu's code.")
        f_out.close()


class Statistic():
    def viennad_rnatype_number(self, viennad_file,microRNA_database): ## In this ViennaRNA file analysis to determine the number of RNA types, I attempted to modify certain genes. If their sequences had a 100% match with microRNA, I relabeled them as microRNA, such as changing MIR17HG to miR-17 microRNA.
        miRNA_seq_name_dict = Database().microRNA_sequence_to_name_database_1st_18nt(input=microRNA_database)
        print(viennad_file)
        with open(viennad_file, 'r+') as f1:
            list_hyb = []
            dict_hyb = {}
            for line1 in f1:
                line2 = line1.strip().split('_')
                if len(line2) == 15:  # microRNA on the left of hybrid, type=mim pref=mim
                    if (len(list_hyb) == 6) and (list_hyb.count('')==1):
                        hyb_each_num = int(list_hyb[0].split('_')[1])

                        RNA1_seq = list_hyb[2].split('\t')[0].strip('-')
                        RNA1_type = 'microRNA' if 'microRNA' in list_hyb[2].split('\t')[1] else list_hyb[2].split('\t')[1].split('#')[-1].split('&')[-1].replace('_mRNA', '')
                        RNA1_name = list_hyb[2].split('\t')[1].split('#')[0]

                        RNA2_seq = list_hyb[3].split('\t')[0].strip('-')
                        RNA2_type = 'microRNA' if 'microRNA' in list_hyb[3].split('\t')[1] else list_hyb[3].split('\t')[1].split('#')[-1].split('&')[-1].replace('_mRNA', '')
                        RNA2_name = list_hyb[3].split('\t')[1].split('#')[0]
                        
                        RNA1_type = 'microRNA' if RNA1_type != 'microRNA' and RNA1_seq[:18] in miRNA_seq_name_dict else RNA1_type
                        RNA2_type = 'microRNA' if RNA2_type != 'microRNA' and RNA2_seq[:18] in miRNA_seq_name_dict else RNA2_type
                        RNA_two_types = '_'.join(sorted([RNA1_type, RNA2_type]))

                        if (RNA1_name == RNA2_name) and ('microRNA' not in [RNA1_type,RNA2_type]):
                            pass ## remove single transcript RNA
                        else:
                            dict_hyb[RNA_two_types] = dict_hyb.get(RNA_two_types,0)
                            dict_hyb[RNA_two_types] += hyb_each_num

                            if RNA_two_types in ['lncRNA_microRNA', 'mRNA_microRNA']:
                                if RNA1_type == 'microRNA':
                                    RC_rna_Seq = Seq(RNA1_seq[1:7]).reverse_complement()
                                    if str(RC_rna_Seq) in RNA2_seq:
                                        dict_hyb['seed_mRNA_lncRNA'] = dict_hyb.get('seed_mRNA_lncRNA', 0) + hyb_each_num
                                elif RNA2_type == 'microRNA':
                                    RC_rna_Seq = Seq(RNA2_seq[1:7]).reverse_complement()
                                    if str(RC_rna_Seq) in RNA1_seq:
                                        dict_hyb['seed_mRNA_lncRNA'] = dict_hyb.get('seed_mRNA_lncRNA', 0) + hyb_each_num
                    list_hyb = []
                list_hyb.append(line1.strip())
        df = pd.DataFrame(list(dict_hyb.items()), columns=['RNA_type', 'Count'])
        df_sorted = df.sort_values(by='Count', ascending=False)
        print(df_sorted)
        df_sorted.to_csv(f'{viennad_file}.RNAtype.csv', index=False)

    def miRNA_count_from_whole_nt(self, input1,microRNA_database):  ## i only count miRNA on the left, right RNA could be any RNA type, like miRNA, lncRNA, rRNA, mRNA etc.
        miRNA_seq_name_dict = Database().microRNA_sequence_to_name_database_whole_nt(input=microRNA_database)
        print(input1)
        dict_miRNA_count = {}
        with open(input1, 'r+') as f1:
            for line1 in f1:
                miRNA_seq = line1.strip()
                if miRNA_seq in miRNA_seq_name_dict:
                    miRNA_name = miRNA_seq_name_dict[line1.strip()]
                    miRNA_name_seq = miRNA_name + '_' + miRNA_seq
                    dict_miRNA_count[miRNA_name_seq] = dict_miRNA_count.get(miRNA_name_seq,0)
                    dict_miRNA_count[miRNA_name_seq] += 1
        df = pd.DataFrame(list(dict_miRNA_count.items()), columns=['miRNA', 'Count'])
        df_sorted = df.sort_values(by='Count', ascending=False)
        print(df_sorted)
        df_sorted.to_csv(f'{input1}.miRNA_count_from_single_reads.csv', index=False)

    def miRNA_count_from_Viennad(self, viennad_file,microRNA_database):  ## i only count miRNA on the left, right RNA could be any RNA type, like miRNA, lncRNA, rRNA, mRNA etc.
        miRNA_seq_name_dict = Database().microRNA_sequence_to_name_database_whole_nt(input=microRNA_database)
        print(viennad_file)
        with open(viennad_file, 'r+') as f1:
            list_hyb = []
            dict_hyb = {}
            for line1 in f1:
                line2 = line1.strip().split('_')
                if len(line2) == 15:  # microRNA on the left of hybrid, type=mim pref=mim
                    if (len(list_hyb) == 6) and (list_hyb.count('') == 1):
                        hyb_each_num = int(list_hyb[0].split('_')[1])
                        RNA1_seq = list_hyb[2].split('\t')[0].strip('-')
                        if RNA1_seq in miRNA_seq_name_dict:
                            RNA1_name = miRNA_seq_name_dict[RNA1_seq]
                            miRNA_name_seq = RNA1_name + '_' + RNA1_seq
                            dict_hyb[miRNA_name_seq] = dict_hyb.get(miRNA_name_seq,0)
                            dict_hyb[miRNA_name_seq] += hyb_each_num
                    list_hyb = []
                list_hyb.append(line1.strip())
        df = pd.DataFrame(list(dict_hyb.items()), columns=['miRNA', 'Count'])
        df_sorted = df.sort_values(by='Count', ascending=False)
        print(df_sorted)
        df_sorted.to_csv(f'{viennad_file}.miRNA_count_from_viennad.csv'.replace('_comp_martquery_0228172604_181_human_unique1_add18S28S_hybrids_ua.viennad',''), index=False)

    def insert_between_miRNA_target(self, input1, str1, str2):
        with open(input1, 'r+') as f1:
            dict1 = {}
            for read in f1:
                try:  # find_sequence != []
                    find_sequence = re.findall(str(str1) + r".{0,200}" + str(str2), read)[0]
                    insert_sequence = (find_sequence[len(str1):-len(str2)])
                    dict1[insert_sequence] = dict1.get(insert_sequence, [0])
                    dict1[insert_sequence][0] += 1
                except:
                    pass
            f1 = pd.DataFrame.from_dict(dict1, orient='index')
            f1.sort_values(by=[0], ascending=False, inplace=True)
            f1.rename(columns={0: f"{str1}_{str2}"}, inplace=True)
            print(f1)
            f1.to_csv(input1.replace('fasta', 'insert.csv').replace('fastq', 'insert.csv'))

    def dG_energy_bar_chart(self, input1):
        dict_dG = {}
        list1 = []
        with open(input1, 'r+') as f1:
            for line1 in f1:
                line2 = line1.strip().split('_')
                if (len(line2) == 15) and (list1 != []):
                    hyb_num1 = int(list1[0].split('_')[1])
                    RNA_type1 = list1[0].split('_')[-4].split('#')[-1].split('&')[-1]
                    dG1 = float(list1[4].split('\t')[1][1:-1])
                    dict_dG[RNA_type1] = dict_dG.get(RNA_type1, {'dG_sum': 0, 'total_hybrids': 0})
                    dict_dG[RNA_type1]['dG_sum'] += dG1 * hyb_num1
                    dict_dG[RNA_type1]['total_hybrids'] += hyb_num1
                    list1 = []
                list1.append(line1.strip())
        f2 = pd.DataFrame.from_dict(dict_dG, orient='index')
        f2['dG_mean'] = (f2['dG_sum'] / f2['total_hybrids'])
        f_dG = f2['dG_mean'].copy()
        f_dG.sort_values(inplace=True)
        f_dG = f_dG.reset_index()
        f_dG.rename(columns={"index": "RNA_type"}, inplace=True)
        f_dG = f_dG.head(20)
        fig = px.bar(f_dG, x='RNA_type', y='dG_mean')
        fig.update_layout(barmode='group', xaxis_tickangle=-45)
        # fig.show()
        fig.write_image(f"{input1}_dG.pdf")

    def rna_type_count(self, input1):
        dict1 = {}
        with open(input1, 'r+') as f1:
            for line1 in f1:
                line2 = line1.strip().split('\t')
                rna_type = line2[0].split('#')[-1].replace('_mRNA', '')
                if 'microRNA' in rna_type:
                    rna_type = 'microRNA'
                if rna_type != '':
                    dict1[rna_type] = dict1.get(rna_type, 0)
                    dict1[rna_type] += int(line2[1])
        f1 = pd.DataFrame.from_dict(dict1, orient='index')
        f1 = f1.head(10)
        print(list(f1[0]))
        print(list(f1.index))
        colors = sns.color_palette('pastel')[0:5]
        plt.pie(x=list(f1[0]), labels=list(f1.index), autopct='%.0f%%', colors=colors)
        plt.title(f'{input1}')
        plt.show()

    def fastq_length_distribution(self, input1):
        with open(input1, 'r+') as f1:
            list1 = []
            dict_read_length = {}
            for line1 in f1:
                if line1[0] == '@':
                    if len(list1) == 4:
                        read_length = len(list1[1])
                        # print(read_length)
                        dict_read_length[read_length] = dict_read_length.get(read_length, 0)
                        dict_read_length[read_length] += 1
                    list1 = []
                    list1.append(line1.strip())
                else:
                    list1.append(line1.strip())
        f2 = pd.DataFrame.from_dict(dict_read_length, orient='index')
        f2.rename(columns={0: input1}, inplace=True)
        f2.reset_index(inplace=True)
        f2.sort_values(by=['index'], ascending=True, inplace=True)
        f2.to_csv(input1.replace('fastq', 'csv').replace('fq', 'csv'), index=False)
        sns.set(style="darkgrid", rc={'figure.figsize': (20, 10)})
        sns.set(font_scale=2)
        sns.lineplot(x='index', y=f'{input1}', data=f2)
        plt.savefig(f"{input1}.pdf")

    def fasta_reads_number(self, input1, min_len, max_len):
        with open(input1, 'r+') as f1:
            total_reads_num = 0
            for line1 in f1:
                if line1[0] != '>':
                    read_length = len(line1.strip())
                    if int(min_len) <= read_length <= int(max_len):
                        total_reads_num += 1
            print(f"{input1} total reads number is {total_reads_num}")

    def miR16_sfGFP_count(self,input1):
        with open(input1, 'r+') as f1:
            total_num = 0
            for line1 in f1:
                if ('miR-16-5p' in line1) and ('sfGFP' in line1):
                    count_num = int(line1.split('\t')[0].split('_')[1])
                    total_num += count_num
        print(f'total miR-16-5p:sfGFP number in {input1} is {total_num}\n')

    def fasta_length_distribution(self, input1):
        with open(input1, 'r+') as f1:
            dict_read_length = {}
            for line1 in f1:
                if line1[0] != '>':
                    read_length = len(line1.strip())
                    dict_read_length[read_length] = dict_read_length.get(read_length, 0)
                    dict_read_length[read_length] += 1
        f2 = pd.DataFrame.from_dict(dict_read_length, orient='index')
        f2.reset_index(inplace=True)
        f2.rename(columns={0: input1, 'index': 'length'}, inplace=True)
        f2.sort_values(by=['length'], ascending=True, inplace=True)
        f2.to_csv(input1.replace('fasta', 'csv'), index=False)
        print(f2)
        sns.set(style="darkgrid", rc={'figure.figsize': (20, 10)})
        sns.set(font_scale=2)
        sns.lineplot(x='length', y=f'{input1}', data=f2)
        plt.savefig(f"{input1}.pdf")

    def hyb_number(self, input1, microRNA_database):
        miRNA_database1 = Database().microRNA_database(input=microRNA_database)
        total_hybrids = 0
        miRNA_hybrids = 0
        miRNA_longer_hybrids = 0
        miRNA_mRNA_longer_hybrids_left = 0
        miRNA_mRNA_longer_hybrids_right = 0
        miRNA_mRNA_longer_hybrids = 0
        miRNA_mRNA_longer_hybrids_seed = 0  # 2nd to 7th match
        with open(input1, 'r+') as f1:
            for line1 in f1:
                line2 = line1.strip().split('\t')
                hybrid_num = int(line2[0].split('_')[1])
                hybrid_sequence = line2[1]
                RNA1 = line2[3]
                RNA2 = line2[9]
                total_hybrids += hybrid_num
                if ('microRNA' in line1) and ('mRNA' in line1) and ('MIMAT' in line1):
                    miRNA_hybrids += hybrid_num
                    if len(hybrid_sequence) >= 50:
                        miRNA_longer_hybrids += hybrid_num
                        if ('mRNA_mRNA' in line1):
                            miRNA_mRNA_longer_hybrids += hybrid_num
                            if 'microRNA' in RNA1:
                                miRNA_mRNA_longer_hybrids_left += hybrid_num
                                miRNA_name = RNA1.split('_')[2]
                                miRNA_seed = miRNA_database1[miRNA_name][1:7]
                                miRNA_seed_RC = str(Seq(miRNA_seed).reverse_complement())
                                if miRNA_seed_RC in hybrid_sequence[int(line2[10]) - 1:int(line2[11])]:
                                    miRNA_mRNA_longer_hybrids_seed += hybrid_num
                            if 'microRNA' in RNA2:
                                miRNA_mRNA_longer_hybrids_right += hybrid_num
                                miRNA_name = RNA2.split('_')[2]
                                miRNA_seed = miRNA_database1[miRNA_name][1:7]
                                miRNA_seed_RC = str(Seq(miRNA_seed).reverse_complement())
                                if miRNA_seed_RC in hybrid_sequence[int(line2[4]) - 1:int(line2[5])]:
                                    miRNA_mRNA_longer_hybrids_seed += hybrid_num
            print(f"total hybrids are {total_hybrids}")
            print(f"microRNA hybrids are {miRNA_hybrids}")
            print(f"! microRNA longer hybrids are {miRNA_longer_hybrids}")
            print(f"left microRNA-mRNA longer hybrids are {miRNA_mRNA_longer_hybrids_left}")
            print(f"right microRNA-mRNA longer hybrids are {miRNA_mRNA_longer_hybrids_right}")
            print(f"! microRNA-mRNA longer hybrids are {miRNA_mRNA_longer_hybrids}")
            print(f"! microRNA-mRNA longer hybrids with seed 2nd-7th match are {miRNA_mRNA_longer_hybrids_seed}")

    def hyb_length(self, hyb_file):
        with open(hyb_file, 'r+') as f1:
            dict_hyb_length = {}
            dict_microRNA_hyb_length = {}
            total_reads = 0
            for line1 in f1:
                line2 = line1.strip().split('\t')
                hyb_num = int(line2[0].split('_')[1])
                hyb_sequence = line2[1]
                hyb_length1 = len(hyb_sequence)
                total_reads += hyb_num
                RNA_type1 = line2[3].split('#')[-1].split('_')[0].split('&')[-1]
                RNA_type2 = line2[9].split('#')[-1].split('_')[0].split('&')[-1]
                RNA_type = '_'.join(sorted([RNA_type1, RNA_type2]))
                RNA_type = re.sub('MIMAT\d{7}', 'microRNA', RNA_type)
                dict_hyb_length[hyb_length1] = dict_hyb_length.get(hyb_length1, {})
                dict_hyb_length[hyb_length1][RNA_type] = dict_hyb_length[hyb_length1].get(RNA_type, 0)
                dict_hyb_length[hyb_length1][RNA_type] += hyb_num
                if 'microRNA' in RNA_type:
                    dict_microRNA_hyb_length[hyb_length1] = dict_microRNA_hyb_length.get(hyb_length1, {})
                    dict_microRNA_hyb_length[hyb_length1][RNA_type] = dict_microRNA_hyb_length[hyb_length1].get(
                        RNA_type, 0)
                    dict_microRNA_hyb_length[hyb_length1][RNA_type] += hyb_num
        ###### RNA-RNA distribution
        f_data = pd.DataFrame.from_dict(dict_hyb_length, orient='index')
        f_data.index = f_data.index.astype(int)
        f_data.sort_index(inplace=True)
        f_data.to_csv(f'{hyb_file}.csv')
        f_sum = f_data.sum()
        f_sum.sort_values(ascending=False, inplace=True)
        top10_RNA_RNA = f_sum.head(10)
        print(f"The top 10 RNA-RNA hybrids types are \n{top10_RNA_RNA}")
        top10_RNA_RNA_list1 = list(top10_RNA_RNA.index)
        f_data = (f_data[top10_RNA_RNA_list1])
        sns.set(style="darkgrid", rc={'figure.figsize': (30, 15)})
        sns.set(font_scale=2)
        f_data.plot(kind='bar', stacked=True, ).set(title=f'{hyb_file}')
        # plt.show()
        plt.savefig(f'{hyb_file}_RNA_RNA.pdf')
        ###### microRNA-RNA distribution
        f_data = pd.DataFrame.from_dict(dict_microRNA_hyb_length, orient='index')
        f_data.index = f_data.index.astype(int)
        f_data.sort_index(inplace=True)
        # f_data.to_csv(f'{hyb_file}.csv')
        f_sum = f_data.sum()
        f_sum.sort_values(ascending=False, inplace=True)
        top10_RNA_RNA = f_sum.head(10)
        print(f"The top 10 microRNA-RNA hybrids types are \n{top10_RNA_RNA}")
        top10_RNA_RNA_list1 = list(top10_RNA_RNA.index)
        f_data = (f_data[top10_RNA_RNA_list1])
        sns.set(style="darkgrid", rc={'figure.figsize': (30, 15)})
        sns.set(font_scale=2)
        f_data.plot(kind='bar', stacked=True, ).set(title=f'{hyb_file}')
        # plt.show()
        plt.savefig(f'{hyb_file}_microRNA_RNA.pdf')

    def single_miRNA_target_average_length(self, hyb_file, miRNA_seq, target_seq):
        with open(hyb_file, 'r+') as f1:
            dict_CLASH_length = {}
            for line1 in f1:
                line2 = line1.strip().split('\t')
                hyb_num = int(line2[0].split('_')[1])
                hyb_sequence = line2[1]
                hyb_length1 = (len(hyb_sequence))
                if (miRNA_seq in line1) and (target_seq in line1):
                    dict_CLASH_length[hyb_length1] = dict_CLASH_length.get(hyb_length1, 0)
                    dict_CLASH_length[hyb_length1] += hyb_num
        # dict_CLASH_length = collections.OrderedDict(sorted(dict_CLASH_length.items()))
        dict_CLASH_length = pd.DataFrame.from_dict(dict_CLASH_length, orient='index')
        sns.barplot(x=dict_CLASH_length.index, y=dict_CLASH_length[0])
        plt.xticks(rotation=45, fontsize=5)
        print(dict_CLASH_length)
        plt.show()

    def miRNA_CSV_merge(self, string1):
        print(sorted(glob.glob(f'*{string1}*.csv')))
        f_all = pd.DataFrame()
        for file1 in sorted(glob.glob(f'*{string1}*.csv')):
            f = pd.read_csv(file1, index_col=[0])
            f.rename(columns={f.columns[0]: file1}, inplace=True)
            f_all = pd.concat([f_all, f], axis=1)
        f_all.to_csv(f'{string1}_Merged_miRNA_raw_count.csv')
        f_all = f_all / f_all.sum() * 1000000
        f_all.to_csv(f'{string1}_Merged_miRNA_normalized_by_total_reads.csv')

    def table_merge(self, string1):
        print(sorted(glob.glob(f'*{string1}*')))
        f_all = pd.DataFrame()
        for file1 in sorted(glob.glob(f'*{string1}*')):
            f = pd.read_table(file1, index_col=[0])
            print(f)
            f.rename(columns={f.columns[0]: file1}, inplace=True)
            f_all = pd.concat([f_all, f], axis=1)
        f_all.to_csv(f'{string1}_Merged.count')


    def miRNA_CSV_two_index_merge(self, string1):
        print(sorted(glob.glob(f"*{string1}*.csv")))
        f_all = pd.DataFrame()
        for file1 in sorted(glob.glob(f"*{string1}*.csv")):
            f1 = pd.read_csv(file1, index_col=[0, 1])
            f1.rename(columns={f1.columns[0]: file1}, inplace=True)
            f_all = pd.concat([f_all, f1], axis=1)
        f_all.to_csv(f'{string1}_Merged_miRNA_raw_count.csv')
        f_all = f_all / f_all.sum() * 1000000
        f_all.to_csv(f'{string1}_Merged_miRNA_normalized_by_total_reads.csv')

    def gene_table_two_index_merge(self, string1):
        list1 = sorted(glob.glob(f"*{string1}*"))
        print(list1)
        f_all = pd.DataFrame()
        for file1 in list1:
            print(file1)
            f1 = pd.read_table(file1, index_col=[0, 1], header=None)
            f1.rename(columns={2: f"{file1}"}, inplace=True)
            f_all = pd.concat([f_all, f1], axis=1)
        f_all.to_csv(f'{string1}_Merged_gene_RAWcount.csv')

    def all_miRNA_isoform_number_table_1st_18th(self, input_file, microRNA_database, longest_isoform_length):
        miRNA_seq_name_dict = Database().microRNA_sequence_to_name_database_1st_18nt(input=microRNA_database)
        output_file = input_file + f'_miRNA_isoform_number_table.csv'
        miRNA_isoform_dict1 = dict()
        miRNA_isoform_dict1['count'] = miRNA_isoform_dict1.get('count', {})
        with open(input_file) as f1:
            for line1 in f1:
                miRNA_sequence_from_reads = line1.strip()
                if (line1[:18] in miRNA_seq_name_dict) and (
                        len(line1.strip()) <= int(longest_isoform_length)):  ## miRNA identification
                    each_miRNA_name1_rawseq = miRNA_seq_name_dict[line1[:18]]
                    each_miRNA_name1 = each_miRNA_name1_rawseq.split('&')[0]
                    each_miRNA_rawseq1 = each_miRNA_name1_rawseq.split('&')[1]
                    if each_miRNA_rawseq1 == miRNA_sequence_from_reads:  ## raw miRNA
                        each_miRNA_name1 = f"{each_miRNA_name1}_Raw"
                    miRNA_name_isoformseq_tuble = (f"{each_miRNA_name1}", f"{miRNA_sequence_from_reads}")
                    miRNA_isoform_dict1['count'][miRNA_name_isoformseq_tuble] = miRNA_isoform_dict1['count'].get(
                        miRNA_name_isoformseq_tuble, 0)
                    miRNA_isoform_dict1['count'][miRNA_name_isoformseq_tuble] += 1
        f1 = pd.DataFrame.from_dict(miRNA_isoform_dict1)
        f1.to_csv(output_file)
        return f1

    def all_miRNA_isoform_kind_count_1st_18th(self, input_file, microRNA_database, longest_isoform_length):
        miRNA_seq_name_dict = Database().microRNA_sequence_to_name_database_1st_18nt(input=microRNA_database)
        dict_miRNA_isoform_kind_sequence = {}  ## miRNA kind name and sequence
        dict_miRNA_isoform_kind_number = {}  ## miRNA kind name and number
        output_file = input_file + '_miRNA_isoform_kind.csv'
        with open(input_file) as f1:
            for line1 in f1:
                if len(line1.strip()) <= int(longest_isoform_length):
                    read_seq = line1.strip()[:18]
                    if read_seq in miRNA_seq_name_dict:
                        read_name_Rawseq = miRNA_seq_name_dict[read_seq]
                        read_name = read_name_Rawseq.split('&')[0]
                        dict_miRNA_isoform_kind_sequence[read_name] = dict_miRNA_isoform_kind_sequence.get(read_name,
                                                                                                           set())
                        dict_miRNA_isoform_kind_sequence[read_name].add(line1.strip())
        for miRNA_name, miRNA_different_reads in dict_miRNA_isoform_kind_sequence.items():
            dict_miRNA_isoform_kind_number[miRNA_name] = len(miRNA_different_reads)
        dict_miRNA_isoform_kind_number = collections.OrderedDict(sorted(dict_miRNA_isoform_kind_number.items()))
        dict_miRNA_isoform_kind_number = pd.DataFrame.from_dict(dict_miRNA_isoform_kind_number, orient='index')
        dict_miRNA_isoform_kind_number.rename(columns={0: input_file}, inplace=True)
        dict_miRNA_isoform_kind_number.to_csv(output_file)
        print('Done!')
        return dict_miRNA_isoform_kind_number

    def all_miRNA_length_distribution_1st_18th(self, input_file, microRNA_database, longest_isoform_length):
        miRNA_seq_name_dict = Database().microRNA_sequence_to_name_database_1st_18nt(input=microRNA_database)
        dict_miRNA_length_abundance = {}  ## miRNA length abundance
        dict_miRNA_length_abundance['miRNA_length'] = dict_miRNA_length_abundance.get('miRNA_length', {})
        output_file = input_file + '_miRNA_length_distribution.csv'
        with open(input_file) as f1:
            for line1 in f1:
                miRNA_isoform_length = len(line1.strip())
                if miRNA_isoform_length <= int(longest_isoform_length):
                    read_seq = line1.strip()[:18]
                    if read_seq in miRNA_seq_name_dict:
                        read_name_Rawseq = miRNA_seq_name_dict[read_seq]
                        miRNA_name_length_tuble = (f"{read_name_Rawseq}", f"{miRNA_isoform_length}")
                        dict_miRNA_length_abundance['miRNA_length'][miRNA_name_length_tuble] = \
                        dict_miRNA_length_abundance['miRNA_length'].get(
                            miRNA_name_length_tuble, 0)
                        dict_miRNA_length_abundance['miRNA_length'][miRNA_name_length_tuble] += 1
        f1 = pd.DataFrame.from_dict(dict_miRNA_length_abundance)
        f1.to_csv(output_file)
        return f1

    def miRNA_length_distribution(self, input_file, microRNA_database, miRNA_sequence,
                                  primiRNA_sequence='xxx'):  ## miRNA isoform lengrh distribution range  from 16nt to 30nt
        miRNA_seq_name_dict_16nt = Database().microRNA_sequence_to_name_database_16nt(microRNA_database)
        miRNA_sequence_16nt = miRNA_sequence[:16]
        if miRNA_sequence_16nt in miRNA_seq_name_dict_16nt:
            print(f"miRNA sequence is {miRNA_sequence}")
            print(f"miRNA name is {miRNA_seq_name_dict_16nt[miRNA_sequence_16nt]}")
        elif miRNA_sequence_16nt not in miRNA_seq_name_dict_16nt:
            print(f"miRNA sequence is not exist in database")

        output_file = input_file + f'_length_distribution_{miRNA_seq_name_dict_16nt[miRNA_sequence_16nt]}.csv'
        miRNA_length_distribution_dict1 = dict()
        with open(input_file) as f1:
            for line1 in f1:
                read_sequence = line1.strip()
                read_length = len(read_sequence)
                read_length_name = f"{read_length}_nt"
                if (read_length <= 30) and (miRNA_sequence_16nt == line1[:16]):  ## miRNA identification
                    miRNA_length_distribution_dict1[f'{read_length_name}'] = miRNA_length_distribution_dict1.get(
                        f'{read_length_name}', {"match": 0, "unmatch": 0})
                    if read_sequence in primiRNA_sequence:
                        miRNA_length_distribution_dict1[f'{read_length_name}']['match'] += 1
                    elif read_sequence not in primiRNA_sequence:
                        miRNA_length_distribution_dict1[f'{read_length_name}']['unmatch'] += 1
        miRNA_length_distribution_dict1 = collections.OrderedDict(sorted(miRNA_length_distribution_dict1.items()))
        f = pd.DataFrame.from_dict(miRNA_length_distribution_dict1)
        f = f.T
        print(f"miRNA total number is {f.sum().sum()}")
        f.to_csv(output_file)
        print('done!')
        return miRNA_length_distribution_dict1

    def miRNA_length_distribution_merge(self):
        print(sorted(glob.glob('*length_distribution*csv')))
        f_all = pd.DataFrame()
        for file1 in sorted(glob.glob('*length_distribution*csv')):
            f = pd.read_csv(file1, index_col=[0])
            f.rename(columns={'match': f'match_{file1}'}, inplace=True)
            f.rename(columns={'unmatch': f'unmatch_{file1}'}, inplace=True)
            f_all = pd.concat([f_all, f], axis=1)
        f_all.to_csv('Merged_miRNA_length_distribution.csv')

    def RNA_RNA_hybrid_distribution(self, hyb_file):  # This hyb including all different kinds of RNA-RNA
        with open(hyb_file, 'r+') as f1:
            dict1 = {}
            for line1 in f1:
                line2 = line1.strip().split('\t')
                RNA1 = line2[3].split('_')[2].split('#')[0]
                RNA2 = line2[9].split('_')[2].split('#')[0]
                hybrid_name = '&'.join(sorted([RNA1, RNA2]))
                number1 = int(line2[0].split('_')[1])
                dict1[hybrid_name] = dict1.get(hybrid_name, 0)
                dict1[hybrid_name] += number1
        f2 = pd.DataFrame.from_dict(dict1, orient='index')
        f2.sort_values(by=[0], inplace=True, ascending=False)
        hyb_file_name = hyb_file.split('_')[0]
        f2.rename(columns={0: hyb_file_name}, inplace=True)
        f2.reset_index(inplace=True)
        f2 = (f2.head(50))
        print(f2)
        fig = px.bar(f2, x='index', y=f'{hyb_file_name}')
        fig.update_layout(barmode='group', xaxis_tickangle=-45)
        fig.show()
        fig.write_image(f"{hyb_file_name}.pdf")

    def Viennad_table_element_region_statistic(self, input_table):
        f1 = pd.read_table(input_table, header=0)
        f1 = f1[f1['element_region'].notna()]
        dict_region = {}
        total_element_region_abundance = 0
        for index1 in f1.index:
            element_region_abundance = f1.loc[index1]['abundance']
            element_region_type = f1.loc[index1]['element_region']
            dict_region[element_region_type] = dict_region.get(element_region_type, 0)
            dict_region[element_region_type] += element_region_abundance
            total_element_region_abundance += element_region_abundance
        for element1, abundance1 in dict_region.items():
            print(f"{element1}\t\t{abundance1}\t{round(abundance1 / total_element_region_abundance, 3)}")

    def miRNA_abundance_in_hybrids(self, hyb_file):  # focus on miRNA start with the 1st nt in each hybrid
        with open(hyb_file, 'r+') as f1:
            dict_miRNA = {}
            for line1 in f1:
                line2 = line1.strip().split('\t')
                hybrid_number = int(line2[0].split('_')[1])
                if ('microRNA' in line2[3]) and ('1' == line2[4]):  ## microRNA on the left and start with the 1st nt
                    if ('18S' in line1) or ('28S' in line1):  # microRNA ligased with rRNA
                        miRNA_name = line2[3]
                        dict_miRNA[miRNA_name] = dict_miRNA.get(miRNA_name, 0)
                        dict_miRNA[miRNA_name] += hybrid_number
        dict_miRNA = pd.DataFrame.from_dict(dict_miRNA, orient='index')
        dict_miRNA.sort_values(by=[0], inplace=True, ascending=False)
        dict_miRNA.rename(columns={0: hyb_file}, inplace=True)
        dict_miRNA.to_csv(hyb_file + '_miRNA_abundance_in_hyb.csv')
        return dict_miRNA

    def Bowtie2_SAM_singleTranscript_shortReads_statistic(self,input1): # 22bp, 29bp
        samFile = input1
        with open(f"{samFile}", "r") as f1:
            dict_reads = {}  # for each reads, summarize all mapped RNA type
            for line1 in f1:
                if line1.startswith('@'):
                    continue  # Skip header lines
                # Process alignment data lines
                columns = line1.split('\t')
                read_name = columns[0]
                mapped_RNA_type = 'microRNA' if 'microRNA' in columns[2] else columns[2].split("#")[-1].split("_")[0]
                # Initialize key with empty set if not present, then add mapped_RNA_type to set for read_name
                if mapped_RNA_type != '*':
                    dict_reads.setdefault(read_name, set()).add(mapped_RNA_type)

        # Initialize a dictionary to store the summary of RNA types
        dict_RNA_type_summarize = {}

        # Iterate over each read and its mapped RNA type
        for read_name, mapped_rna_type in dict_reads.items():
            # Check for microRNA
            if "microRNA" in mapped_rna_type: ## some times microRNA and lincRNA belong to one reads, if microRNA be mapped, it is microRNA but not other RNA like lincRNA
                dict_RNA_type_summarize["microRNA"] = dict_RNA_type_summarize.get("microRNA", 0) + 1
            # Check for rRNA types RN18S1 and RN28S1
            elif any(rRNA in mapped_rna_type for rRNA in ["RN18S1&rRNA", "RN28S1&rRNA"]):
                dict_RNA_type_summarize["rRNA"] = dict_RNA_type_summarize.get("rRNA", 0) + 1
            # Handle other single RNA types
            elif len(mapped_rna_type) == 1:
                single_RNA_type = next(iter(mapped_rna_type))  # Extract the single RNA type
                dict_RNA_type_summarize[single_RNA_type] = dict_RNA_type_summarize.get(single_RNA_type, 0) + 1

        # Assuming dict_RNA_type_summarize is your existing dictionary
        table = pd.DataFrame.from_dict(dict_RNA_type_summarize, orient='index',
                                       columns=[f"{samFile.replace('.sam', '')}"])

        table.to_csv(f"{samFile.replace('sam', 'csv')}")

    def Bowtie2_SAM_singleTranscript_longReads_statistic(self,input1): # 44bp
        samFile = input1
        with open(f"{samFile}", "r") as f1:
            dict_reads = {}  # for each reads, summarize all mapped RNA type
            for line1 in f1:
                if line1.startswith('@'):
                    continue  # Skip header lines
                # Process alignment data lines
                columns = line1.split('\t')
                read_name = columns[0]
                mapped_RNA_type = 'microRNA' if 'microRNA' in columns[2] else columns[2].split("#")[-1].split("_")[0]
                # Initialize key with empty set if not present, then add mapped_RNA_type to set for read_name
                if mapped_RNA_type != '*':
                    dict_reads.setdefault(read_name, set()).add(mapped_RNA_type)

        # Initialize a dictionary to store the summary of RNA types
        dict_RNA_type_summarize = {}

        # Iterate over each read and its mapped RNA type
        for read_name, mapped_rna_type in dict_reads.items():
            if (len(mapped_rna_type) == 1) :
                single_RNA_type = next(iter(mapped_rna_type)).replace("RN28S1&", "").replace("RN18S1&",
                                                                                             "")  # Extract the single RNA type
                dict_RNA_type_summarize[single_RNA_type] = dict_RNA_type_summarize.get(single_RNA_type, 0) + 1

        # Assuming dict_RNA_type_summarize is your existing dictionary
        table = pd.DataFrame.from_dict(dict_RNA_type_summarize, orient='index',
                                       columns=[f"{samFile.replace('.sam', '')}"])
        table.to_csv(f"{samFile.replace('sam', 'csv')}")
class FASTQ():
    def deduplicate(self, input_file):
        output_file = input_file.replace('.fastq', '_deduplicated.fastq').replace('.fq', '_deduplicated.fa')
        f2 = open(output_file, 'w+')
        with open(input_file, 'r+') as f1:
            set1 = set()
            list1 = []
            num1 = 0
            for line1 in f1:
                if (line1[0] == '@') and (list1 != []) and ('UMI' in line1):
                    UMI1 = list1[0].split('UMI:')[1]
                    seq1 = list1[1]
                    raw_seq1 = seq1 + '_' + UMI1
                    if raw_seq1 not in set1:
                        set1.add(raw_seq1)
                        num1 += 1
                        f2.write(f">{num1}\n{seq1}\n")
                    list1 = []
                list1.append(line1.strip())
        f2.close()

    def fastq_each_read_count(self, input1):
        with open(input1, 'r+') as f1:
            list1 = []
            dict_read_count = {}
            for line1 in f1:
                if line1[0] == '@':
                    if len(list1) == 4:
                        read_seq = (list1[1])
                        dict_read_count[read_seq] = dict_read_count.get(read_seq, 0)
                        dict_read_count[read_seq] += 1
                    list1 = []
                    list1.append(line1.strip())
                else:
                    list1.append(line1.strip())
        df = pd.DataFrame(list(dict_read_count.items()), columns=['miRNA_sequence', 'Count'])
        df.to_csv(input1.replace('fastq','each_read.csv').replace('fq','each_read.csv'),index=False)
        print(df)

class FASTA():
    def reverse_complementary_sequence(self, input):
        with open(input, 'r+') as f1:
            for line1 in f1:
                if line1[0] == '>':
                    print(line1.strip())
                else:
                    read_sequence = Seq(line1.strip().upper())
                    read_sequence_RC = read_sequence.reverse_complement()
                    print(read_sequence_RC)

    def concatenated_reads(self, input1, fragment_length1):
        with open(input1, 'r+') as f1:
            for line1 in f1:
                if line1[0] == '>':
                    name1 = line1.strip()
                else:
                    read1 = line1.strip()
                    for position1 in range(len(read1)):
                        fragment_seq1 = read1[position1:position1 + fragment_length1]
                        if (len(fragment_seq1) == fragment_length1) and (read1.count(fragment_seq1) > 1):
                            print(f'{name1}\n{read1}')
                            break


class SAM():
    def bowtie2_count(self, input1):
        dict_RNA = {}
        with open(input1, 'r+') as f1:
            for line1 in f1:
                line2 = line1.strip().split('\t')
                if (line2[2] != '*') and len(line2) > 3:
                    RNA_name = line2[2]
                    dict_RNA[RNA_name] = dict_RNA.get(RNA_name, 0)
                    dict_RNA[RNA_name] += 1
        dict2_RNA = dict(sorted(dict_RNA.items(), key=lambda x: x[1], reverse=True))
        with open(input1.replace('sam', 'csv'), 'w+') as f2:
            f2.write(f"Name,{input1}\n")
            for x, y in dict2_RNA.items():
                f2.write(f"{x},{y}\n")

    def bowtie2_hOH7_single_Transcripts_count(self, input1):
        print("Criteria: unmapped nucleotides of each read does not exceed 5!\n")
        dict_RNA = {}
        with open(input1, 'r+') as f1:
            for line1 in f1:
                line2 = line1.strip().split('\t')
                if (line2[2] != '*') and len(line2) >= 9:
                    gene_info = line2[2]
                    gene_id = gene_info.split('_')[0]
                    gene_name = gene_info.split('_')[2].split('#')[0]
                    gene_type = gene_info.split('_')[2].split('#')[-1]
                    gene_id_name_type = gene_id + '_' + gene_name + '_' + gene_type
                    gene_length = len(line2[9])
                    gene_mapped_length = sum(int(num1) for num1 in re.findall(r'(\d+)M', line2[5]))
                    if gene_length - gene_mapped_length <= 5:
                        dict_RNA[gene_id_name_type] = dict_RNA.get(gene_id_name_type, 0)
                        dict_RNA[gene_id_name_type] += 1
        dict2_RNA = dict(sorted(dict_RNA.items(), key=lambda x: x[1], reverse=True))
        with open(input1.replace('.sam', 'single_Transcripts.csv'), 'w+') as f2:
            f2.write(f"Name,{input1}\n")
            for x, y in dict2_RNA.items():
                f2.write(f"{x},{y}\n")

    def Hisat2_SAM_gene_Table(self, input1):
        dict_gene = {}
        with open(input1, 'r+') as f1:
            for line1 in f1:
                line2 = line1.strip().split('\t')
                if (line2[2] != '*') and (line1[0] != '@'):
                    gene_info = line2[2]
                    dict_gene[gene_info] = dict_gene.get(gene_info, 0)
                    dict_gene[gene_info] += 1
        df = pd.DataFrame(list(dict_gene.items()), columns=['gene', 'count'])
        df.to_csv(f"{input1.replace('sam', 'csv')}")
        print('Done!')

    def hisat2_SAM_mRNA_count(self, input1):
        total_read_num = 0
        mapped_read_num = 0
        mRNA_read_num = 0
        dict_RNA_type_30bp = {}
        dict_RNA_type_18bp = {}
        with open(input1, 'r+') as f1:
            for line1 in f1:
                line2 = line1.strip().split('\t')
                if len(line2) > 9:
                    mapped_name = line2[2]
                    mapped_read = line2[9]
                    if (len(mapped_read) >= 18):
                        total_read_num += 1
                    if (mapped_name != '*') and (len(mapped_read) >= 18):
                        mapped_read_num += 1
                    if ('RNA' in mapped_name) and (len(mapped_read) >= 30):
                        if 'mRNA' in mapped_name:
                            RNA_type = mapped_name.split('#')[-1].replace('RN18S1&', '').replace('RN28S1&', '').replace(
                                'Mt&', '')
                        elif 'microRNA' in mapped_name:
                            RNA_type = 'microRNA'
                        dict_RNA_type_30bp[RNA_type] = dict_RNA_type_30bp.get(RNA_type, 0)
                        dict_RNA_type_30bp[RNA_type] += 1
                    if ('RNA' in mapped_name) and (30 > len(mapped_read) >= 18):
                        if 'mRNA' in mapped_name:
                            RNA_type = mapped_name.split('#')[-1].replace('RN18S1&', '').replace('RN28S1&', '').replace(
                                'Mt&', '')
                        elif 'microRNA' in mapped_name:
                            RNA_type = 'microRNA'
                        dict_RNA_type_18bp[RNA_type] = dict_RNA_type_18bp.get(RNA_type, 0)
                        dict_RNA_type_18bp[RNA_type] += 1
        print(f"Total reads >= 18nt are {total_read_num}")
        print(f"Mapped reads >= 18nt are {mapped_read_num}")
        sorted_dictionary = sorted(dict_RNA_type_30bp.items(), key=lambda x: x[1], reverse=True)
        # Convert the sorted list of tuples back into a dictionary
        sorted_dictionary = dict(sorted_dictionary)
        for RNA_type_name, number1 in sorted_dictionary.items():
            print(f"Mapped {RNA_type_name.split('_')[0]} reads >= 30nt are {number1}")
        sorted_dictionary = sorted(dict_RNA_type_18bp.items(), key=lambda x: x[1], reverse=True)
        # Convert the sorted list of tuples back into a dictionary
        sorted_dictionary = dict(sorted_dictionary)
        for RNA_type_name, number1 in sorted_dictionary.items():
            print(f"Mapped {RNA_type_name.split('_')[0]} reads < 30nt and >=18nt are {number1}")


class Target_analysis:
    def __init__(self, Deseq2_file, baseMean_threshold):
        self.Deseq2_file_name = Deseq2_file
        self.baseMean_threshold = baseMean_threshold
        f1 = pd.read_csv(Deseq2_file, index_col=[0], header=0)
        f1 = f1.groupby(f1.index).first()  ## deduplicate index
        f1.dropna(inplace=True, axis='index')
        f1 = f1[f1['baseMean'] > int(baseMean_threshold)]
        self.Deseq2_file = f1

    def predicted_targets(self, targets_file):
        targets_pd = pd.read_table(f'{targets_file}', index_col=[0], header=0)
        target_list = (targets_pd.index)
        return target_list  ## whole targets from targetScan, ~10% targets are not in Htseq-count output.Thereby,
        # intersection between mRNA from Htseq-count and targets from TargetScam are necessary

    def output_target_Nontarget_list(self, miR_all_target_file, conserved_target_file):  # list file,
        mRNA_list = set(self.Deseq2_file.index)
        all_target_list_intersection = set(self.predicted_targets(miR_all_target_file)) & mRNA_list
        non_target_list = (mRNA_list - all_target_list_intersection)
        conserved_target_list_intersection = set(self.predicted_targets(conserved_target_file)) & mRNA_list
        return all_target_list_intersection, \
            conserved_target_list_intersection, \
            non_target_list  # target_list_intersection is tragets name list

    def output_target_Nontarget_list_CLASH_targetScan(self, miR_all_target_file, conserved_target_file, CLASH_file,
                                                      CLASH_targetScan_file):  # list file,
        mRNA_list = set(self.Deseq2_file.index)
        all_target_list_intersection = set(self.predicted_targets(miR_all_target_file)) & mRNA_list
        conserved_target_list_intersection = set(self.predicted_targets(conserved_target_file)) & mRNA_list
        CLASH_target_list_intersection = set(self.predicted_targets(CLASH_file)) & mRNA_list
        CLASH_targetScan_target_list_intersection = set(self.predicted_targets(CLASH_targetScan_file)) & mRNA_list
        non_target_list = (mRNA_list - all_target_list_intersection - CLASH_target_list_intersection)
        return all_target_list_intersection, \
            conserved_target_list_intersection, \
            CLASH_target_list_intersection, \
            CLASH_targetScan_target_list_intersection, \
            non_target_list  # target_list_intersection is tragets name list

    def table_target_Nontraget(self, all_target_list, conserved_target_list, Nontarget_list):
        all_targets = list(all_target_list)
        conserved_target = list(conserved_target_list)
        Nontarget = list(Nontarget_list)
        all_gene_pd = self.Deseq2_file
        all_gene_length = len(all_gene_pd.index)  ## all gene target matrix, sorted by values
        print(f'all gene length: {all_gene_length}')
        try:
            all_gene_pd['order'] = np.arange(0, 1, 1 / all_gene_length)
        except:
            all_gene_pd['order'] = np.arange(0, 1, 1 / all_gene_length)[:-1]

        all_target_pd = (self.Deseq2_file.loc[all_targets])  ## all target matrix, sort by value
        all_target_pd = all_target_pd.sort_values(by=['log2FoldChange'])
        target_length = len(all_target_pd.index)
        print(f'all target length: {target_length}')
        try:
            all_target_pd['order'] = np.arange(0, 1, 1 / target_length)
        except:
            all_target_pd['order'] = np.arange(0, 1, (1 / (target_length)))[:-1]

        conserved_target_pd = (self.Deseq2_file.loc[conserved_target])  ## conserved target matrix, sort by value
        conserved_target_pd = conserved_target_pd.sort_values(by=['log2FoldChange'])
        conserved_target_length = len(conserved_target_pd.index)
        print(f'conserved target length: {conserved_target_length}')
        try:
            conserved_target_pd['order'] = np.arange(0, 1, 1 / conserved_target_length)
        except:
            conserved_target_pd['order'] = np.arange(0, 1, 1 / conserved_target_length)[:-1]

        Nontarget_pd = (self.Deseq2_file.loc[Nontarget])  ## Nontarget matrix, sort by value
        Nontarget_pd = (Nontarget_pd.sort_values(by=['log2FoldChange']))
        Nontarget_length = len(Nontarget_pd.index)
        print(f'Nontarget length: {Nontarget_length}')
        try:
            Nontarget_pd['order'] = np.arange(0, 1, 1 / Nontarget_length)
        except:
            Nontarget_pd['order'] = np.arange(0, 1, (1 / (Nontarget_length)))[:-1]
        return all_gene_pd, all_target_pd, conserved_target_pd, Nontarget_pd

    def table_target_Nontraget_CLASH(self, all_target_list, conserved_target_list, Nontarget_list, CLASH_list,
                                     CLASH_targetScan_list):
        all_targets = list(all_target_list)
        conserved_target = list(conserved_target_list)
        CLASH_list_target = list(CLASH_list)
        CLASH_targetScan_target = list(CLASH_targetScan_list)
        Nontarget = list(Nontarget_list)

        all_gene_pd = self.Deseq2_file
        all_gene_length = len(all_gene_pd.index)  ## all gene target matrix, sorted by values
        print(f'all gene length: {all_gene_length}')
        try:
            all_gene_pd['order'] = np.arange(0, 1, 1 / all_gene_length)
        except:
            all_gene_pd['order'] = np.arange(0, 1, 1 / all_gene_length)[:-1]

        all_target_pd = (self.Deseq2_file.loc[all_targets])  ## all target matrix, sort by value
        all_target_pd = all_target_pd.sort_values(by=['log2FoldChange'])
        target_length = len(all_target_pd.index)
        print(f'all target length: {target_length}')
        try:
            all_target_pd['order'] = np.arange(0, 1, 1 / target_length)
        except:
            all_target_pd['order'] = np.arange(0, 1, (1 / (target_length)))[:-1]

        conserved_target_pd = (self.Deseq2_file.loc[conserved_target])  ## conserved target matrix, sort by value
        conserved_target_pd = conserved_target_pd.sort_values(by=['log2FoldChange'])
        conserved_target_length = len(conserved_target_pd.index)
        print(f'conserved target length: {conserved_target_length}')
        try:
            conserved_target_pd['order'] = np.arange(0, 1, 1 / conserved_target_length)
        except:
            conserved_target_pd['order'] = np.arange(0, 1, 1 / conserved_target_length)[:-1]

        CLASH_target_pd = (self.Deseq2_file.loc[CLASH_list_target])  ## CLASH target matrix, sort by value
        CLASH_target_pd = CLASH_target_pd.sort_values(by=['log2FoldChange'])
        CLASH_target_length = len(CLASH_target_pd.index)
        print(f'conserved target length: {CLASH_target_length}')
        try:
            CLASH_target_pd['order'] = np.arange(0, 1, 1 / CLASH_target_length)
        except:
            CLASH_target_pd['order'] = np.arange(0, 1, 1 / CLASH_target_length)[:-1]

        CLASH_targetScan_target_pd = (
            self.Deseq2_file.loc[CLASH_targetScan_target])  ## CLASH target matrix, sort by value
        CLASH_targetScan_target_pd = CLASH_targetScan_target_pd.sort_values(by=['log2FoldChange'])
        CLASH_targetScan_target_length = len(CLASH_targetScan_target_pd.index)
        print(f'conserved target length: {CLASH_targetScan_target_length}')
        try:
            CLASH_targetScan_target_pd['order'] = np.arange(0, 1, 1 / CLASH_targetScan_target_length)
        except:
            CLASH_targetScan_target_pd['order'] = np.arange(0, 1, 1 / CLASH_targetScan_target_length)[:-1]

        Nontarget_pd = (self.Deseq2_file.loc[Nontarget])  ## Nontarget matrix, sort by value
        Nontarget_pd = (Nontarget_pd.sort_values(by=['log2FoldChange']))
        Nontarget_length = len(Nontarget_pd.index)
        print(f'Nontarget length: {Nontarget_length}')
        try:
            Nontarget_pd['order'] = np.arange(0, 1, 1 / Nontarget_length)
        except:
            Nontarget_pd['order'] = np.arange(0, 1, (1 / (Nontarget_length)))[:-1]
        return all_gene_pd, all_target_pd, conserved_target_pd, CLASH_target_pd, CLASH_targetScan_target_pd, Nontarget_pd

    def picture(self, TargetScan_all_targets='', TargetScan_conserved_targets=''):

        targets_list = self.output_target_Nontarget_list(
            miR_all_target_file=TargetScan_all_targets,
            conserved_target_file=TargetScan_conserved_targets)
        targets_de_table = self.table_target_Nontraget(
            all_target_list=targets_list[0],
            conserved_target_list=targets_list[1],
            Nontarget_list=targets_list[2])

        all_targets = targets_de_table[1]
        conserved_targets = targets_de_table[2]
        nontargets = targets_de_table[3]

        list_FC_all_targets = (list(all_targets['log2FoldChange']))
        list_FC_conserved_targets = (list(conserved_targets['log2FoldChange']))
        list_FC_nontargets = (list(nontargets['log2FoldChange']))
        U_all_targets, p_all_targets = mannwhitneyu(list_FC_all_targets, list_FC_nontargets, )
        U_conserved_targets, p_conserved_targets = mannwhitneyu(list_FC_conserved_targets, list_FC_nontargets, )
        print(f'all targets/ non targets: {p_all_targets}')
        print(f'conserved targets/ non targets: {p_conserved_targets}')
        fig, ax = plt.subplots()
        fig.set_size_inches(10, 10)
        sns.lineplot(data=all_targets, x='log2FoldChange', y='order', color='orange', linewidth=3,
                     linestyle='solid',
                     label=f'all predicted targets ({len(all_targets.index)})')
        sns.lineplot(data=conserved_targets, x='log2FoldChange', y='order', color='red',
                     linewidth=3,
                     linestyle='solid',
                     label=f'conserved predicted targets ({len(conserved_targets.index)})')
        sns.lineplot(data=nontargets, x='log2FoldChange', y='order', color='black', linewidth=3,
                     linestyle='solid',
                     label=f'Non target ({len(nontargets.index)})')
        ax.set_xlim(-1.5, 1.5)
        ax.set_ylim(0, 1.1)
        ax.set_xlabel('Fold change (log2)', fontsize=15)
        ax.set_ylabel('Cumulative fraction', fontsize=15)
        ax.set_title(f'{self.Deseq2_file_name}_baseMean_{self.baseMean_threshold}', fontsize=15)
        plt.savefig(f'{self.Deseq2_file_name}_baseMean_{self.baseMean_threshold}.svg')
        print('Done!')

    def picture_CLASH_targetScan(self, TargetScan_all_targets='', TargetScan_conserved_targets='', CLASH_targets='',
                                 CLASH_interacted_TargetScan_targets=''):

        targets_list = self.output_target_Nontarget_list_CLASH_targetScan(
            miR_all_target_file=TargetScan_all_targets,
            conserved_target_file=TargetScan_conserved_targets,
            CLASH_file=CLASH_targets,
            CLASH_targetScan_file=CLASH_interacted_TargetScan_targets)
        targets_de_table = self.table_target_Nontraget_CLASH(
            all_target_list=targets_list[0],
            conserved_target_list=targets_list[1],
            CLASH_list=targets_list[2],
            CLASH_targetScan_list=targets_list[3],
            Nontarget_list=targets_list[4])

        all_targets = targets_de_table[1]
        conserved_targets = targets_de_table[2]
        CLASH_targets = targets_de_table[3]
        CLASH_targetScan_targets = targets_de_table[4]
        nontargets = targets_de_table[5]

        list_FC_all_targets = (list(all_targets['log2FoldChange']))
        list_FC_conserved_targets = (list(conserved_targets['log2FoldChange']))
        list_FC_CLASH_targets = (list(CLASH_targets['log2FoldChange']))
        list_FC_CLASH_targetScan_targets = (list(CLASH_targetScan_targets['log2FoldChange']))
        list_FC_nontargets = (list(nontargets['log2FoldChange']))

        U_all_targets, p_all_targets = mannwhitneyu(list_FC_all_targets, list_FC_nontargets, )
        U_conserved_targets, p_conserved_targets = mannwhitneyu(list_FC_conserved_targets, list_FC_nontargets, )
        U_CLASH_targets, p_CLASH_targets = mannwhitneyu(list_FC_CLASH_targets, list_FC_nontargets, )
        U_CLASH_targetScan_targets, p_CLASH_targetScan_targets = mannwhitneyu(list_FC_CLASH_targetScan_targets,
                                                                              list_FC_nontargets, )
        print(f'all targets/ non targets: {U_all_targets}')
        print(f'conserved targets/ non targets: {U_conserved_targets}')
        print(f'\n')
        print(f'all targets/ non targets: {p_all_targets}')
        print(f'conserved targets/ non targets: {p_conserved_targets}')
        print(f'CLASH targets/ non targets: {p_CLASH_targets}')
        print(f'CLASH&targetScan targets/ non targets: {p_CLASH_targetScan_targets}')
        fig, ax = plt.subplots()
        fig.set_size_inches(10, 10)
        sns.lineplot(data=all_targets, x='log2FoldChange', y='order', color='orange', linewidth=3,
                     linestyle='solid',
                     label=f'all predicted targets ({len(all_targets.index)})')
        sns.lineplot(data=conserved_targets, x='log2FoldChange', y='order', color='red',
                     linewidth=3,
                     linestyle='solid',
                     label=f'conserved predicted targets ({len(conserved_targets.index)})')
        sns.lineplot(data=CLASH_targets, x='log2FoldChange', y='order', color='blue', linewidth=3,
                     linestyle='solid',
                     label=f'CLASH target ({len(CLASH_targets.index)})')
        sns.lineplot(data=CLASH_targetScan_targets, x='log2FoldChange', y='order', color='green', linewidth=3,
                     linestyle='solid',
                     label=f'CLASH_targetScan_targets ({len(CLASH_targetScan_targets.index)})')
        sns.lineplot(data=nontargets, x='log2FoldChange', y='order', color='black', linewidth=3,
                     linestyle='solid',
                     label=f'Non target ({len(nontargets.index)})')
        ax.set_xlim(-1.5, 1.5)
        ax.set_ylim(0, 1.1)
        ax.set_xlabel('Fold change (log2)', fontsize=15)
        ax.set_ylabel('Cumulative fraction', fontsize=15)
        ax.set_title(f'{self.Deseq2_file_name}_baseMean_{self.baseMean_threshold}', fontsize=15)
        plt.savefig(f'{self.Deseq2_file_name}_baseMean_{self.baseMean_threshold}.svg')
        print('Done!')


class MicroRNA_isoform:
    def __init__(self, ):
        pass

    def second_18th(self, miRNA_database, input_file1):
        print("miRNA length distribution criteria:\n1. allow two extension at 5' of miRNA\n2. 2nd to 18th of miRNA must match reads\n3.miRNAs with a 3' tail can extend up to a maximum length of 30 nucleotides.\n4.miRNA tail can be Templated and Non-templated tail\n")
        dict_miRNA_database = {}
        with open(miRNA_database, 'r+') as f1:
            for line1 in f1:
                if '>' in line1:
                    miRNA_name = line1.strip().split(' ')[0][1:]
                else:
                    miRNA_seq = line1.strip().upper()
                    miRNA_seq_2nd_18th = miRNA_seq[1:18]
                    dict_miRNA_database[miRNA_seq_2nd_18th] = dict_miRNA_database.get(miRNA_seq_2nd_18th, '')
                    dict_miRNA_database[miRNA_seq_2nd_18th] = (dict_miRNA_database[miRNA_seq_2nd_18th] + '_' + miRNA_name).strip('_')
        dict_miRNA_isoform_count = {}
        with open(input_file1, 'r+') as f1:
            dict_miRNA_isoform_count['count'] = dict_miRNA_isoform_count.get('count', {})
            for line1 in f1:
                if ('>' not in line1) and (len(line1.strip()) <= 30):
                    miRNA_isoform_read = line1.strip()
                    candidate_miRNA_read_0_17 = miRNA_isoform_read[0:17]
                    candidate_miRNA_read_1_18 = miRNA_isoform_read[1:18]
                    candidate_miRNA_read_2_19 = miRNA_isoform_read[2:19]
                    candidate_miRNA_read_3_20 = miRNA_isoform_read[3:20]
                    for each_candidate_miRNA_read in [candidate_miRNA_read_0_17, candidate_miRNA_read_1_18,
                                                      candidate_miRNA_read_2_19, candidate_miRNA_read_3_20]:
                        if each_candidate_miRNA_read in dict_miRNA_database:
                            each_miRNA_name1 = dict_miRNA_database[each_candidate_miRNA_read]
                            miRNA_name_isoformSeq_table = (f"{each_miRNA_name1}", f"{miRNA_isoform_read}")
                            dict_miRNA_isoform_count['count'][miRNA_name_isoformSeq_table] = dict_miRNA_isoform_count[
                                'count'].get(miRNA_name_isoformSeq_table, 0)
                            dict_miRNA_isoform_count['count'][miRNA_name_isoformSeq_table] += 1
        f1 = pd.DataFrame.from_dict(dict_miRNA_isoform_count)
        f1.to_csv(input_file1.replace('fasta', 'csv').replace('fastq', 'csv').replace('fq', 'csv').replace('fa', 'csv'))


class Sprite():
    def fasta_Barcode_identification(self, Barcode_file, input_fasta):
        print(
            'Criteria for barcode identidication\n1. 6 barcode must in the reads\n2. After removing Adapter sequence ATCAGCACTTAGCGTCAG, the rest read length must >=16nt\n')
        output_fasta = input_fasta.replace('fasta', 'barcoded.fasta')
        with open(Barcode_file, 'r+') as f1:
            tag1, tag2, tag3, tag4, tag5, tag6 = [], [], [], [], [], []
            dict_barcode = {}
            for line1 in f1:
                line2 = line1.strip().split('\t')
                name1 = line2[0]
                seq1 = line2[1]
                seq1_RC = str(Seq(seq1).reverse_complement())
                dict_barcode[seq1_RC] = name1
                if 'Tag1' in name1:
                    tag1.append(seq1_RC)
                if 'Tag2' in name1:
                    tag2.append(seq1_RC)
                if 'Tag3' in name1:
                    tag3.append(seq1_RC)
                if 'Tag4' in name1:
                    tag4.append(seq1_RC)
                if 'Tag5' in name1:
                    tag5.append(seq1_RC)
                if 'Tag6' in name1:
                    tag6.append(seq1_RC)
        tag1 = ('|'.join(tag1))
        tag2 = ('|'.join(tag2))
        tag3 = ('|'.join(tag3))
        tag4 = ('|'.join(tag4))
        tag5 = ('|'.join(tag5))
        tag6 = ('|'.join(tag6))
        f2 = open(output_fasta, 'w+')
        with open(input_fasta, 'r+') as f1:
            for line1 in f1:
                if '>' in line1:
                    name1 = line1.strip()
                else:
                    seq1 = line1.strip()
                    candidate_code = re.findall(
                        rf'({tag1})' + '.{5,10}' + f'({tag2})' + '.{5,10}' + f'({tag3})' + '.{5,10}' + f'({tag4})' + '.{5,10}' + f'({tag5})' + '.{5,10}' + f'({tag6})',
                        seq1)
                    if candidate_code != []:
                        barcode = '_'.join(candidate_code[0])
                        new_name1 = name1 + 'Barcode:' + barcode
                        new_seq1 = seq1.split('ATCAGCACTTAGCGTCAG')[0]
                        if len(new_seq1) >= 16:
                            f2.write(f"{new_name1}\n{new_seq1}\n")
        print('done!')
        f2.close()

    def sprite_SAM_distribution_RNA_type_statistic(self, SAM_file):
        print("Reminder:\n1.MIR**HG_lncRNA and VMP1_mRNA identified as microRNA")
        with open(SAM_file, 'r+') as f1:
            dict_reads = {}
            for line1 in f1:
                line2 = line1.strip().split('\t')
                if len(line2) > 10 and 'AS:i:' in line1:
                    raw_name = line2[0]
                    raw_transcript = line2[2]
                    read_name = raw_name.split('Barcode:')[0]
                    Barcode_name = raw_name.split('Barcode:')[1]
                    alignment_string = (re.search(r'AS:i:-*\d+', line1).group())
                    alignment_score = int(re.search(r'\d+', alignment_string).group())
                    aligned_pos = int(line2[3])
                    CIGAR = line2[5]
                    total_length = sum(list(int(re.search(r'\d+', x).group()) for x in (re.findall(r'\d+M', CIGAR))))
                    if 'microRNA' in raw_transcript:
                        mapped_transcript = raw_transcript.split('_')[2] + '_microRNA'
                    if 'mRNA' in raw_transcript:
                        RNA_type = raw_transcript.split('#')[-1].split('_')[0]
                        mapped_transcript = raw_transcript.split('_')[2].split('#')[0] + '_' + RNA_type
                        mapped_transcript = mapped_transcript.replace('VMP1_mRNA', 'miRNA_microRNA')
                        if re.search(r'MIR.{0,20}HG_lncRNA', mapped_transcript) != None:
                            mapped_transcript = mapped_transcript.replace(
                                re.search(r'MIR.{0,20}HG_lncRNA', mapped_transcript).group(), 'miRNA_microRNA')

                    dict_reads[Barcode_name] = dict_reads.get(Barcode_name, {})
                    dict_reads[Barcode_name][read_name] = dict_reads[Barcode_name].get(read_name, [mapped_transcript,
                                                                                                   alignment_score])
                    if ('microRNA' in mapped_transcript) and (
                            'microRNA' not in dict_reads[Barcode_name][read_name][0]):
                        dict_reads[Barcode_name][read_name] = [mapped_transcript, alignment_score]
                    if ('microRNA' in mapped_transcript) and (
                            'microRNA' in dict_reads[Barcode_name][read_name][0]):
                        if alignment_score > dict_reads[Barcode_name][read_name][1]:
                            dict_reads[Barcode_name][read_name] = [mapped_transcript, alignment_score]
                    if ('microRNA' not in mapped_transcript) and (
                            'microRNA' in dict_reads[Barcode_name][read_name][0]):
                        pass
                    if ('microRNA' not in mapped_transcript) and (
                            'microRNA' not in dict_reads[Barcode_name][read_name][0]):
                        if alignment_score > dict_reads[Barcode_name][read_name][1]:
                            dict_reads[Barcode_name][read_name] = [mapped_transcript, alignment_score]

        dict_barcode_distribution = {}
        dict_RNA_type = {}
        for barcode1, each_read_info in dict_reads.items():
            RNA_num = len(each_read_info)
            dict_barcode_distribution[RNA_num] = dict_barcode_distribution.get(RNA_num, 0)
            dict_barcode_distribution[RNA_num] += 1
            for read_ID in each_read_info:
                gene_name = each_read_info[read_ID][0]
                RNA_type = gene_name.split('_')[-1].split('&')[-1]
                dict_RNA_type[RNA_type] = dict_RNA_type.get(RNA_type, 0)
                dict_RNA_type[RNA_type] += 1
        sorted_dict_RNA_type = sorted(dict_RNA_type.items(), key=lambda x: x[1], reverse=True)
        sorted_dict_barcode_distribution = sorted(dict_barcode_distribution.items(), key=lambda x: x[0], reverse=False)
        print(sorted_dict_RNA_type)
        print(sorted_dict_barcode_distribution)

    def sprite_SAM_table(self, Transcript_database, SAM_file):
        print(
            f'Reminder:\n1.The reads sequence have been corrected based on transcript sequence\n2. each cluster only has two/three/four RNAs')
        with open(Transcript_database, 'r+') as f1:
            dict_database = {}
            name1, seq1 = '', ''
            for line1 in f1:
                if '>' in line1:
                    if name1 != '':
                        dict_database[name1] = seq1
                    name1 = line1.strip()[1:]
                    seq1 = ''
                else:
                    seq1 += line1.strip()
            dict_database[name1] = seq1

        with open(SAM_file, 'r+') as f1:
            dict_reads = {}
            for line1 in f1:
                line2 = line1.strip().split('\t')
                if len(line2) > 10 and 'AS:i:' in line1:
                    raw_name = line2[0]
                    raw_transcript = line2[2]
                    read_name = raw_name.split('Barcode:')[0]
                    Barcode_name = raw_name.split('Barcode:')[1]
                    alignment_string = (re.search(r'AS:i:-*\d+', line1).group())
                    alignment_score = int(re.search(r'\d+', alignment_string).group())
                    aligned_pos = int(line2[3])
                    CIGAR = line2[5]
                    total_length = sum(list(int(re.search(r'\d+', x).group()) for x in (re.findall(r'\d+M', CIGAR))))
                    if 'microRNA' in raw_transcript:
                        corrected_sequence = dict_database[raw_transcript]
                        mapped_transcript = raw_transcript.split('_')[2] + '_microRNA'
                    if 'mRNA' in raw_transcript:
                        corrected_sequence = dict_database[raw_transcript][aligned_pos - 1:aligned_pos + total_length]
                        RNA_type = raw_transcript.split('#')[-1].split('_')[0]
                        mapped_transcript = raw_transcript.split('_')[2].split('#')[0] + '_' + RNA_type
                    dict_reads[Barcode_name] = dict_reads.get(Barcode_name, {})
                    transcript_seq1 = mapped_transcript + '_' + corrected_sequence
                    dict_reads[Barcode_name][read_name] = dict_reads[Barcode_name].get(read_name, [transcript_seq1,
                                                                                                   alignment_score])
                    if read_name in dict_reads[Barcode_name]:
                        if ('microRNA' in mapped_transcript) and (
                                'microRNA' not in dict_reads[Barcode_name][read_name][0]):
                            dict_reads[Barcode_name][read_name] = [transcript_seq1, alignment_score]
                        elif ('microRNA' in mapped_transcript) and (
                                'microRNA' in dict_reads[Barcode_name][read_name][0]):
                            if alignment_score > dict_reads[Barcode_name][read_name][1]:
                                dict_reads[Barcode_name][read_name] = [transcript_seq1, alignment_score]
                        elif ('microRNA' not in mapped_transcript) and (
                                'microRNA' in dict_reads[Barcode_name][read_name][0]):
                            pass
                        elif ('microRNA' not in mapped_transcript) and (
                                'microRNA' not in dict_reads[Barcode_name][read_name][0]):
                            if alignment_score > dict_reads[Barcode_name][read_name][1]:
                                dict_reads[Barcode_name][read_name] = [transcript_seq1, alignment_score]
        dict_hyb2, dict_hyb3, dict_hyb4 = {}, {}, {}  ## These dictionaries show the number of RNAs in each cluster.
        for barcode, read_info in dict_reads.items():
            list_hyb = []
            for read_name, mapped_transcript_info in read_info.items():
                list_hyb += [mapped_transcript_info[0]]
            if len(list_hyb) == 2:
                name1 = ';'.join(sorted(list_hyb))
                dict_hyb2[name1] = dict_hyb2.get(name1, 0)
                dict_hyb2[name1] += 1
            elif len(list_hyb) == 3:
                name1 = ';'.join(sorted(list_hyb))
                dict_hyb3[name1] = dict_hyb3.get(name1, 0)
                dict_hyb3[name1] += 1
            elif len(list_hyb) == 4:
                name1 = ';'.join(sorted(list_hyb))
                dict_hyb4[name1] = dict_hyb4.get(name1, 0)
                dict_hyb4[name1] += 1
        f = pd.DataFrame.from_dict(dict_hyb2, orient='index')
        f.rename(columns={0: 'abundance'}, inplace=True)
        f.sort_values(by=['abundance'], inplace=True, ascending=False)
        f.to_csv(SAM_file.replace('sam', 'cluster2.csv'))

        f = pd.DataFrame.from_dict(dict_hyb3, orient='index')
        f.rename(columns={0: 'abundance'}, inplace=True)
        f.sort_values(by=['abundance'], inplace=True, ascending=False)
        f.to_csv(SAM_file.replace('sam', 'cluster3.csv'))

        f = pd.DataFrame.from_dict(dict_hyb4, orient='index')
        f.rename(columns={0: 'abundance'}, inplace=True)
        f.sort_values(by=['abundance'], inplace=True, ascending=False)
        f.to_csv(SAM_file.replace('sam', 'cluster4.csv'))
        print('Done!')

    def sprite_table_filter(self, table_csv):
        out_file = table_csv.replace('csv', 'filtered.csv')
        print(
            'Criteria:\n'
            '1. each cluster has 1 microRNA and 1 other RNA\n'
            '2. VMP1_mRNA .....TATCAGACTGATG.... as miR-21-5p\n'
            '3. VMP1_mRNA CAACACCAGTCGATGGGCTGT as miR-21-3p\n'
            '4. HOXA9_mRNA TAGGTAGTTTCCTGTTGTTGGG as miR-196b-3p\n'
            '5. remove tRNA and rRNA\n'
            '6. keep 1 microRNA and 1 other RNA\n'
            '7. The cluster should contain a minimum of one microRNA and one additional RNA.\n')

        f1 = pd.read_csv(table_csv, index_col=[0])
        dict_new_filter = {}  # new filtered cluster
        for index1 in f1.index:
            abundance1 = int(f1.loc[index1])
            list_need_to_delete = []
            index1_split = index1.split(';')
            for each_RNA_pos in range(len(index1_split)):
                if (re.search(r'VMP1_mRNA.{0,10}TATCAGACTGATG.{0,20}', index1_split[each_RNA_pos]) != None):
                    index1_split[each_RNA_pos] = 'miR-21-5p_microRNA_TAGCTTATCAGACTGATGTTGA'
                if (re.search(r'VMP1_mRNA.{0,10}CAACACCAGTCGATGGGCTGT.{0,20}', index1_split[each_RNA_pos]) != None):
                    index1_split[each_RNA_pos] = 'miR-21-3p_microRNA_CAACACCAGTCGATGGGCTGT'
                if (re.search(r'HOXA9_mRNA.{0,10}TAGGTAGTTTCCTGTTGTTGGG.{0,20}', index1_split[each_RNA_pos]) != None):
                    index1_split[each_RNA_pos] = 'miR-196b-5p_microRNA_TAGGTAGTTTCCTGTTGTTGGG'
                if 'rRNA' in index1_split[each_RNA_pos]:
                    list_need_to_delete.append(index1_split[each_RNA_pos])
                if 'tRNA' in index1_split[each_RNA_pos]:
                    list_need_to_delete.append(index1_split[each_RNA_pos])
            index1_split = sorted(set(x for x in index1_split if x not in list_need_to_delete))
            new_index1 = ';'.join(index1_split)
            if 'microRNA' not in new_index1:
                pass
            elif new_index1.count('microRNA') == new_index1.count(';') + 1:
                pass
            else:
                dict_new_filter[new_index1] = dict_new_filter.get(new_index1, 0)
                dict_new_filter[new_index1] += abundance1

        df = pd.DataFrame(list(dict_new_filter.items()), columns=['cluster', 'count'])
        df = df.sort_values('count', ascending=False)
        df.to_csv(out_file, index=False)

    def sprite_filtered_table_to_fasta(self, input_csv):
        f_mi = open(f'{input_csv}_miRNA_read.fa', 'w+')
        f_other_RNA = open(f'{input_csv}_other_RNA_read.fa', 'w+')
        f1 = pd.read_csv(input_csv, index_col=[0])

        for index1 in f1.index:
            index_split = index1.split(';')
            abundance1 = int(f1.loc[index1])
            microRNA_list = [x for x in index_split if 'microRNA' in x]
            other_RNA_list = [x for x in index_split if 'microRNA' not in x]
            for miRNA_info1 in microRNA_list:
                for other_RNA_info1 in other_RNA_list:
                    miRNA_name = '_'.join(miRNA_info1.split('_')[:-1]) + f'_abundance_{abundance1};'
                    miRNA_seq = miRNA_info1.split('_')[-1]
                    other_RNA_name = '_'.join(other_RNA_info1.split('_')[:-1])
                    other_RNA_seq = other_RNA_info1.split('_')[-1]
                    f_mi.write(f">{miRNA_name}\n{miRNA_seq}\n")
                    f_other_RNA.write(f">{other_RNA_name}\n{other_RNA_seq}\n")
        f_mi.close()
        f_other_RNA.close()

    def cluster_info1(self, list_nt1, list_pos1, list_basepattern1, dG1, cluster_name1):
        cluster_seq = ''.join(list_nt1)
        second_RNA_pos = [i for i, x in enumerate(list_pos1) if x == '1'][-1]
        string_basepattern = ''.join(['b' if int(element) > 0 else '.' for element in list_basepattern1])
        miRNA_name = re.search(r'.{0,30}microRNA|let.{0,20}microRNA', cluster_name1).group()
        cluster_abundance = int(re.search(r'abundance_\d+', cluster_name1).group().split('_')[-1])
        other_RNA_name = cluster_name1.split(';-')[1]
        miRNA_pos = cluster_name1.find(miRNA_name)
        dG1 = '-' + str(re.search(r'\d+.*\d*', dG1).group())
        if miRNA_pos == 0:  # miRNA on the left
            miRNA_seq = cluster_seq[:second_RNA_pos]
            miRNA_pattern = string_basepattern[:second_RNA_pos].replace('b', '(')
            other_RNA_seq = cluster_seq[second_RNA_pos:]
            other_RNA_pattern = string_basepattern[second_RNA_pos:].replace('b', ')')
            each_cluster_info = [miRNA_name, miRNA_seq, miRNA_pattern, other_RNA_name, other_RNA_seq, other_RNA_pattern,
                                 dG1, str(cluster_abundance)]
        elif miRNA_pos != 0:  # miRNA on the right
            print('Error, microRNA shouldn\'t on the left, check your input')
        return each_cluster_info

    def ct_format_to_basepattern_table(self, input_ct):
        f_out = open(input_ct.replace('.ct', '.ct.basepattern.txt'), 'w+')
        list_header = ['miRNA_name', 'miRNA_seq', 'miRNA_pattern', 'target_name', 'target_seq', 'target_pattern', 'dG',
                       'abundance']
        f_out.write('\t'.join(list_header) + '\n')
        with open(input_ct, 'r+') as f1:
            list_pos = []
            list_nt = []
            list_basepattern = []
            for line1 in f1:
                line2 = line1.strip().split('\t')
                if len(line2) == 4:
                    if list_nt != []:
                        each_cluster_info = self.cluster_info1(list_nt1=list_nt, list_pos1=list_pos,
                                                               list_basepattern1=list_basepattern, dG1=dG,
                                                               cluster_name1=cluster_name)
                        f_out.write('\t'.join(each_cluster_info) + '\n')
                    dG = line2[1]
                    cluster_name = line2[3]
                    list_nt = []
                    list_pos = []
                    list_basepattern = []
                elif len(line2) == 8:
                    list_nt.append(line2[1])
                    list_pos.append(line2[5])
                    list_basepattern.append(line2[4])
            each_cluster_info = self.cluster_info1(list_nt1=list_nt, list_pos1=list_pos,
                                                   list_basepattern1=list_basepattern, dG1=dG,
                                                   cluster_name1=cluster_name)
            f_out.write('\t'.join(each_cluster_info) + '\n')
        f_out.close()

    def sprite_final_table_statistic(self, input_csv):
        dG_16 = 0
        dG_12TO16 = 0
        dG_12 = 0
        seed_match_num = 0
        f1 = pd.read_table(input_csv)
        for indexNUM in f1.index:
            microRNA_seed = f1.iloc[indexNUM]['miRNA_seq'][1:8]
            microRNA_seed_RC = str(Seq(microRNA_seed).reverse_complement())
            target_seq = f1.iloc[indexNUM]['target_seq']
            dG_num = float(f1.iloc[indexNUM]['dG'])
            seed_match_num += microRNA_seed_RC in target_seq
            dG_16 += (dG_num <= -16)
            dG_12TO16 += (-16 < dG_num <= -12)
            dG_12 += (dG_num > -12)
        print(
            f"In the {input_csv}\nTotal cluster number is: \t{len(f1)}\nSeed match number is:    \t{seed_match_num}\ndG less than -16 is:    \t{dG_16}\ndG between -12 to -16 is: \t{dG_12TO16}\ndG larger than -12 is:    \t{dG_12}")

class AQ_seq():
    def cutadapt(self,fastq1,fastq2,fastq_out):
        print('The minimum size for adapter trimming in AQ-seq cutadapt is 24 base pairs (bp).')
        fastq1_out = f"{fastq_out}_1.fastq"
        fastq2_out = f"{fastq_out}_2.fastq"
        command = f'module load cutadapt && cutadapt -a TGGAATTCTCGGGTGCCAAG -A GATCGTCGGACTGTAGAACT -o {fastq1_out} -p {fastq2_out} {fastq1} {fastq2} --minimum-length 24 -j 10'
        # Run the command and capture the output
        completed_process = subprocess.run(command, shell=True, capture_output=True, text=True)
        # Get the output and error messages
        output = completed_process.stdout
        error = completed_process.stderr
        # Save the output to a file
        with open("cutadapt_out.txt", "a+") as output_file:
            output_file.write(output)
        with open("cutadapt_out_error.txt", "a+") as output_file:
            output_file.write(error)

    def pear(self,fastq_in):
        print('The minimum size for peared reads is 24 base pairs (bp).')
        fastq1_in = f"{fastq_in}_1.fastq"
        fastq2_in = f"{fastq_in}_2.fastq"
        command = f'module load pear && pear -f {fastq1_in} -r {fastq2_in} -o {fastq_in} -n 24 -j 10'
        # Run the command and capture the output
        completed_process = subprocess.run(command, shell=True, capture_output=True, text=True)
        # Get the output and error messages
        output = completed_process.stdout
        error = completed_process.stderr
        # Save the output to a file
        with open("pear_out.txt", "a+") as output_file:
            output_file.write(output)
        with open("pear_error.txt", "a+") as output_file:
            output_file.write(error)

    def collapse(self,fastq_in):
        print('The reads are currently being collapsed using a collapser.')
        command = f'module load fastx_toolkit && fastx_collapser -i {fastq_in}.assembled.fastq > {fastq_in}.collapsed.fasta'
        # Run the command and capture the output
        completed_process = subprocess.run(command, shell=True, capture_output=True, text=True)
        # Get the output and error messages
        output = completed_process.stdout
        error = completed_process.stderr
        # Save the output to a file
        with open("collapser_out.txt", "a+") as output_file:
            output_file.write(output)
        with open("collapser_error.txt", "a+") as output_file:
            output_file.write(error)

    def cut_UMI(self,fastq_in):
        print('The minimum size for UMI trimming in AQ-seq cutadapt is 16 base pairs (bp).')
        fastq_collapsed = f"{fastq_in}.collapsed.fasta"
        Cut_UMI_fastq = f"{fastq_in}.cut8N.fasta"
        command = f'module load cutadapt && cutadapt -u 4 -u -4  {fastq_collapsed} -o {Cut_UMI_fastq} --minimum-length 16 -j 10'
        # Run the command and capture the output
        completed_process = subprocess.run(command, shell=True, capture_output=True, text=True)
        # Get the output and error messages
        output = completed_process.stdout
        error = completed_process.stderr
        # Save the output to a file
        with open("cutadapt_UMI_out.txt", "a+") as output_file:
            output_file.write(output)
        with open("cutadapt_UMI_out_error.txt", "a+") as output_file:
            output_file.write(error)
    def isformCSV_totalAbundanceCSV(self,input1):
        df = pd.read_csv(f'{input1}', header=0)
        grouped_df = df.groupby(df.columns[0])[df.columns[2]].sum().reset_index()
        grouped_df.to_csv(f'{input1}.abundance.csv', index=False, header=True)

if __name__ == "__main__":
    # print('This script is edited by Lu Li from Mingyi Xie Lab, University of Florida')

    def command_line():
        argv_step_command = sys.argv[1]
        argv = sys.argv[2:]
        try:
            if argv_step_command == 'AQ-seq':
                opts, args = getopt.getopt(argv, '1:2:o:d:',['fastq_1=', 'fastq_2=', 'fastq_out=','miR_database='])
                for opt in opts:
                    if opt[0] == '-1':
                        input_file1 = opt[1]
                    elif opt[0] == '-2':
                        input_file2 = opt[1]
                    elif opt[0] == '-o':
                        out1 = opt[1]
                    elif opt[0] == '-d':
                        database1 = opt[1]
                print(opts)
                AQ_seq().cutadapt(fastq1= input_file1,fastq2=input_file2,fastq_out=out1)
                AQ_seq().pear(fastq_in=out1)
                AQ_seq().collapse(fastq_in=out1)
                AQ_seq().cut_UMI(fastq_in=out1)
                MicroRNA_isoform().second_18th(input_file1=f"{out1}", miRNA_database=database1)
                AQ_seq().isformCSV_totalAbundanceCSV(input1=out1)
                print(f'Done!\n')
            elif argv_step_command == 'making_unique_redundant_database_hg38':
                opts, args = getopt.getopt(argv, 'i:', ['input='])  # getopt.getopt(args, shortopts, longopts=[])
                input_martquery_database = opts[0][1]
                print(input_martquery_database)
                Database().making_unique_redundant_database_hg38(input=input_martquery_database)  # step 1
            elif argv_step_command == 'making_unique_redundant_database_mm39':
                opts, args = getopt.getopt(argv, 'i:', ['input='])  # getopt.getopt(args, shortopts, longopts=[])
                input_martquery_database = opts[0][1]
                print(input_martquery_database)
                Database().making_unique_redundant_database_mm39(input=input_martquery_database)  # step

            elif argv_step_command == 'making_unique_redundant_database_WBcel235':
                opts, args = getopt.getopt(argv, 'i:', ['input='])  # getopt.getopt(args, shortopts, longopts=[])
                input_martquery_database = opts[0][1]
                print(input_martquery_database)
                Database().making_unique_redundant_database_WBcel235(input=input_martquery_database)  # step 1

            elif argv_step_command == 'element_region_statistic':
                opts, args = getopt.getopt(argv, 'i:', ['input=', ])  # getopt.getopt(args, shortopts, longopts=[])
                for opt in opts:
                    if opt[0] == '-i':
                        input_file1 = opt[1]
                print(opts)
                Statistic().Viennad_table_element_region_statistic(input_table=input_file1)

            elif argv_step_command == 'miRNA_abundance_in_hybrids':
                opts, args = getopt.getopt(argv, 'i:', ['input=', ])  # getopt.getopt(args, shortopts, longopts=[])
                for opt in opts:
                    if opt[0] == '-i':
                        input_file1 = opt[1]
                print(opts)
                Statistic().miRNA_abundance_in_hybrids(hyb_file=input_file1)

            elif argv_step_command == 'dG':
                opts, args = getopt.getopt(argv, 'i:', ['input=', ])  # getopt.getopt(args, shortopts, longopts=[])
                for opt in opts:
                    if opt[0] == '-i':
                        input_file1 = opt[1]
                print(opts)
                Statistic().dG_energy_bar_chart(input1=input_file1)

            elif argv_step_command == 'SAM_bowtie2_count':
                opts, args = getopt.getopt(argv, 'i:', ['input=', ])  # getopt.getopt(args, shortopts, longopts=[])
                for opt in opts:
                    if opt[0] == '-i':
                        input_file1 = opt[1]
                print(opts)
                SAM().bowtie2_count(input1=input_file1)

            elif argv_step_command == 'SAM_bowtie2_hOH7_single_Transcripts_count':
                opts, args = getopt.getopt(argv, 'i:', ['input=', ])  # getopt.getopt(args, shortopts, longopts=[])
                for opt in opts:
                    if opt[0] == '-i':
                        input_file1 = opt[1]
                print(opts)
                SAM().bowtie2_hOH7_single_Transcripts_count(input1=input_file1)

            elif argv_step_command == 'SAM_hisat2_mRNA_count':
                opts, args = getopt.getopt(argv, 'i:', ['input=', ])  # getopt.getopt(args, shortopts, longopts=[])
                for opt in opts:
                    if opt[0] == '-i':
                        input_file1 = opt[1]
                print(opts)
                SAM().hisat2_SAM_mRNA_count(input1=input_file1)

            elif argv_step_command == 'Hisat2_SAM_gene_Table':
                opts, args = getopt.getopt(argv, 'i:', ['input=', ])  # getopt.getopt(args, shortopts, longopts=[])
                for opt in opts:
                    if opt[0] == '-i':
                        input_file1 = opt[1]
                print(opts)
                SAM().Hisat2_SAM_gene_Table(input1=input_file1)

            elif argv_step_command == 'FASTQ_length_distribution':
                opts, args = getopt.getopt(argv, 'i:', ['input=', ])  # getopt.getopt(args, shortopts, longopts=[])
                for opt in opts:
                    if opt[0] == '-i':
                        input_file1 = opt[1]
                        print(input_file1)
                print(opts)
                Statistic().fastq_length_distribution(input1=input_file1)
            elif argv_step_command == 'FASTA_reads_number':
                opts, args = getopt.getopt(argv, 'i:m:M:',
                                           ['input=', 'minmium_length=', 'maxmium_length='])
                for opt in opts:
                    if opt[0] == '-i':
                        input_file1 = opt[1]
                    elif opt[0] == '-m':
                        min1 = opt[1]
                    elif opt[0] == '-M':
                        max1 = opt[1]
                print(opts)
                Statistic().fasta_reads_number(input1=input_file1, min_len=min1, max_len=max1)

            elif argv_step_command == 'FASTA_length_distribution':
                opts, args = getopt.getopt(argv, 'i:', ['input=', ])  # getopt.getopt(args, shortopts, longopts=[])
                for opt in opts:
                    if opt[0] == '-i':
                        input_file1 = opt[1]
                        print(input_file1)
                print(opts)
                Statistic().fasta_length_distribution(input1=input_file1)

            elif argv_step_command == 'Viennad_Rnatype_number':
                opts, args = getopt.getopt(argv, 'i:d:', ['input=','microRNA_database'])  # getopt.getopt(args, shortopts, longopts=[])
                for opt in opts:
                    if opt[0] == '-i':
                        input_file1 = opt[1]
                    elif opt[0] == '-d':
                        database1 = opt[1]
                print(opts)
                Statistic().viennad_rnatype_number(viennad_file=input_file1, microRNA_database=database1)

            elif argv_step_command == 'miRNA_count_from_Viennad':
                opts, args = getopt.getopt(argv, 'i:d:', ['input=','microRNA_database'])  # getopt.getopt(args, shortopts, longopts=[])
                for opt in opts:
                    if opt[0] == '-i':
                        input_file1 = opt[1]
                    elif opt[0] == '-d':
                        database1 = opt[1]
                print(opts)
                Statistic().miRNA_count_from_Viennad(viennad_file=input_file1, microRNA_database=database1)

            elif argv_step_command == 'miRNA_count_from_whole_nt':
                opts, args = getopt.getopt(argv, 'i:d:', ['input=','microRNA_database'])  # getopt.getopt(args, shortopts, longopts=[])
                for opt in opts:
                    if opt[0] == '-i':
                        input_file1 = opt[1]
                    elif opt[0] == '-d':
                        database1 = opt[1]
                print(opts)
                Statistic().miRNA_count_from_whole_nt(input1=input_file1, microRNA_database=database1)

            elif argv_step_command == 'hybrid_number_statistic':
                opts, args = getopt.getopt(argv, 'i:d:', ['input=',
                                                          'microRNA_database'])  # getopt.getopt(args, shortopts, longopts=[])
                for opt in opts:
                    print(opt)
                    if opt[0] == '-i':
                        input_file1 = opt[1]
                    elif opt[0] == '-d':
                        database1 = opt[1]
                print(opts)
                Statistic().hyb_number(input1=input_file1, microRNA_database=database1)
            elif argv_step_command == 'reverse_complement_FASTA':
                opts, args = getopt.getopt(argv, 'i:', ['input=', ])  # getopt.getopt(args, shortopts, longopts=[])
                for opt in opts:
                    if opt[0] == '-i':
                        input_file1 = opt[1]
                FASTA().reverse_complementary_sequence(input=input_file1)

            elif argv_step_command == 'concatenated_reads':
                opts, args = getopt.getopt(argv, 'i:l:', ['input=',
                                                          'fragment_length='])  # getopt.getopt(args, shortopts, longopts=[])
                for opt in opts:
                    if opt[0] == '-i':
                        input_file1 = opt[1]
                    elif opt[0] == '-l':
                        fragment_len = opt[1]
                FASTA().concatenated_reads(input1=input_file1, fragment_length1=int(fragment_len))

            elif argv_step_command == 'deduplicate_BGI':
                opts, args = getopt.getopt(argv, 'i:', ['input='])  # getopt.getopt(args, shortopts, longopts=[])
                for opt in opts:
                    if opt[0] == '-i':
                        input_file1 = opt[1]
                print(opts)
                FASTQ().deduplicate(input_file=input_file1)
            elif argv_step_command == 'miRNA_CSV_merge':
                opts, args = getopt.getopt(argv, 's:', ['string='])  # getopt.getopt(args, shortopts, longopts=[])
                for opt in opts:
                    if opt[0] == '-s':
                        string_info = opt[1]
                Statistic().miRNA_CSV_merge(string1=string_info)

            elif argv_step_command == 'table_merge':
                opts, args = getopt.getopt(argv, 's:', ['string='])  # getopt.getopt(args, shortopts, longopts=[])
                for opt in opts:
                    if opt[0] == '-s':
                        string_info = opt[1]
                Statistic().table_merge(string1=string_info)

            elif argv_step_command == 'miRNA_CSV_two_index_merge':
                opts, args = getopt.getopt(argv, 's:', ['string='])  # getopt.getopt(args, shortopts, longopts=[])
                for opt in opts:
                    if opt[0] == '-s':
                        string_info = opt[1]
                Statistic().miRNA_CSV_two_index_merge(string1=string_info)

            elif argv_step_command == 'gene_table_two_index_merge':
                opts, args = getopt.getopt(argv, 's:', ['string='])  # getopt.getopt(args, shortopts, longopts=[])
                for opt in opts:
                    if opt[0] == '-s':
                        string_info = opt[1]
                Statistic().gene_table_two_index_merge(string1=string_info)

            elif argv_step_command == 'all_miRNA_isoform_number_table_1st_18th':
                opts, args = getopt.getopt(argv, 'i:d:l:', ['input=', 'miRNA_database=',
                                                            'longest_isoform_length'])  # getopt.getopt(args, shortopts, longopts=[])
                for opt in opts:
                    if opt[0] == '-i':
                        input_file1 = opt[1]
                    elif opt[0] == '-d':
                        database1 = opt[1]
                    elif opt[0] == '-l':
                        longest_isoform_length1 = opt[1]
                print(opts)
                print(
                    f"miRNA abundance criteria:\n1. The top 18 nts of miRNA in the beginning of reads\n2. The reads length are between 18 nts to {longest_isoform_length1} nts\n")
                Statistic().all_miRNA_isoform_number_table_1st_18th(input_file=input_file1, microRNA_database=database1,
                                                                    longest_isoform_length=longest_isoform_length1)

            elif argv_step_command == 'all_miRNA_isoform_table_2nd_18th':
                opts, args = getopt.getopt(argv, 'i:d:', ['input=',
                                                          'miRNA_database=', ])  # getopt.getopt(args, shortopts, longopts=[])
                for opt in opts:
                    if opt[0] == '-i':
                        input_file1 = opt[1]
                    elif opt[0] == '-d':
                        database1 = opt[1]
                print(opts)
                MicroRNA_isoform().second_18th(input_file1=input_file1, miRNA_database=database1)
                AQ_seq().isformCSV_totalAbundanceCSV(input1=input_file1.replace('.fasta','.csv').replace('.fastq','.csv').replace('.fa','.csv'))
                print(f'Done!\n')

            elif argv_step_command == 'all_miRNA_isoform_kind_count_1st_18th':
                opts, args = getopt.getopt(argv, 'i:d:l:', ['input=',
                                                            'miRNA_database=',
                                                            'longest_isoform_length'])  # getopt.getopt(args, shortopts, longopts=[])
                for opt in opts:
                    if opt[0] == '-i':
                        input_file1 = opt[1]
                    elif opt[0] == '-d':
                        database1 = opt[1]
                    elif opt[0] == '-l':
                        longest_isoform_length1 = opt[1]
                print(opts)
                print(
                    f"miRNA length distribution criteria:\n1. The top 18 nts of miRNA in the beginning of reads\n2. The reads length are between 18 nts to {longest_isoform_length1} nts\n")
                Statistic().all_miRNA_isoform_kind_count_1st_18th(input_file=input_file1, microRNA_database=database1,
                                                                  longest_isoform_length=longest_isoform_length1)

            elif argv_step_command == 'all_miRNA_length_distribution_1st_18th':
                opts, args = getopt.getopt(argv, 'i:d:l:', ['input=',
                                                            'miRNA_database=',
                                                            'longest_isoform_length'])  # getopt.getopt(args, shortopts, longopts=[])
                for opt in opts:
                    if opt[0] == '-i':
                        input_file1 = opt[1]
                    elif opt[0] == '-d':
                        database1 = opt[1]
                    elif opt[0] == '-l':
                        longest_isoform_length1 = opt[1]
                print(opts)
                print(
                    f"miRNA length distribution criteria:\n1. The top 18 nts of miRNA in the beginning of reads\n2. The reads length are between 18 nts to {longest_isoform_length1} nts\n")
                Statistic().all_miRNA_length_distribution_1st_18th(input_file=input_file1, microRNA_database=database1,
                                                                   longest_isoform_length=longest_isoform_length1)
            elif argv_step_command == 'miRNA_length_distribution':
                opts, args = getopt.getopt(argv, 'i:d:m:p:', ['input=', 'miRNA_database=',
                                                              'miRNA_sequence',
                                                              'primiRNA_sequence'])  # getopt.getopt(args, shortopts, longopts=[])
                for opt in opts:
                    if opt[0] == '-i':
                        input_file1 = opt[1]
                    elif opt[0] == '-d':
                        database1 = opt[1]
                    elif opt[0] == '-m':
                        miRNA_sequence1 = opt[1]
                    elif opt[0] == '-p':
                        primiRNA_sequence1 = opt[1]
                print(opts)
                print(
                    "miRNA length distribution criteria:\n1. The top 16 nts of miRNA in the beginning of reads\n2. The reads length are between 16 nts to 30 nts\n")
                Statistic().miRNA_length_distribution(input_file=input_file1, microRNA_database=database1,
                                                      miRNA_sequence=miRNA_sequence1,
                                                      primiRNA_sequence=primiRNA_sequence1)

            elif argv_step_command == 'miRNA_length_distribution_merge':
                Statistic().miRNA_length_distribution_merge()

            elif argv_step_command == 'hyb_length_distribution':  # show stacked barplot of hybrid length and abundance
                opts, args = getopt.getopt(argv, 'i:', ['input=', ])  # getopt.getopt(args, shortopts, longopts=[])
                for opt in opts:
                    if opt[0] == '-i':
                        input_file1 = opt[1]
                print(opts)
                Statistic().hyb_length(hyb_file=input_file1)
            elif argv_step_command == 'one_kind_hybrid_distribution':  # show one specific miRNA-target length and abundance
                opts, args = getopt.getopt(argv, 'i:m:t:', ['input=', 'miRNA_sequence=',
                                                            'target_sequence='])  # getopt.getopt(args, shortopts, longopts=[])
                for opt in opts:
                    if opt[0] == '-i':
                        input_file1 = opt[1]
                    elif opt[0] == '-m':
                        miRNA_sequence1 = opt[1]
                    elif opt[0] == '-t':
                        target_sequence1 = opt[1]
                print(opts)
                Statistic().single_miRNA_target_average_length(hyb_file=input_file1, miRNA_seq=miRNA_sequence1,
                                                               target_seq=target_sequence1)
            elif argv_step_command == 'RNA_RNA_hybrid_distribution':  # show all RNA-RNA hybrids and abundance
                opts, args = getopt.getopt(argv, 'i:', ['input=', ])  # getopt.getopt(args, shortopts, longopts=[])
                for opt in opts:
                    if opt[0] == '-i':
                        input_file1 = opt[1]
                print(opts)
                Statistic().RNA_RNA_hybrid_distribution(hyb_file=input_file1)
            elif argv_step_command == 'making_simple_BedGraph':
                opts, args = getopt.getopt(argv, 'i:c:',
                                           ['input=', 'chromosome='])  # getopt.getopt(args, shortopts, longopts=[])
                for opt in opts:
                    if opt[0] == '-i':
                        input_bedGraph = opt[1]
                    elif opt[0] == '-c':
                        chr_name = opt[1]
                print(input_bedGraph)
                print(opts)
                BedGraph().making_simple_version(file=input_bedGraph, chromosome_name=chr_name)  # step 3
            elif argv_step_command == 'making_transcript_sequence_genomeposition_conservation_database':
                opts, args = getopt.getopt(argv, 'g:t:c:', ['genome_database=', 'transcript_database=',
                                                            'chromosome='])  # getopt.getopt(args, shortopts, longopts=[])
                for opt in opts:
                    if opt[0] == '-g':
                        genome_database = opt[1]
                    elif opt[0] == '-t':
                        transcript_database = opt[1]
                    elif opt[0] == '-c':
                        chromosome = opt[1]
                print(opts)
                Database().making_transcript_sequence_genomeposition_conservation_database(
                    genome_file=genome_database, transcript_file=transcript_database, chr_name=chromosome)  # step 4
            elif argv_step_command == 'Viennad_to_Table':
                opts, args = getopt.getopt(argv, 't:c:i:n:',
                                           ['transcirpt_database=', 'transcript_ConservationScore_database=',
                                            'viennad=', 'name_database'])  # getopt.getopt(args, shortopts, longopts=[])
                for opt in opts:
                    if opt[0] == '-c':
                        transcript_CS = opt[1]
                    elif opt[0] == '-t':
                        transcript_DB = opt[1]
                    elif opt[0] == '-i':
                        input1 = opt[1]
                    if opt[0] == '-n':
                        Name_database1 = opt[1]
                print(opts)
                print("step I: convert Viennad to table and calculate conservation score, ~5 minutes")
                Viennad_to_table(transcript_ConservationScore_database=transcript_CS,
                                 transcript_only_database=transcript_DB).input_viennad(file=input1)
                print("\nstep II: Compressed table, ~2 minutes")
                Compressed_table().files(input=input1)
                print("\nstep III: Calculate miRNA AU extension in hybrids, ~2 minutes")
                MiRNA_CLASH_extension(infile=input1).miRNA_AU_ex_ratio(hyb_viennad_names=input1)
                print("\nStep IV: Calculate miRNA binding site position, ~2 minutes")
                Gene_region(ensemble_database=Name_database1).table(input=input1)
                print("\nStep V: TDMD analyzer, ~1 minutes")
                Combined_table().TDMD_candidate(input_file=input1 + '_compressed_AUex_region.txt')
                print('\ndone!')

            elif argv_step_command == 'Combined_table_normalize_value':
                opts, args = getopt.getopt(argv, 'i:r:o:', ['input_tables=', 'replicate_number=',
                                                            'output_name='])  # getopt.getopt(args, shortopts, longopts=[])
                for opt in opts:
                    if opt[0] == '-i':
                        viennad_tables = opt[1]
                    elif opt[0] == '-o':
                        output_name = opt[1]
                    elif opt[0] == '-r':
                        replicate_num = opt[1]
                print(opts)
                Combined_table(outfile=output_name, replicates=replicate_num).input_multiple_clash(
                    viennad_tables)
            elif argv_step_command == 'Combined_table_raw_value':
                opts, args = getopt.getopt(argv, 'i:o:', ['input_tables=', 'output_name='])
                for opt in opts:
                    if opt[0] == '-i':
                        file_names = opt[1]
                    elif opt[0] == '-o':
                        output_name = opt[1]
                print(opts)
                Combined_table(outfile=output_name, ).input_multiple_raw_clash(
                    viennad_files=file_names)
            elif argv_step_command == 'TDMD_analyzer':
                opts, args = getopt.getopt(argv, 'i:',
                                           ['input_table='])  # getopt.getopt(args, shortopts, longopts=[])
                for opt in opts:
                    if opt[0] == '-i':
                        input1 = opt[1]
                print(opts)
                Combined_table().TDMD_candidate(input_file=input1)
            elif argv_step_command == 'Cumulative_fraction_curve_targetScan':
                opts, args = getopt.getopt(argv, 'i:a:c:b:',
                                           ['input_tables=', 'all_targets=', 'conserved_targets', 'baseMean'])
                for opt in opts:
                    if opt[0] == '-i':
                        deseq_CSV = opt[1]
                    elif opt[0] == '-a':
                        all_file = opt[1]
                    elif opt[0] == '-c':
                        Conserved_file = opt[1]
                    elif opt[0] == '-b':
                        baseMean1 = opt[1]
                print(opts)
                Target_analysis(Deseq2_file=deseq_CSV, baseMean_threshold=baseMean1).picture(
                    TargetScan_all_targets=all_file,
                    TargetScan_conserved_targets=Conserved_file)
            elif argv_step_command == 'Cumulative_fraction_curve_targetScan_CLASH':
                opts, args = getopt.getopt(argv, 'i:a:c:b:',
                                           ['input_tables=', 'all_targets=', 'conserved_targets', 'baseMean'])
                for opt in opts:
                    if opt[0] == '-i':
                        deseq_CSV = opt[1]
                    elif opt[0] == '-a':
                        all_file = opt[1]
                    elif opt[0] == '-c':
                        Conserved_file = opt[1]
                    elif opt[0] == '-l':
                        CLASH_file = opt[1]
                    elif opt[0] == '-t':
                        CLASH_interacted_targetScan_file = opt[1]
                    elif opt[0] == '-b':
                        baseMean1 = opt[1]
                print(opts)
                Target_analysis(Deseq2_file=deseq_CSV, baseMean_threshold=baseMean1).picture_CLASH_targetScan(
                    TargetScan_all_targets=all_file, TargetScan_conserved_targets=Conserved_file,
                    CLASH_targets=CLASH_file, CLASH_interacted_TargetScan_targets=CLASH_interacted_targetScan_file)
            elif argv_step_command == 'insert_miRNA_target':
                opts, args = getopt.getopt(argv, 'i:m:t:', ['input=', 'miRNA_sequence=', 'target_sequence'])
                for opt in opts:
                    if opt[0] == '-i':
                        input_file = opt[1]
                    elif opt[0] == '-m':
                        miRNA_sequence = opt[1]
                    elif opt[0] == '-t':
                        target_sequence = opt[1]
                print(opts)
                Statistic().insert_between_miRNA_target(input1=input_file, str1=miRNA_sequence, str2=target_sequence)

            elif argv_step_command == 'fasta_Barcode_identification':
                opts, args = getopt.getopt(argv, 'b:i:', ['barcode=', 'input_fasta=', ])
                for opt in opts:
                    if opt[0] == '-b':
                        barcode1 = opt[1]
                    elif opt[0] == '-i':
                        fasta1 = opt[1]
                print(opts)
                Sprite().fasta_Barcode_identification(Barcode_file=barcode1, input_fasta=fasta1)
            elif argv_step_command == 'Sprite_SAM_to_Table':
                opts, args = getopt.getopt(argv, 't:i:', ['Transcript_database=', 'SAM_file=', ])
                for opt in opts:
                    if opt[0] == '-t':
                        t1 = opt[1]
                    elif opt[0] == '-i':
                        s1 = opt[1]
                print(opts)
                Sprite().sprite_SAM_table(Transcript_database=t1, SAM_file=s1)
            elif argv_step_command == 'Sprite_SAM_distribution_RNA_type_statistic':
                opts, args = getopt.getopt(argv, 'i:', ['SAM_file=', ])
                for opt in opts:
                    if opt[0] == '-i':
                        s1 = opt[1]
                print(opts)
                Sprite().sprite_SAM_distribution_RNA_type_statistic(SAM_file=s1)
            elif argv_step_command == 'Sprite_Table_filter':
                opts, args = getopt.getopt(argv, 'i:', ['table_csv=', ])
                for opt in opts:
                    if opt[0] == '-i':
                        s1 = opt[1]
                print(opts)
                Sprite().sprite_table_filter(table_csv=s1)

            elif argv_step_command == 'Sprite_filtered_table_to_fasta':
                opts, args = getopt.getopt(argv, 'i:', ['table_csv=', ])
                for opt in opts:
                    if opt[0] == '-i':
                        s1 = opt[1]
                print(opts)
                Sprite().sprite_filtered_table_to_fasta(input_csv=s1)

            elif argv_step_command == 'Sprite_ct_format_to_basepattern_table':
                opts, args = getopt.getopt(argv, 'i:', ['ct_file=', ])
                for opt in opts:
                    if opt[0] == '-i':
                        s1 = opt[1]
                print(opts)
                Sprite().ct_format_to_basepattern_table(input_ct=s1)

            elif argv_step_command == 'Sprite_final_table_statistic':
                opts, args = getopt.getopt(argv, 'i:', ['table_csv=', ])
                for opt in opts:
                    if opt[0] == '-i':
                        s1 = opt[1]
                print(opts)
                Sprite().sprite_final_table_statistic(input_csv=s1)

            elif argv_step_command == 'fastq_each_read_count':
                opts, args = getopt.getopt(argv, 'i:', ['input=', ])
                for opt in opts:
                    if opt[0] == '-i':
                        s1 = opt[1]
                print(opts)
                FASTQ().fastq_each_read_count(input1=s1)

            elif argv_step_command == 'miR16_sfGFP_count':
                opts, args = getopt.getopt(argv, 'i:', ['input=', ])
                for opt in opts:
                    if opt[0] == '-i':
                        s1 = opt[1]
                print(opts)
                Statistic().miR16_sfGFP_count(input1=s1)

            elif argv_step_command == 'Bowtie2_SAM_singleTranscript_shortReads_statistic':
                opts, args = getopt.getopt(argv, 'i:', ['input=', ])
                for opt in opts:
                    if opt[0] == '-i':
                        s1 = opt[1]
                print(opts)
                Statistic().Bowtie2_SAM_singleTranscript_shortReads_statistic(input1=s1)

            elif argv_step_command == 'Bowtie2_SAM_singleTranscript_longReads_statistic':
                opts, args = getopt.getopt(argv, 'i:', ['input=', ])
                for opt in opts:
                    if opt[0] == '-i':
                        s1 = opt[1]
                print(opts)
                Statistic().Bowtie2_SAM_singleTranscript_longReads_statistic(input1=s1)

            else:
                print('USAGE')
                print('-h [--help]')
                print('python3 CLASH.py insert_miRNA_target -i [--input] -m [--miRNA_sequence] -t [--target_sequence]')
                print('python3 CLASH.py dG -i [--input] <Viennad_file>')
                print('python3 CLASH.py element_region_statistic -i [--input] <Viennad_table>')
                print('python3 CLASH.py SAM_bowtie2_count -i [--input] <SAM>')
                print('python3 CLASH.py SAM_bowtie2_hOH7_single_Transcripts_count -i [--input] <SAM>')
                print('python3 CLASH.py SAM_hisat2_mRNA_count -i [--input] <SAM>')
                print('python3 CLASH.py Hisat2_SAM_gene_Table -i [--input] <SAM>')
                print('python3 CLASH.py FASTQ_length_distribution -i [--input] <FASTQ>')
                print('python3 CLASH.py FASTA_length_distribution -i [--input] <FASTA>')
                print(
                    'python3 CLASH.py FASTA_reads_number -i [--input] <FASTA> -m [--minmium_length] -M [--maxmium_length]')
                print('python3 CLASH.py reverse_complement_FASTA -i [--input] <FASTA>')
                print('python3 CLASH.py concatenated_reads -i [--input] <FASTA> -l [--fragment_length]')
                print('python3 CLASH.py hybrid_number_statistic -i [--input] <HYB> -d [--microRNA_database]')
                print('python3 CLASH.py hyb_length_distribution -i [--input] <hyb_file>')
                print('python3 CLASH.py Viennad_Rnatype_number -i [--input] <Viennad> -d [--microRNA_database]')
                print('python3 CLASH.py miRNA_count_from_Viennad -i [--input] <Viennad> -d [--microRNA_database]')
                print('python3 CLASH.py miRNA_count_from_whole_nt -i [--input] <fasta> -d [--microRNA_database]')
                print('python3 CLASH.py deduplicate_BGI -i [--input] <fastq>')
                print('python3 CLASH.py miRNA_CSV_merge -s [--string_info]')
                print ('python3 CLASH.py table_merge -s [--string_info] ')
                print('python3 CLASH.py miRNA_CSV_two_index_merge -s [--string_info]')
                print('python3 CLASH.py gene_table_two_index_merge -s [--string_info]')
                print(
                    'python3 CLASH.py all_miRNA_isoform_number_table_1st_18th -i [--input] <fasta/fastq> -d [--miRNA_database] -l [--longest_isoform_length]')
                print(
                    'python3 CLASH.py all_miRNA_isoform_table_2nd_18th -i [--input] <fasta/fastq> -d [--miRNA_database]')
                print(
                    'python3 CLASH.py all_miRNA_isoform_kind_count_1st_18th -i [--input] <fasta/fastq> -d [--miRNA_database] -l [--longest_isoform_length]')
                print(
                    'python3 CLASH.py miRNA_length_distribution -i [--input] <fasta/fastq> -d [--miRNA_database] -m [--miRNA_sequence] -p [--primiRNA_sequence]')
                print(
                    'python3 CLASH.py all_miRNA_length_distribution_1st_18th -i [--input] <fasta/fastq> -d [--miRNA_database]')
                print('python3 CLASH.py miRNA_abundance_in_hybrids -i [--input] <hybrid> ')
                print('python3 CLASH.py RNA_RNA_hybrid_distribution -i [--input] <hyb_file>')
                print(
                    'python3 CLASH.py one_kind_hybrid_distribution -i [--input] <hyb_file> -m [--miRNA_sequence] <CSV_file> -t [--target_sequence]')
                print('python3 CLASH.py making_unique_redundant_database_hg38 -i [--input] <TXT>')
                print('python3 CLASH.py making_unique_redundant_database_mm39 -i [--input] <TXT>')
                print('python3 CLASH.py making_unique_redundant_database_WBcel235 -i [--input] <TXT>')
                print('python3 CLASH.py making_simple_BedGraph -i [--input] <TXT> -c chromosome_name')
                print(
                    'python3 CLASH.py making_transcript_sequence_genomeposition_conservation_database -g [--genome_database] -t [--transcript_database] -c [--chromosome]')
                print(
                    'python3 CLASH.py Viennad_to_Table -i [--viennad] <TXT> -c [--transcript_ConservationScore_database] -t [--transcript_database] -n [--name_database]')
                print(
                    'python3 CLASH.py Combined_table_normalize_value -o [--output_name] -r [--replicate_number] -i [--input_table]')
                print(
                    'python3 CLASH.py Combined_table_raw_value -o [--output_name] -i [--input_table]')
                print('python3 CLASH.py TDMD_analyzer -i [--input]')
                print(
                    'python3 CLASH.py Cumulative_fraction_curve_targetScan -i [--DEseq_file] -a [--all_targets] -c [--conserved_targets] -b [--baseMean]')
                print(
                    'python3 CLASH.py Cumulative_fraction_curve_targetScan_CLASH -i [--DEseq_file] -a [--all_targets] -c [--conserved_targets] -l [--clash_targets] -t [--clash_targetScan_interacted] -b [--baseMean]')
                print('python3 CLASH.py fasta_Barcode_identification -b [--barcode] -i [--input_fasta]')
                print('python3 CLASH.py Sprite_SAM_to_Table -t [--Transcript_database] -i [--SAM_file]')
                print('python3 CLASH.py Sprite_SAM_distribution_RNA_type_statistic -i [--SAM_file]')
                print('python3 CLASH.py Sprite_Table_filter -i [--table_csv]')
                print('python3 CLASH.py Sprite_filtered_table_to_fasta -i [--table_csv]')
                print('python3 CLASH.py Sprite_ct_format_to_basepattern_table -i [--ct_file]')
                print('python3 CLASH.py Sprite_final_table_statistic -i [--table_csv]')
                print ('python3 CLASH.py AQ-seq -1 [fastq_1] -2 [fastq_2] -o [fastq_out] -d [miR_database]')
                print ('python3 CLASH.py fastq_each_read_count -i [--input] ')
                print ('python3 CLASH.py miR16_sfGFP_count -i [--input] ')
                print ('python3 CLASH.py Bowtie2_SAM_singleTranscript_shortReads_statistic -i [--input] ') ## if the reads do not contain hybrids,mainly count microRNA , like reads 22bp or 29bp
                print ('python3 CLASH.py Bowtie2_SAM_singleTranscript_longReads_statistic -i [--input] ')  ## if the reads contain hybrids, mainly count RNA type except miRNA, like reads 44bp

                sys.exit()
        except:
            sys.exit()


    command_line()
# os.chdir("/Users/luli/Dropbox (UFL)/Tianqi/Lu_Tianqi/CLASH all exp/Exp 371. 2024-03-06 CLASH 293T test/table")
# Combined_table(outfile='293T_TDMD_Test_20240306.txt').input_multiple_raw_clash(viennad_files='S1_S1cutUMI_comp_martquery_0228172604_181_human_unique1_add18S28S_hybrids_ua_compressed_AUex_region_TDMD.txt,'
#                                                                                              'S1_S5cutUMI_comp_martquery_0228172604_181_human_unique1_add18S28S_hybrids_ua_compressed_AUex_region_TDMD.txt,'
#                                                                                              'S2_S2cutUMI_comp_martquery_0228172604_181_human_unique1_add18S28S_hybrids_ua_compressed_AUex_region_TDMD.txt,'
#                                                                                              'S2_S6cutUMI_comp_martquery_0228172604_181_human_unique1_add18S28S_hybrids_ua_compressed_AUex_region_TDMD.txt,'
#                                                                                              'S3_S3cutUMI_comp_martquery_0228172604_181_human_unique1_add18S28S_hybrids_ua_compressed_AUex_region_TDMD.txt,'
#                                                                                              'S3_S7cutUMI_comp_martquery_0228172604_181_human_unique1_add18S28S_hybrids_ua_compressed_AUex_region_TDMD.txt,'
#                                                                                              'S4_S4cutUMI_comp_martquery_0228172604_181_human_unique1_add18S28S_hybrids_ua_compressed_AUex_region_TDMD.txt,'
#                                                                                              'S4_S8cutUMI_comp_martquery_0228172604_181_human_unique1_add18S28S_hybrids_ua_compressed_AUex_region_TDMD.txt')

