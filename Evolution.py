import sys
import os
import regex
import pandas as pd
import openpyxl as op
import numpy as np

from collections import Counter
import multiprocessing

from Bio.SeqIO.QualityIO import FastqGeneralIterator
from Bio.Seq import Seq
from Bio.SeqIO import write
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from Bio.Seq import reverse_complement

from openpyxl import Workbook
from openpyxl.utils.dataframe import dataframe_to_rows
from openpyxl import load_workbook
from openpyxl.styles import Font, Alignment





ind = os.getcwd()
f_list = open('target_list.txt').readlines()


first_line_list = f_list[0][:-1].split('\t')

minumum_num=int(first_line_list[1])
comparison_range=int(first_line_list[3])
indicator_range=int(first_line_list[5])
window_range=int(first_line_list[7])


WT_seq="CCTAGTTATTCACCAGTAAAAATGAATCCCTTAGCGTGCTTGGACGCGAAAATGGAATCCGTTCCGCAGATTTTCGCCGATCACGTGGCAACGCTGATTAACAAGAAGTTTGGTGAACTGATCAAAAACAACCCGGCTTATAGCCGTCGTCGTGTTTTGGCGGGTATTGTTATGTCAGAAGATTATAACTTGGAGTCTTGCAAGGTTCTGACCGTTAGCACCGGCACCAAGTGCATTAGCGGTGATCGTATTTCCGCTAATGGCAGCGTCCTGAATGATTCCCATGCTGAAATCATTTCTCGTCGCTGCCTGTGTAATTACTTCTATAGCGAGCTGGAAAAATATGTTAAAGACGACACCGGTCGTAATAGCATCTTCCTGAAAAGCCCGAACAAGGGTTTCCTGATCGATCCGAAATACAAATTCCACCTGTACATTACGAGCGCCCCTTGTGGTGACAGCCGTTTGTTCAGCCCGCATGATGTGTGTATCTTTGATAAACACCCGAATCGTAAAGCGAGAGGCCAGCTGCGTACTAAGATCGAGCGTGGTCAAGGCACGATCCCAGTTAATAGCTCCGACCTTATCCAAACCTGGGATGGTGTGTTGGTTGGCGAGCGCCTGTTGACCATGAGTTGCAGCGACAAACTGGCGAAGTGGAACGTGGTGGGCGTACAGGGTGCGTTGCTCTCTAACTTTATTCAACCGGTCTATCTGGACAGTATCATTCTGGGTAGCCTGTTCCATCCGAGCCACCTGTACCGCGCGGTGTACGGCCGCATTGAGAACGAGATCTACGACTTGACCCCACCGTATCGTCTGAACAAACCGAAGTTAGGTGCGATCATCAACACCGAAGCACGTCAAGTGGGAAAGTCGCCTTCATTTTCTATCAACTGGAATTACGGCGAGCCGGAAGCTGAGATCGTGAATGCAATAACCGGTCGTACCGATCTGGACCACATTTCCCGCATCTCGAAATGCAAGCTGTTCTACCGCTTTCTGCGCTTACTGGGCAAGATCCCGACCATTACGAACATTCATGCCACCAATGAACCGCTGATGTATAACAAGTGCAAAGAAAAGGCTACTGGTTTTCAGGCGGCGAAGATGGCACTGCTGAAAGGCTTCAGCAAAGCGAAGTTGGGTGAGTGGATGAAAAAGCCGCCGGAGGTGGATCAGTTTGAATGGGACGAGGACTCCCTGCTTAAGCACATTATCTACAACTCTTCCTATGGCTCCATTAGC"

ref_seq0="CCTAGTTATTCACCAGTAAAAATGAATCCCTTAGCGTGCTTGGACGCGAAAATGGAATCCGTTCCGCAGATTTTCGCCGATCACGTGGCAACGCTGATTAACAAGAAGTTTGGTGAACTGATCAAAAACAACCCGGCTTATAGCCGTCGTCGTGTTTTGGCGGGTATTGTTATGTCAGAAGATTAT"
ref_seq0_indicator_L="AAGAAGCGGAAAGTC"
ref_seq0_indicator_R="AACTTGGAGTCTTGC"

ref_seq1="TTATAACTTGGAGTCTTGCAAGGTTCTGACCGTTAGCACCGGCACCAAGTGCATTAGCGGTGATCGTATTTCCGCTAATGGCAGCGTCCTGAATGATTCCCATGCTGAAATCATTTCTCGTCGCTGCCTGTGTAATTACTTCTATAGCGAGCTGGAAAAATATGTTAAAGACGACACCGGTCGTAATAGCATCTTCCTGAAAAGCC"
ref_seq1_indicator_L="TGTTATGTCAGAAGA"
ref_seq1_indicator_R="CGAACAAGGGTTTCC"

ref_seq2="AAAAGCCCGAACAAGGGTTTCCTGATCGATCCGAAATACAAATTCCACCTGTACATTACGAGCGCCCCTTGTGGTGACAGCCGTTTGTTCAGCCCGCATGATGTGTGTATCTTTGATAAACACCCGAATCGTAAAGCGAGAGGCCAGCTGCGTACTAAGATCGAGCGTGGTCAAGGCACGATCCCAGTTAATAGCTCCGACCTTAT"
ref_seq2_indicator_L="AATAGCATCTTCCTG"
ref_seq2_indicator_R="CCAAACCTGGGATGG"

ref_seq3="GTGTGTATCTTTGATAAACACCCGAATCGTAAAGCGAGAGGCCAGCTGCGTACTAAGATCGAGCGTGGTCAAGGCACGATCCCAGTTAATAGCTCCGACCTTATCCAAACCTGGGATGGTGTGTTGGTTGGCGAGCGCCTGTTGACCATGAGTTGCAGCGACAAACTGGCGAAGTGGAACGTGGTGGGCGTACAGGGTGCGTTGCTCT"
ref_seq3_indicator_L="TTCAGCCCGCATGAT"
ref_seq3_indicator_R="CTAACTTTATTCAAC"

ref_seq4="CTCTCTAACTTTATTCAACCGGTCTATCTGGACAGTATCATTCTGGGTAGCCTGTTCCATCCGAGCCACCTGTACCGCGCGGTGTACGGCCGCATTGAGAACGAGATCTACGACTTGACCCCACCGTATCGTCTGAACAAACCGAAGTTAGGTGCGATCATCAACACCGAAGCACGTCAAGTGGGAAAGTCGCCTTCATTTTCTATCAACT"
ref_seq4_indicator_L="GTACAGGGTGCGTTG"
ref_seq4_indicator_R="GGAATTACGGCGAGC"

ref_seq5="GGAAAGTCGCCTTCATTTTCTATCAACTGGAATTACGGCGAGCCGGAAGCTGAGATCGTGAATGCAATAACCGGTCGTACCGATCTGGACCACATTTCCCGCATCTCGAAATGCAAGCTGTTCTACCGCTTTCTGCGCTTACTGGGCAAGATCCCGACCATTACGAACATTCATGCCACCAATGAACCGCTGATGTATAACAAGTGCAAAGAAA"
ref_seq5_indicator_L="GAAGCACGTCAAGTG"
ref_seq5_indicator_R="AGGCTACTGGTTTTC"

ref_seq6="CGCTGATGTATAACAAGTGCAAAGAAAAGGCTACTGGTTTTCAGGCGGCGAAGATGGCACTGCTGAAAGGCTTCAGCAAAGCGAAGTTGGGTGAGTGGATGAAAAAGCCGCCGGAGGTGGATCAGTTTGAATGGGACGAGGACTCCCTGCTTAAGCACATTATCTACAACTCTTCCTATGGCTCCATTAGCTCTGGAGGATCTAG"
ref_seq6_indicator_L="ATGCCACCAATGAAC"
ref_seq6_indicator_R="CGGAGGATCCTCTGG"




# minimum_frequency : Time Complexity O(1)

def minimum_frequency(records, n=2):
    seq_form_fastq = [record[1] for record in records]  
    counts = Counter(seq_form_fastq)    
    minimum_records = [record for record in records if counts[record[1]] >= n] 
    return minimum_records 



# indicator selection : w/ mismatch

def pick_indicator_with_mismatch(seq_input, left_indicator, right_indicator): # seq_input type = minimum_records
    indicator_records=[]
    for seq in seq_input:
        is_WT_L = regex.search(rf"({left_indicator}){{s<={2}}}", seq[1]) # regex fuzzy matching LI
        is_WT_R = regex.search(rf"({right_indicator}){{s<={2}}}", seq[1]) # regex fuzzy matching RI
        if is_WT_L is not None and is_WT_R is not None:
            nomean=['A','B','C']
            nomean[1]=seq[1][is_WT_L.end():is_WT_R.start()]
            indicator_records.append(nomean)

    return indicator_records 





# count mutagen

def count_mutagen(seq_input, ref_seq): # seq_input type = indicator_records
    sigle_read_muts = {f'{i} mut': 0 for i in range(1, 11)}
    sub_rate = [[0,0,0,0, ref_seq[i]] for i in range(len(ref_seq))] #ATGC
    cnt = {'dimer' : 0, 'ins' : 0, 'del' : 0, 'WT' :0, 'mutants': 0, 'mutation(nucleotide)_number' : 0}
    
    
    for seq in seq_input: 
        if len(seq[1]) < len(ref_seq)*0.4:
            cnt['dimer'] += 1

        else:
            alignments = pairwise2.align.globalms(ref_seq, seq[1], 1, -1, -10, -.5)
            split_format = format_alignment(*alignments[0]).split('\n')
            alignments_list = [s.strip() for s in split_format]

            if '-' in alignments_list[0]:
                cnt['ins'] += 1
                
            elif '-' in alignments_list[2]:
                cnt['del'] += 1
                
            elif '.' not in alignments_list[1]:
                cnt['WT'] += 1
                
            else:
                cnt['mutants'] += 1
                num_mutations = alignments_list[1].count(".")
                cnt['mutation(nucleotide)_number'] += num_mutations
                mut_positions = tuple((index) for index, char in enumerate(alignments_list[1]) if char == '.')
                if num_mutations >= 10:
                    sigle_read_muts['10 mut'] += 1
                else:
                    sigle_read_muts[f'{num_mutations} mut'] += 1
                
                for i in mut_positions:
                    if alignments_list[2][i] == 'A':
                        sub_rate[i][0] += 1
                    elif alignments_list[2][i] == 'T':
                        sub_rate[i][1] += 1
                    elif alignments_list[2][i] == 'G':
                        sub_rate[i][2] += 1
                    elif alignments_list[2][i] == 'C':
                        sub_rate[i][3] += 1
    return sub_rate, cnt, sigle_read_muts
    
   








### 기록 시작부분 ###








wb = Workbook()

thisis1=0

with pd.ExcelWriter('sub.xlsx', engine='openpyxl', mode='w') as writer:
    writer.book = wb
    
    
    for t in f_list[1:]:    
        each_line_list0 = t.split('\t')
        each_line_list = [i.strip() for i in each_line_list0]
        
            









    

        direc_j = each_line_list[0]       
        f0_index = each_line_list[1]    
        f1_index = each_line_list[2]
        f2_index = each_line_list[3]
        f3_index = each_line_list[4]
        f4_index = each_line_list[5]
        f5_index = each_line_list[6]
        f6_index = each_line_list[7]
        


        direc = str(direc_j)




    






    
    
    

        f0_file = str(f0_index) + '.fastqjoin'
        f1_file = str(f1_index) + '.fastqjoin'
        f2_file = str(f2_index) + '.fastqjoin'
        f3_file = str(f3_index) + '.fastqjoin'  
        f4_file = str(f4_index) + '.fastqjoin'
        f5_file = str(f5_index) + '.fastqjoin'
        f6_file = str(f6_index) + '.fastqjoin'


        os.chdir(ind + '/' + direc)


    




    



        file_names = [f0_file, f1_file, f2_file, f3_file, f4_file, f5_file, f6_file]
        print(file_names)
        ref_seq_list = [ref_seq0, ref_seq1, ref_seq2, ref_seq3, ref_seq4, ref_seq5, ref_seq6]
        
        

        indicator_group = (
            (ref_seq0_indicator_L, ref_seq0_indicator_R), 
            (ref_seq1_indicator_L, ref_seq1_indicator_R), 
            (ref_seq2_indicator_L, ref_seq2_indicator_R), 
            (ref_seq3_indicator_L, ref_seq3_indicator_R), 
            (ref_seq4_indicator_L, ref_seq4_indicator_R), 
            (ref_seq5_indicator_L, ref_seq5_indicator_R), 
            (ref_seq6_indicator_L, ref_seq6_indicator_R)
        )
        
        records = [[] for _ in range(7)]



        for file_name, empty_record in zip(file_names, records):
            with open(file_name, "r") as handle:
                for record in FastqGeneralIterator(handle):
                    empty_record.append(record)
                    

                        

        os.chdir(ind)



        record_minimum = [[] for _ in range(7)]
        record_indicator = [[] for _ in range(7)]
        record_mutagen = [[] for _ in range(7)]
        cnt_list = [{} for _ in range(7)]
        single_read_muts_list = [{} for _ in range(7)]
        raw_sub_rate_list = [[] for _ in range(7)]
        
        A_sub_rate_percentage_list =[[] for _ in range(7)]
        T_sub_rate_percentage_list =[[] for _ in range(7)]
        G_sub_rate_percentage_list =[[] for _ in range(7)]
        C_sub_rate_percentage_list =[[] for _ in range(7)]
        total_sub_rate_percentage_list =[[] for _ in range(7)]
                
        
        
        

        # 0's calculator
        
        
        for idx, (minimum_X, indicator_Y, mutagen_Z, input_K, file_name_N, (left_indicator, right_indicator), ref_seq_M, raw_sub, A_sub, T_sub, G_sub, C_sub, total_sub) in enumerate(zip(
                record_minimum, record_indicator, record_mutagen, records, file_names, indicator_group, ref_seq_list, raw_sub_rate_list,
                A_sub_rate_percentage_list, T_sub_rate_percentage_list, G_sub_rate_percentage_list, C_sub_rate_percentage_list, total_sub_rate_percentage_list)):
            
            
            minimum_X=minimum_frequency(input_K, n=minumum_num)
            indicator_Y=pick_indicator_with_mismatch(minimum_X, left_indicator, right_indicator)
            
   
            print(f"#{file_name_N} : indicator matched = {len(indicator_Y)}")
            
            if len(indicator_Y) == 0:
                error_message = f"\u25B6 # {file_name_N} has no read."
                print(error_message)
                print("--------------------------------------------------------------------------")
                
                alternative_sigle_read_muts = {f'{i} mut': 0 for i in range(1, 11)}
                alternative_sub_rate = [[0,0,0,0, ref_seq_M[i]] for i in range(len(ref_seq_M))] #ATGC
                alternative_cnt = {'dimer' : 0, 'ins' : 0, 'del' : 0, 'WT' :0, 'mutants': 0, 'mutation(nucleotide)_number' : 0}
                A_sub.extend([0 for i in range(len(sub_rate))]) # A_rate
                T_sub.extend([0 for i in range(len(sub_rate))]) # T_rate
                G_sub.extend([0 for i in range(len(sub_rate))]) # G_rate
                C_sub.extend([0 for i in range(len(sub_rate))]) # C_rate
                total_sub.extend([0 for i in range(len(sub_rate))]) # %
                
                cnt_list[idx] = alternative_cnt
                single_read_muts_list[idx] = alternative_sigle_read_muts
                sub_rate = alternative_sub_rate
                
                
                continue
                ##raise RuntimeError(error_message)

            
            else:        
                mutagen_Z = count_mutagen(indicator_Y, ref_seq_M)
                sub_rate = mutagen_Z[0]
                raw_sub.extend(sub_rate)
                
                cnt_list[idx] = mutagen_Z[1]
                
                single_read_muts_list[idx] = mutagen_Z[2]
                
                mutation_number = mutagen_Z[1]['mutants']+mutagen_Z[1]['WT']
            
                if mutation_number != 0:
                    A_sub.extend([sub_rate[i][0] / mutation_number * 100 for i in range(len(sub_rate))]) # A_rate
                    T_sub.extend([sub_rate[i][1] / mutation_number * 100 for i in range(len(sub_rate))]) # T_rate
                    G_sub.extend([sub_rate[i][2] / mutation_number * 100 for i in range(len(sub_rate))]) # G_rate
                    C_sub.extend([sub_rate[i][3] / mutation_number * 100 for i in range(len(sub_rate))]) # C_rate
                    total_sub.extend([sum(sub_rate[i][:4]) / mutation_number * 100 for i in range(len(sub_rate))]) # %
                    print(mutagen_Z[1])
                    print("--------------------------------------------------------------------------")
                    
                    
                else:
                    print(f"\u25B6 # {mutation_number} is Zero.")
                    A_sub.extend([0 for i in range(len(sub_rate))]) # A_rate
                    T_sub.extend([0 for i in range(len(sub_rate))]) # T_rate
                    G_sub.extend([0 for i in range(len(sub_rate))]) # G_rate
                    C_sub.extend([0 for i in range(len(sub_rate))]) # C_rate
                    total_sub.extend([0 for i in range(len(sub_rate))]) # %
                    print("--------------------------------------------------------------------------")

                    continue



        ref_seqs = [
            "CCTAGTTATTCACCAGTAAAAATGAATCCCTTAGCGTGCTTGGACGCGAAAATGGAATCCGTTCCGCAGATTTTCGCCGATCACGTGGCAACGCTGATTAACAAGAAGTTTGGTGAACTGATCAAAAACAACCCGGCTTATAGCCGTCGTCGTGTTTTGGCGGGTATTGTTATGTCAGAAGATTAT",
            "TTATAACTTGGAGTCTTGCAAGGTTCTGACCGTTAGCACCGGCACCAAGTGCATTAGCGGTGATCGTATTTCCGCTAATGGCAGCGTCCTGAATGATTCCCATGCTGAAATCATTTCTCGTCGCTGCCTGTGTAATTACTTCTATAGCGAGCTGGAAAAATATGTTAAAGACGACACCGGTCGTAATAGCATCTTCCTGAAAAGCC",
            "AAAAGCCCGAACAAGGGTTTCCTGATCGATCCGAAATACAAATTCCACCTGTACATTACGAGCGCCCCTTGTGGTGACAGCCGTTTGTTCAGCCCGCATGATGTGTGTATCTTTGATAAACACCCGAATCGTAAAGCGAGAGGCCAGCTGCGTACTAAGATCGAGCGTGGTCAAGGCACGATCCCAGTTAATAGCTCCGACCTTAT",
            "GTGTGTATCTTTGATAAACACCCGAATCGTAAAGCGAGAGGCCAGCTGCGTACTAAGATCGAGCGTGGTCAAGGCACGATCCCAGTTAATAGCTCCGACCTTATCCAAACCTGGGATGGTGTGTTGGTTGGCGAGCGCCTGTTGACCATGAGTTGCAGCGACAAACTGGCGAAGTGGAACGTGGTGGGCGTACAGGGTGCGTTGCTCT",
            "CTCTCTAACTTTATTCAACCGGTCTATCTGGACAGTATCATTCTGGGTAGCCTGTTCCATCCGAGCCACCTGTACCGCGCGGTGTACGGCCGCATTGAGAACGAGATCTACGACTTGACCCCACCGTATCGTCTGAACAAACCGAAGTTAGGTGCGATCATCAACACCGAAGCACGTCAAGTGGGAAAGTCGCCTTCATTTTCTATCAACT",
            "GGAAAGTCGCCTTCATTTTCTATCAACTGGAATTACGGCGAGCCGGAAGCTGAGATCGTGAATGCAATAACCGGTCGTACCGATCTGGACCACATTTCCCGCATCTCGAAATGCAAGCTGTTCTACCGCTTTCTGCGCTTACTGGGCAAGATCCCGACCATTACGAACATTCATGCCACCAATGAACCGCTGATGTATAACAAGTGCAAAGAAA",
            "CGCTGATGTATAACAAGTGCAAAGAAAAGGCTACTGGTTTTCAGGCGGCGAAGATGGCACTGCTGAAAGGCTTCAGCAAAGCGAAGTTGGGTGAGTGGATGAAAAAGCCGCCGGAGGTGGATCAGTTTGAATGGGACGAGGACTCCCTGCTTAAGCACATTATCTACAACTCTTCCTATGGCTCCATTAGC"
        ]


        full_length = len(WT_seq)
        full_A_sub_rate_percentage = [0 for _ in range(full_length)]
        full_T_sub_rate_percentage = [0 for _ in range(full_length)]
        full_G_sub_rate_percentage = [0 for _ in range(full_length)]
        full_C_sub_rate_percentage = [0 for _ in range(full_length)]
        full_sum_sub_rate_percentage = [0 for _ in range(full_length)]
        full_sub_rate = [[0,0,0,0, WT_seq[i]] for i in range(len(WT_seq))]

        counts = [0 for _ in range(full_length)]  

        




        full_A_sub_rate_sum = [0 for _ in range(full_length)]
        full_T_sub_rate_sum = [0 for _ in range(full_length)]
        full_G_sub_rate_sum = [0 for _ in range(full_length)]
        full_C_sub_rate_sum = [0 for _ in range(full_length)]
        full_sum_sub_rate_sum = [0 for _ in range(full_length)]
        full_sub_rate_sum = [[0, 0, 0, 0, WT_seq[i]] for i in range(full_length)]
        counts = [0 for _ in range(full_length)]  






        
        
        avg_full_A_sub_rate_percentage = [0 for _ in range(full_length)]
        avg_full_T_sub_rate_percentage = [0 for _ in range(full_length)]
        avg_full_G_sub_rate_percentage = [0 for _ in range(full_length)]
        avg_full_C_sub_rate_percentage = [0 for _ in range(full_length)]
        avg_full_sum_sub_rate_percentage = [0 for _ in range(full_length)]
        avg_full_sub_rate = [[0,0,0,0, WT_seq[i]] for i in range(len(WT_seq))]
        
        
        
        
        full_cnt = {'dimer': 0, 'ins': 0, 'del': 0, 'WT': 0, 'mutants': 0, 'mutation(nucleotide)_number': 0}
        full_single_read_muts = {f'{i} mut': 0 for i in range(1, 11)}
        
                
        start_positions = [
            WT_seq.find(ref_seqs[0][:15]),
            WT_seq.find(ref_seqs[1][:15]),
            WT_seq.find(ref_seqs[2][:15]),
            WT_seq.find(ref_seqs[3][:15]),
            WT_seq.find(ref_seqs[4][:15]),
            WT_seq.find(ref_seqs[5][:15]),
            WT_seq.find(ref_seqs[6][:15])
        ]




        for k, (fragment, start_pos, cnt, single_read_muts, A_sub, T_sub, G_sub, C_sub, total_sub, raw_sub) in enumerate(zip(
                ref_seqs, start_positions, cnt_list, single_read_muts_list,
                A_sub_rate_percentage_list, T_sub_rate_percentage_list, G_sub_rate_percentage_list, C_sub_rate_percentage_list, total_sub_rate_percentage_list, raw_sub_rate_list)):
            
            
            if start_pos < 0 or start_pos >= full_length:
                print(f"Invalid start position for fragment {k}: {start_pos}")
                continue



            for i, (sub_rate, a_rate, t_rate, g_rate, c_rate, total_rate) in enumerate(zip(raw_sub, A_sub, T_sub, G_sub, C_sub, total_sub)):
                if start_pos + i < full_length:
                    counts[start_pos + i] += 1  
                    full_A_sub_rate_sum[start_pos + i] += a_rate
                    full_T_sub_rate_sum[start_pos + i] += t_rate
                    full_G_sub_rate_sum[start_pos + i] += g_rate
                    full_C_sub_rate_sum[start_pos + i] += c_rate
                    full_sum_sub_rate_sum[start_pos + i] += total_rate

                    for j in range(4):
                        full_sub_rate_sum[start_pos + i][j] += sub_rate[j]


            
            for key in full_cnt:
                full_cnt[key] += cnt[key]

            for key in full_single_read_muts:
                full_single_read_muts[key] += single_read_muts[key]
                        

        
        for i in range(full_length):
            if counts[i] > 0:  
                full_A_sub_rate_sum[i] /= counts[i]
                full_T_sub_rate_sum[i] /= counts[i]
                full_G_sub_rate_sum[i] /= counts[i]
                full_C_sub_rate_sum[i] /= counts[i]
                full_sum_sub_rate_sum[i] /= counts[i]

                for j in range(4):
                    full_sub_rate_sum[i][j] /= counts[i]

        
        full_A_sub_rate_percentage = full_A_sub_rate_sum
        full_T_sub_rate_percentage = full_T_sub_rate_sum
        full_G_sub_rate_percentage = full_G_sub_rate_sum
        full_C_sub_rate_percentage = full_C_sub_rate_sum
        full_sum_sub_rate_percentage = full_sum_sub_rate_sum
        full_sub_rate = full_sub_rate_sum

        
        
        
        
        
        
        
        Most_chaned_base = []
        
        for a, t, g , c in zip(full_A_sub_rate_percentage, full_T_sub_rate_percentage, full_G_sub_rate_percentage, full_C_sub_rate_percentage):
            if a > t and a > g and a > c:
                Most_chaned_base.append('A')
            elif t > a and t > g and t > c:
                Most_chaned_base.append('T')
            elif g > a and g > t and g > c:
                Most_chaned_base.append('G')
            elif c > a and c > t and c > g:
                Most_chaned_base.append('C')
            else:
                Most_chaned_base.append('N')
        
        
        changed_base_data = [f'{wt_base}({i+1}){changed_base}' for i, (changed_base, wt_base)  in enumerate(zip(Most_chaned_base, WT_seq))]

                    
                    




        changed_aa = []

        for mutation in changed_base_data:
            original_base = mutation[0]
            changed_base = mutation[-1]
            
            start = mutation.find('(') + 1
            end = mutation.find(')')
            position_str = mutation[start:end]

            position = int(position_str) - 1  
            
            codon_start = position - (position % 3)
            codon_end = codon_start + 3
            aa_position = position // 3 + 254

            WT_seq_codon = WT_seq[codon_start:codon_end]

            
            if position % 3 == 0:
                new_seq = changed_base + WT_seq_codon[1:]
            elif position % 3 == 1:
                new_seq = WT_seq_codon[0] + changed_base + WT_seq_codon[2]
            else:
                new_seq = WT_seq_codon[:2] + changed_base
                
            old_aa = Seq(WT_seq_codon).translate()
            if 'N' == changed_base:
                new_aa = 'X'
            else:
                new_aa = Seq(new_seq).translate()
            
            new_info = old_aa + str(aa_position) + new_aa
            
            changed_aa.append(new_info)













        data_name= {
            'Files' : [item for item in file_names]
        }
            

            
        data0 = {
            'title': list(full_cnt.keys()),
            '1st read': list(cnt_list[0].values()),
            '2nd read': list(cnt_list[1].values()),
            '3rd read': list(cnt_list[2].values()),
            '4th read': list(cnt_list[3].values()),
            '5th read': list(cnt_list[4].values()),
            '6th read': list(cnt_list[5].values()),
            '7th read': list(cnt_list[6].values()),              
            'total': list(full_cnt.values())
        }





        data = {
            'Seq': [f'{base}({i+1})' for i, base in enumerate(WT_seq)],
            'A rates': full_A_sub_rate_percentage,
            'T rates': full_T_sub_rate_percentage,
            'G rates': full_G_sub_rate_percentage,
            'C rates': full_C_sub_rate_percentage,
            'total(sum)': full_sum_sub_rate_percentage,
            'changed_base_data' : changed_base_data,
            'changed_aa' : changed_aa

        }                
    

        data1 = {
            'title': list(full_single_read_muts.keys()),
            '1st read': list(single_read_muts_list[0].values()),
            '2nd read': list(single_read_muts_list[1].values()),
            '3rd read': list(single_read_muts_list[2].values()),
            '4th read': list(single_read_muts_list[3].values()),
            '5th read': list(single_read_muts_list[4].values()),
            '6th read': list(single_read_muts_list[5].values()),
            '7th read': list(single_read_muts_list[6].values()),
            'contents': list(full_single_read_muts.values())
        }



        df_name = pd.DataFrame(data_name)
        df0 = pd.DataFrame(data0)
        df = pd.DataFrame(data)
        df1 = pd.DataFrame(data1) 

        df['A rates'] = df['A rates'].round(2)
        df['T rates'] = df['T rates'].round(2)
        df['G rates'] = df['G rates'].round(2)
        df['C rates'] = df['C rates'].round(2)
        df['total(sum)'] = df['total(sum)'].round(2)

        df_name_transposed = df_name.transpose()
        df0_transposed = df0.transpose()
        df_transposed = df.transpose()
        df1_transposed = df1.transpose()





        # raw_sub_rate 
        
        

        base_index = {'A': 0, 'T': 1, 'G': 2, 'C': 3}

        # 4x4 matrix
        matrix = np.zeros((4, 4), dtype=float)

        
        ## print (full_sub_rate_percentage)
        
        # sub_rate
        for i in range(len(full_sub_rate)):
            base_from = full_sub_rate[i][4]
            for j in range(4):
                base_to = 'ATGC'[j]
                matrix[base_index[base_from]][j] += full_sub_rate[i][j]

        # mutation_number

        matrix_percentage = matrix / full_cnt['mutation(nucleotide)_number'] * 100


        # Pandas DataFrame
        columns = ['A', 'T', 'G', 'C']
        index = ['\u25B6A', '\u25B6T', '\u25B6G', '\u25B6C']

        df_matrix = pd.DataFrame(matrix_percentage.T, columns=columns, index=index).round(2)
        df_matrix_transposed=df_matrix.transpose()




        thisis1 += 1
        
        df_name_transposed.to_excel(writer, sheet_name=str(f0_file), startrow=0, startcol=0, header=False)
        df_matrix_transposed.to_excel(writer, sheet_name=str(f0_file), startrow=2, startcol=0, header=True)
        df0_transposed.to_excel(writer, sheet_name=str(f0_file), startrow=8, startcol=0, header=False)
        df1_transposed.to_excel(writer, sheet_name=str(f0_file), startrow=18, startcol=0, header=False)
        df.to_excel(writer, sheet_name=str(f0_file), startrow=28, startcol=0, header=True, index=False)
        
        

        
        
        print(f"\u25B6 Job No. #{thisis1} analysis completed.")
        print("--------------------------------------------------------------------------")

print('Jobs done! mady by HW')




                        
                    