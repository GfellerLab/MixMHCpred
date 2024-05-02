#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 23 10:56:18 2023

@author: danieltadros
"""

import pandas as pd
import numpy as np

AA = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y']
def creation_PWM_dict_spec(PWMs,AAs=AA,PL=9):
    pos = [f'{i}' for i in range(PL)]
    PWMs_dict = []
    
    for PWM  in PWMs:
        PWM_D = []
        for PW in PWM:
            PWM_d = dict({})
            for i in range (len(AAs)):
                for j in range (len(pos)):
                    PWM_d[AAs[i]+pos[j]] = PW.iloc[i,j]
            PWM_D.append(PWM_d)
        PWMs_dict.append(PWM_D)

    return PWMs_dict




def interpolate_v2(x_values, y_values, x_new_values):
    x_values_reversed = x_values[::-1]
    y_values_reversed = y_values[::-1]

    idx = np.searchsorted(x_values_reversed, x_new_values, side='left') - 1
    idx = np.clip(idx, 0, len(x_values_reversed) - 2)

    x0 = x_values_reversed[idx]
    x1 = x_values_reversed[idx + 1]
    y0 = y_values_reversed[idx]
    y1 = y_values_reversed[idx + 1]

    y_new = y0 + (x_new_values - x0) * (y1 - y0) / (x1 - x0)

    # Set out-of-range values to the minimum or maximum of y_values
    y_new[x_new_values < x_values_reversed[0]] = y_values_reversed[0]
    y_new[x_new_values > x_values_reversed[-1]] = y_values_reversed[-1]

    return y_new



def ligands_PWMs_scores_spec(PWMs_dict,peptides,alphas,bias,std_dev,k,perRank):
    allele_scores = []
    allele_ranks = []
    LL = len(peptides[0])

    for i in range(len(PWMs_dict)):
        scores_corr=[]


        for x in peptides:
            score = 0
            for z in range(len(alphas[i])):
                s = 1
                for j,a in enumerate(x):
                    
                    s*=PWMs_dict[i][z][f'{a}{j}']
                # score+=alphas[i][z]*s
                score= score+(alphas[i][z]*s)


            score = np.log(score) / LL
            aa = (score - bias[i][k])/(std_dev[i][k])
            # aa = aa / std_dev[i][k] 
            scores_corr.append(aa)
            
        allele_scores.append(scores_corr)

        allele_ranks.append(interpolate_v2(perRank[i][0], perRank[i][1],scores_corr))

    return allele_scores,allele_ranks
    

def Ligands_scores_spec(ligands,ligands_L,Alleles,PWMs_dict,alphas,bias,std_dev,perRank):
    
    ligands = ligands.sort_values('length',ascending = True).reset_index()
   
    SCORES = []
    RANKS = []
    for i in range(len(ligands_L)):
        
        peptides = list(ligands[ligands['length']==ligands_L[i]]['Peptide'])
        L_i = ligands_L[i]-8
        allele_scores,allele_ranks = ligands_PWMs_scores_spec(PWMs_dict[L_i],peptides,alphas[L_i],bias,std_dev,L_i,perRank)

        SCORES.append(allele_scores)
        RANKS.append(allele_ranks)
        
    SCORES_comb = []
    RANKS_comb = []
    # for x in range(len(Alleles)):
        
    #     scores = []
    #     ranks = []
    #     for i in range(len(ligands_L)):
    #         scores = scores + list(SCORES[i][x])
    #         ranks = ranks + list(RANKS[i][x])
    #     ligands[f'Score_{Alleles[x]}'] = scores
    #     ligands[f'%Rank_{Alleles[x]}'] = ranks

    for x in range(len(Alleles)):
        
        scores = []
        ranks = []
        for i in range(len(ligands_L)):
            scores = scores + list(SCORES[i][x])
            ranks = ranks + list(RANKS[i][x])

        temp_df = pd.DataFrame({f'Score_{Alleles[x]}': scores, f'%Rank_{Alleles[x]}': ranks})
        ligands = pd.concat([ligands, temp_df], axis=1)


    ligands.insert(loc = 3,
            column = 'BestAllele',
            value = ligands[[f'%Rank_{a}' for a in Alleles]].idxmin(axis=1))

    ligands['BestAllele'] = [a[6:] for a in ligands['BestAllele']]
    ligands.insert(loc = 2,
            column = 'Score_bestAllele',
            value = ligands[[f'Score_{a}' for a in Alleles]].max(axis=1))
    ligands.insert(loc = 4,
            column = '%Rank_bestAllele',
            value = ligands[[f'%Rank_{a}' for a in Alleles]].min(axis=1))

    ligands= ligands.sort_values('index',ascending=True).reset_index(drop=True)
    ligands=ligands.drop(columns=['index','length'])
        
    return ligands


# Quality Control

def get_Allele_index(file_directory):
    indices_valid=[int(int(val)-2) for val in open(f"{file_directory}/binding_sites.txt").read().strip().split(' ')]
    Allele_Pos_idx = indices_valid
    
    return Allele_Pos_idx 


def pairwise_score(seq1, seq2, matrix):

    score = 0
    for i in range(len(seq2)):
        pair = seq1[i]+seq2[i]
        if pair in matrix:
            score += matrix[pair]
            
        else:
            raise ValueError(str(pair)+' is not in the matrix')
    return score

def best_score(x,K_seq,matrix):
    
    score,index = 0,0

    for pos, y in enumerate(K_seq):
        scores = pairwise_score(x,y,matrix)
        if scores > score :
            score,index = scores,pos
    return score,index

def sim_sequence(Alleles,K_seq,mutant_seq,matrix):
    

    score , index = best_score(mutant_seq,K_seq,matrix)
    Alleles_sim = Alleles[index]
    score =1  - (score/np.sqrt(pairwise_score(mutant_seq, mutant_seq, matrix)*pairwise_score(K_seq[index], K_seq[index], matrix)))
    
    return Alleles_sim, score 
def euclidean_distance(pwm_a, pwm_b):
    return np.sqrt(np.sum((pwm_a - pwm_b) ** 2, axis=0))


def distance_to_training(lib_path,training_Alleles,predicted_Alleles):  

    Allele_Pos_idx = get_Allele_index(f'{lib_path}/Allele_pos')
    # print(training_Alleles)
    seq_data = pd.read_csv(f"{lib_path}/MHC_I_sequences.txt", delim_whitespace= True)
    for x in training_Alleles:
        if x not in seq_data['Allele'].tolist():
            print(f'{x} Fuck')
    Alleles_seq_training = [seq_data[seq_data['Allele']== allele]['Sequence'].iloc[0] for allele in training_Alleles]

    Alleles_seq_predicted = [seq_data[seq_data['Allele']== allele]['Sequence'].iloc[0] for allele in predicted_Alleles]


    dict_load=np.load(f'{lib_path}/blosum62_update.npy', allow_pickle=True)
    blosum62=dict_load.item()    
   


    training_sequences = ["".join([allele[z] for z in Allele_Pos_idx]) for allele in Alleles_seq_training]

    predicted_sequences = ["".join([allele[z] for z in Allele_Pos_idx]) for allele in Alleles_seq_predicted]

    close_Alleles = []
    sim_scores = []


    for i in range(len(predicted_sequences)):

        Alleles_sim, score = sim_sequence(training_Alleles,training_sequences,predicted_sequences[i],blosum62)
        close_Alleles.append(Alleles_sim)
        sim_scores.append(score)

    return close_Alleles,sim_scores
