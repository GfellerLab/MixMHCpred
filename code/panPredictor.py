#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 23 10:56:18 2023

@author: danieltadros
"""

import pandas as pd
import numpy as np
import os

AA = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y']
def load_blosum(bl_path="./lib/blosum62.txt"):
    # load blosum transition matrix
    trans_prob_file=open(os.path.join(bl_path))

    dict_trans=trans_prob_file.readline().strip().split('\t')
    blosum=np.zeros((len(dict_trans),len(dict_trans)))

    for i, line in enumerate(trans_prob_file):
        blosum[i,:]=[float(val) for val in line.strip().split('\t')[1:] ]
        sum_line=sum(blosum[i,:])

        blosum[i,:]=[val/sum_line for val in blosum[i,:]] 
        
    return blosum   


def load_dictionary(file_name):
    
    with open(file_name) as f:
        amino_acids=f.readline().strip().split('\t')
        if '\n' in amino_acids[-1]:
            amino_acids[-1]=amino_acids[-1].replace("\n",'')
    
    dictionary={}
    for i_aa,aa in enumerate(amino_acids):
        dictionary[aa]=i_aa
    return dictionary


def get_Allele_pos_idx(file_directory):
    indices_valid=[int(int(val)-2) for val in open(f"{file_directory}/binding_sites.txt").read().strip().split(' ')]
    Allele_Pos_idx = indices_valid
    
    return Allele_Pos_idx    

def blosum_array(sequence,dictionary,blosum):
    dict_len=len(dictionary)
    blosum_sequence=np.zeros((len(sequence)*(dict_len)))
    for i_seq,aa in enumerate(sequence):
        try:
            blosum_sequence[i_seq*(dict_len):(i_seq+1)*(dict_len)]=blosum[dictionary[aa]]
        except Exception as e:
            print(aa + ' is not in the BLOSUM62 Matrix')

            
    return blosum_sequence

def one_hot_encode(seq):
    mapping = dict(zip('ACDEFGHIKLMNPQRSTVWY', range(20)))  
    seq2 = [mapping[i] for i in seq]
    return np.concatenate(np.eye(20)[seq2]).ravel()


def transform_allele(Allele_Pos_idx,Alleles_seq,blosum,dictionary):
    Allele_pos = []

    for allele in Alleles_seq:
        allele_sequence="".join([allele[z] for z in Allele_Pos_idx])
                    
        blos_seq = blosum_array(allele_sequence,dictionary,blosum)
        OHE_seq = one_hot_encode(allele_sequence)
        Allele_pos.append([x + y for x, y in zip(blos_seq, OHE_seq)])

    
    return Allele_pos



def softmax_manual(x, axis=-2):
    """Compute softmax values for each column"""
    e_x = np.exp(x - np.max(x, axis=axis, keepdims=True))
    return e_x / np.sum(e_x, axis=axis, keepdims=True)

def predict_pwms(X, weights, L=9):
    # Compute the PWM from the weights and the input sequence

    W1 = weights["dense_1"]["kernel"]
    b1 = weights["dense_1"]["bias"]
    W2 = weights["dense_3"]["kernel"]
    b2 = weights["dense_3"]["bias"]

    # Perform the forward pass of the neural network
    z1 = np.dot(X, W1) + b1
    a1 = np.maximum(0, z1)  

    z2 = np.dot(a1, W2) + b2 
    y_pred = softmax_manual(z2.reshape(-1, 20, L))
    pwm = y_pred.reshape(-1, 20, L)
    return pwm



def predict_pwm_list(sequences, weights, L=9):
    pwms = []
    for seq in sequences:
        pwm = predict_pwms(seq, weights, L)
        pwms.append(pwm)
    return np.array(pwms)


def predict_model(lib_path,allele_input,PL=9):
    weights_Path = f'{lib_path}/weights/weights_M{PL}.npy'

    test_allele = np.array(allele_input)

    weights = np.load(weights_Path, allow_pickle=True).item()


    PWMs = predict_pwm_list(test_allele,weights,L=PL)
    PWMs = arrays_to_pwm_dataframes(PWMs)
    
        
    return PWMs


def arrays_to_pwm_dataframes(predicted_arrays, amino_acids=AA):
    pwm_dataframes = []

    for pwm_array in predicted_arrays:
        pwm_array_2d = pwm_array.reshape((20, -1))
        pwm_df = pd.DataFrame(pwm_array_2d, index=amino_acids, columns=range(1, pwm_array_2d.shape[1] + 1))
        pwm_dataframes.append(pwm_df)

    return pwm_dataframes

def Blosum_Corr_pred(blosum_t,PWMs,AAs=AA,N=5000,PL=9):
    
    PWM_c = []
    for q in range(len(PWMs)):
        B = 200 
        PWM = PWMs[q]
        PWM_df = pd.DataFrame({})
        T = B + N 

        for i in range(PL):
            pos = []
            for j in range(20):
                x = (N * PWM[i+1][j]) 
                a = 0
                for z in range(20):
                    a+= (blosum_t[j][z] * PWM[i+1][z])
                x = x + (B*a)

                ps = x/T
                pos.append(ps)
            if (sum(pos)<0.999) | (sum(pos)>1.001):
                pos = pos/sum(pos)
            PWM_df[i+1] = pos
        PWM_df.index = AAs
        PWM_c.append(PWM_df)    
    return PWM_c 

def Freq_Corr_spec(Alleles,AA_Human_freq,PWMs):

    PWM_c = []
    for q in range(len(Alleles)):
        
        PWM = PWMs[q].copy()
            
        for i in range(20):
            PWM.iloc[i,:] = PWM.iloc[i,:]/AA_Human_freq[i]

                    
        PWM_c.append([PWM])
            
            
    return PWM_c

def PWMs_scores_spec(PWMs_dict,peptides,alphas):
    
    allele_scores = []
    LL = len(peptides[0])
    
    for i in range(len(PWMs_dict)):
        Scores=[]

        for x in peptides:
            
            score = 0
            for z in range(len(alphas[i])):
                s = 1
                for j,a in enumerate(x):
                    
                    s = s * PWMs_dict[i][z][f'{a}{j}']

                score = score + (alphas[i][z]*s)
                
            score = np.log(score) / LL
            Scores.append(score)
            
        allele_scores.append(Scores)
    return allele_scores    

def PWMs_proteome_scores_Length_spec(lib_path,Alleles,PWMs_pred_dict,L,alphas):
    proteome = []
    for l in L:
        prot = pd.read_csv(f'{lib_path}/proteome/proteome_{l}mers_100k.txt',header = None)
        proteome.append(prot[0].tolist())
        

    PWMs_pred_scores = [PWMs_scores_spec(PWMs_pred_dict[i],proteome[i],alphas[i]) for i in range(len(L))]
    for i in range(len(L)):
        for j in range(len(Alleles)):
            PWMs_pred_scores[i][j] = np.sort(PWMs_pred_scores[i][j])[::-1]
    return PWMs_pred_scores


def std_deviation(Alleles,PWMs_pred_scores,L):
    

    standard_dev=[[np.std(PWMs_pred_scores[j][i]) for j in range(len(L))] for i in range(len(Alleles))]
    
    return standard_dev

def compute_bias(Alleles,dist_len,PWMs_pred_scores,L):
    bias = []
    
    for i in range(len(Alleles)):
        dist_l = dist_len.iloc[i,1:].tolist()
        b = []
        for j in range(len(L)):
            if int(round(dist_l[j]*700)) != 0.0:
                b.append(PWMs_pred_scores[j][i][int(round(dist_l[j]*700)-1)])
            else:
                b.append(PWMs_pred_scores[j][i][0])
        bias.append(b)
    return bias





def L_predict(x, weights):
    # Extract the weights for each layer
    W1 = weights["dense_1"]["kernel"]
    b1 = weights["dense_1"]["bias"]
    W2 = weights["dense_3"]["kernel"]
    b2 = weights["dense_3"]["bias"]
    
    # Perform the forward pass of the neural network using matrix multiplication
    z1 = np.matmul(x, W1) + b1
    a1 = np.maximum(0, z1)
    z2 = np.matmul(a1, W2) + b2
    y_pred = np.exp(z2) / np.sum(np.exp(z2))
    return y_pred


def L_predict_model(allele_input,weights_Path):
    
    test_allele_Xnd = np.array(allele_input)


    weights = np.load(f'{weights_Path}/Length_weights.npy', allow_pickle=True).item()
    test_scores_p = L_predict(test_allele_Xnd,weights)
        
    return test_scores_p




def L_dist_test(Alleles,test_scores_p):
    L_dist_pred = pd.DataFrame({},columns = ['Allele','8','9','10','11','12','13','14'])
    for i in range(len(Alleles)):

        to_append = [Alleles[i]] + list(test_scores_p[i])

        L_dist_pred.loc[len(L_dist_pred)] = to_append
        
    return L_dist_pred

def prot_scores_correction(Alleles,PWMs_scores,bias,standard_dev,L):
    PWMs_scores_corr = []
    for i in range(len(Alleles)):
        PWMs_scores_corr.append([])
        for j in range(len(L)):
            a = PWMs_scores[j][i] - bias[i][j]
            a = a / standard_dev[i][j]
            PWMs_scores_corr[i] = PWMs_scores_corr[i] + a.tolist()
        PWMs_scores_corr[i] = np.sort(PWMs_scores_corr[i])[::-1]
    return PWMs_scores_corr

def perRank_function(Proteome_scores):

    y = np.concatenate([np.linspace(0.001, 0.009,num=9),np.linspace(0.01, 0.09,num=9),
                            np.linspace(0.1, 0.9,num=9),np.linspace(1, 9,num=9),
                            np.linspace(10, 100,num=10)])
    perRank = []
    for scores in Proteome_scores:
        x = [scores[int(round((i/100)*700000))] for i in y[0:-1]]
        corr = [scores[-1]]
        x = x + corr
        x = np.array(x)
        perRank.append([x, y])

    return perRank


