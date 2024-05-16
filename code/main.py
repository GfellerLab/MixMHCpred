from argparse import ArgumentParser
import numpy as np
import pandas as pd 
from MixMHCpred import *
import warnings
warnings.simplefilter(action='ignore', category=Warning)
import os
import shutil
parser=ArgumentParser()


parser.add_argument("-i", "--InputFilePath",help="file path of the peptidess")
parser.add_argument("-o", "--OutputFilePath",help="file path where to save results")
parser.add_argument("-a","--Alleles",help="Alleles to compute i.e A0101,B0702,C0107")
parser.add_argument("-l","--lib_path",help="Library directory")
parser.add_argument("-m","--output_motifs",help="create the binding motifs",default='0')
parser.add_argument("-p","--score_peptides",help="to compute binding scores",default='0')

args=parser.parse_args()

file_input=args.InputFilePath
output_dir = args.OutputFilePath
file_output = args.OutputFilePath
lib_path = args.lib_path
output_motifs = int(args.output_motifs)
Alleles_to_test=args.Alleles

if output_motifs == 0:
    if os.path.isdir(file_output):
        print(f'''
            ERROR: A folder named '{file_output}' already exists.

            Please choose a different file name.
        ''')
        exit()
    # if file_output[0:2] != './':
    #     file_output = f'./{file_output}'
    # # Extract the directory part from the file_output path
    # directory = os.path.dirname(file_output)

    # # Check if the directory exists
    # if not os.path.exists(directory):
    #     print(f'''cd 
    #     The directory {directory} does not exist. Creating it now.
    #     ''')
    #     os.makedirs(directory)


if output_motifs==1:
    if os.path.isfile(file_output):
        print(f'''
        ERROR: A file named '{file_output}' already exists.
        
        Please choose a different folder name to save multiple files in it.
        ''')
       
        
        exit()

    # Check if the output directory exists
    if os.path.exists(output_dir):
        # Remove it if it exists
        shutil.rmtree(output_dir)

    # Create the output directory
    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(f'{output_dir}/Motifs', exist_ok=True)


    import matplotlib.pyplot as plt
    from matplotlib.ticker import MaxNLocator

    import logomaker 
    plt.style.use('ggplot')   
    chemistry_colors = {
        'G': '#109648', 'S': '#109648', 'T': '#109648', 'Y': '#109648', 'C': '#109648', 
        'N': '#5E239D', 'Q': '#5E239D', 
        'K': '#255C99', 'R': '#255C99', 'H': '#255C99', 
        'D': '#D62839', 'E': '#D62839', 
        'P': '#221E22', 'A': '#221E22', 'W': '#221E22', 'F': '#221E22', 
        'L': '#221E22', 'I': '#221E22', 'M': '#221E22', 'V': '#221E22'
    }


def Convert(string):
    li = list(string.split(","))
    return li


Alleles_input = Convert(Alleles_to_test)

Alleles_to_test = []
for h in Alleles_input:
    if h.startswith('HLA-'):
        h = h[4:]
    if h[1] == '*':
        h = h[0] + h[2:]
    if h[3] == ':':
        h = h[:3] + h[4:]
    if h.startswith('H-2'):
        h = 'H2' + h[3:]
    Alleles_to_test.append(h)

del Alleles_input

Alleles_to_test = [s.replace(':', '') for s in Alleles_to_test]
Alleles_to_test = [s.replace('*', '') for s in Alleles_to_test]


L = [8,9,10,11,12,13,14]
alleles_list = pd.read_csv(f'{lib_path}/alleles_list.txt',sep='\t')

# peptides = pd.read_csv(file_input ,header = None)[0].to_list()
# Check input file
if os.path.isfile(file_input) == False:
    print(f'''
    ERROR: No such file : {file_input} exists.
    ''')
    exit()
else:
    with open(file_input, 'r') as f:
        peptide = [line.strip() for line in f if not line.startswith(">") and line.strip()]

if not peptide:
    print("Empty peptide file")
    exit()

Lmin, Lmax = 8, 14

amino_acids_map = {
    "A": 0, "C": 1, "D": 2, "E": 3, "F": 4, "G": 5, "H": 6, "I": 7,
    "K": 8, "L": 9, "M": 10, "N": 11, "P": 12, "Q": 13, "R": 14, "S": 15,
    "T": 16, "V": 17, "W": 18, "Y": 19
}


peptide_num = []

for p in peptide:
    try:
        peptide_num.append([amino_acids_map[char] for char in p])
    except KeyError as e:
        print(f"Error in input sequence - unknown amino acids {e.args[0]}: {p}")
        exit()

    peptide_len = len(p)
    if not Lmin <= peptide_len <= Lmax:
        print(f"Incompatible peptide length: {p}\t{peptide_len}. \nOnly peptides of length {Lmin}-{Lmax} are supported")
        exit()
    


del peptide_num
del peptide_len


Ligands = pd.DataFrame({'Peptide':peptide})
del peptide
Ligands['length'] = Ligands.Peptide.str.len()
Ligands_L = np.unique(Ligands['length'])

Alleles = alleles_list['Allele'].tolist()




PWMs_pred = []
PWMs_pred_dict = []
PWMs_pred_dict_ = []
Alleles_perRank_f = []
bias = []
standard_dev = []

Alleles_in = np.sort(list(set(Alleles_to_test) & set(Alleles))).tolist()
Alleles_out =  np.sort(list(set(Alleles_to_test) - set(Alleles))).tolist()

alphas = []
alphas_ = []
if len(Alleles_in)>0:


    for l in L:
        ratios = pd.read_csv(f'{lib_path}/pwm/class1_{l}/alphas.txt',sep='\t')
        PWMs_pred.append([])
        alphas.append([])

        for xxx in Alleles_in:
            PWM = []
            alpha = []
            n = alleles_list[alleles_list['Allele']==xxx][f'{l}'].iloc[0]
            for i in range(n):
                PWM.append(pd.read_csv(f'{lib_path}/pwm/class1_{l}/PWM_{xxx}_{i+1}.csv',index_col=0))
                alpha.append(ratios[ratios['Allele']==f'{xxx}_{i+1}']['ratio'].iloc[0])
            PWMs_pred[-1].append(PWM)
            alphas[-1].append(alpha)                                                                                                                              

        PWMs_pred_dict.append(creation_PWM_dict_spec(PWMs_pred[-1],AAs=AA,PL=l))

    Alleles_perRank= [pd.read_csv(f'{lib_path}/PerRank/{zz}.txt',sep='\t') for zz in Alleles_in]

    Alleles_perRank_f = Alleles_perRank_f + [[Alleles_perRank[i]['score'].to_numpy(), Alleles_perRank[i]['rank'].to_numpy()] for i in range(len(Alleles_in))]
    bias = bias + pd.read_csv(f'{lib_path}/shifts/bias.txt',sep='\t',index_col=0).loc[Alleles_in].to_numpy().tolist()
    standard_dev = standard_dev + pd.read_csv(f'{lib_path}/shifts/standard_dev.txt',sep='\t',index_col=0).loc[Alleles_in].to_numpy().tolist()


if len(Alleles_out)>0:
    
    from panPredictor import *
    
    seq_data = pd.read_csv(f"{lib_path}/MHC_I_sequences.txt", delim_whitespace= True)
    Alleles_ = seq_data['Allele'].tolist()
    # Alleles_seq = list(seq_data['Sequence'])
    # Alleles_H = Alleles_

    invalid_alleles = set(Alleles_out) - set(Alleles_)

    if invalid_alleles:
        print("Please make sure you have the correct name for the following alleles:")
        for allele in invalid_alleles:
            print(allele)
        exit()

    blosum = load_blosum(bl_path=f'{lib_path}/blosum62.txt')
    dictionary=load_dictionary(f'{lib_path}/blosum62.txt')
    blosum_t = np.transpose(blosum)
    AA_Human_freq = list(pd.read_csv(f'{lib_path}/proteome/AA_frequency_human_Methionine.csv')['Frequency'])
    
    Alleles_seq_test = [seq_data[seq_data['Allele']== allele]['Sequence'].iloc[0] for allele in Alleles_out]

    
    Allele_Pos_idx = get_Allele_pos_idx(f'{lib_path}/Allele_pos')
    Allele_pos = transform_allele(Allele_Pos_idx,Alleles_seq_test,blosum,dictionary)


    for l in L:
    
        alphas_.append([])

        PWMs_test = predict_model(lib_path,Allele_pos,PL=l)
        PWMS_pred_bl = Blosum_Corr_pred(blosum_t,PWMs_test,AAs=AA,N=5000,PL=l)

        PWMs_pred = Freq_Corr_spec(Alleles_out,AA_Human_freq,PWMS_pred_bl)
        if output_motifs==1:
            print(f'Creating Motifs for length {l}')
            for i in range(len(Alleles_out)):
                
                
                PWMs_pred[i][0].to_csv(f'{output_dir}/Motifs/{Alleles_out[i]}_PWM{l}.txt',sep='\t')

                plt.ioff()

                pwm_df = PWMs_pred[i][0].T
                pwm_df = pwm_df.div(pwm_df.sum(axis=1), axis=0)
                pssm_df = logomaker.transform_matrix(pwm_df, 
                                                    from_type='probability', 
                                                    to_type='information')

                fig, ax = plt.subplots(figsize=(4, 4))
                fig.set_facecolor('white'),ax.set_facecolor('white')
                ax.spines['top'].set_visible(False),ax.spines['right'].set_visible(False)
                ax.spines['left'].set_edgecolor('grey'),ax.spines['bottom'].set_edgecolor('grey')
                logo = logomaker.Logo(pssm_df,color_scheme=chemistry_colors,ax=ax)
                logo.ax.set_ylim([0,4])
                positions = range(1, len(pwm_df) + 1)
                logo.ax.set_xticks(positions),logo.ax.set_xticklabels(positions)
                logo.ax.grid(False)
                ax.yaxis.set_major_locator(MaxNLocator(integer=True))
                fig.savefig(f'{output_dir}/Motifs/{Alleles_out[i]}_PWM{l}.png', bbox_inches='tight',dpi=300)

        alphas_[-1] = [[1.0] for t in range(len(Alleles_out))]
        PWMs_pred_dict_.append(creation_PWM_dict_spec(PWMs_pred,AAs=AA,PL=l))


    prot_scores =  PWMs_proteome_scores_Length_spec(lib_path,Alleles_out,PWMs_pred_dict_,L,alphas_)
    scores = prot_scores

    standard_dev = standard_dev + std_deviation(Alleles_out,scores,L)

    L_dist_pred = L_predict_model(Allele_pos,lib_path)
    L_dist_pred = L_dist_test(Alleles_out,L_dist_pred)
    #### --------------------
    L_dist_pred.iloc[:, 1:] = L_dist_pred.iloc[:, 1:].div(L_dist_pred.iloc[:, 1:].sum(axis=1), axis=0)
    #### --------------------
    if output_motifs==1:
        L_dist_pred.to_csv(f'{output_dir}/Motifs/PLD_Alleles.txt',sep='\t',index=False)

        for allele in Alleles_out:
            fig, ax = plt.subplots(figsize=(4, 3))
            fig.set_facecolor('white'),ax.set_facecolor('white')
            ax.plot(L, L_dist_pred[L_dist_pred['Allele']==allele].iloc[0,1:],color='#047FE5')
            ax.spines['top'].set_visible(False),ax.spines['right'].set_visible(False)
            ax.spines['left'].set_edgecolor('grey'),ax.spines['bottom'].set_edgecolor('grey')
            ax.set_ylim([0,1]),ax.set_xlabel('Peptide length')
            safe_filename = allele.replace(":", "_")
            fig.savefig(f'{output_dir}/Motifs/{allele}_PLD.png', bbox_inches='tight',dpi=300)


    bias = bias + compute_bias(Alleles_out,L_dist_pred,scores,L)
    PWMs_prot_scores_corr = prot_scores_correction(Alleles_out,scores,bias,standard_dev,L)
    Alleles_perRank_f =Alleles_perRank_f + perRank_function(PWMs_prot_scores_corr)

Alleles_to_testt = Alleles_in + Alleles_out

if (len(PWMs_pred_dict) == 0):
    
    PWMs_pred_dict = PWMs_pred_dict_
    alphas = alphas_
elif len(PWMs_pred_dict_) > 1:
    for i in range(len(L)):
        PWMs_pred_dict[i] = PWMs_pred_dict[i] + PWMs_pred_dict_[i]
        alphas[i] = alphas[i] + alphas_[i]

ligands = Ligands_scores_spec(Ligands,Ligands_L,Alleles_to_testt,PWMs_pred_dict,alphas,bias,standard_dev,Alleles_perRank_f)
columns_df = ['Peptide','Score_bestAllele','BestAllele','%Rank_bestAllele']
for x in Alleles_to_test:
    columns_df+=[f'Score_{x}',f'%Rank_{x}']
ligands = ligands[columns_df]

ligands = ligands.round(6)

closest_Alleles,distance_scores = distance_to_training(lib_path,Alleles,Alleles_to_test)
Alleles_quality_info = [f'{closest_Alleles[i]} ({np.round(distance_scores[i],4)})' for i in range(len(closest_Alleles))]
if len(Alleles_out)>0:
    header_comments = [
        "####################",
        "# Output from MixMHCpred (v3.0)",
        f"# Alleles: {', '.join(Alleles_to_test)} - predicted motif for {', '.join(Alleles_out)}",
        f"# Closest Allele (distance score): {' -- '.join(Alleles_quality_info)}",
        f"# Input file: {file_input}",
        "# MixMHCpred is freely available for academic users.",
        "# Private companies should contact Nadette Bulgin (nbulgin@lcr.org) at the Ludwig Institute for Cancer Research Ltd for commercial licenses.",
        "# ",
        "# To cite MixMHCpred3.0, please refer to:",
        "# Tadros et al., Predicting MHC-I ligands across alleles and species: How far can we go?, BioRxiv (2024).",
        "####################"
    ]
else:
    header_comments = [
        "####################",
        "# Output from MixMHCpred (v3.0)",
        f"# Alleles: {', '.join(Alleles_to_test)}",
        f"# Closest Allele (distance score): {' -- '.join(Alleles_quality_info)}",
        f"# Input file: {file_input}",
        "# MixMHCpred is freely available for academic users.",
        "# Private companies should contact Nadette Bulgin (nbulgin@lcr.org) at the Ludwig Institute for Cancer Research Ltd for commercial licenses.",
        "# ",
        "# To cite MixMHCpred3.0, please refer to:",
        "# Tadros et al., Predicting MHC-I ligands across alleles and species: How far can we go?, BioRxiv (2024).",
        "####################"
    ]



if output_motifs==1:
    with open(f'{output_dir}/Binding_predictions.txt', 'w') as f:
        f.write('\n'.join(header_comments))
        f.write('\n')
    print(f'''{output_dir}/Binding_predictions.txt is created''' )




    with open(f'{output_dir}/Binding_predictions.txt', 'a') as f:
        ligands.to_csv(f, sep='\t', index=False)

    for allele in Alleles_to_test:
        if allele in Alleles_in:
            for l in L:
                source_path = f'{lib_path}/pwm/class1_{l}/Motifs/{allele}_PWM{l}.png'
                destination_path = f'{output_dir}/Motifs/{allele}_PWM{l}.png'
                shutil.copy(source_path, destination_path)

            source_path = f'{lib_path}/pwm/PLDs/{allele}_PLD.png'
            destination_path = f'{output_dir}/Motifs/{allele}_PLD.png'
            shutil.copy(source_path, destination_path)

    # Path where your images are stored (change this to your actual path)
    image_path = f'./Motifs'
    # Image display size settings
    image_width = 200
    image_height = 225
   # Start of the HTML file with CSS for spacing
    html_content = """
    <html>
    <head>
        <title>Allele Binding Motifs and Distributions</title>
        <style>
            img {
                margin: 15px;  /* Space around images */
            }
        </style>
    </head>
    <body>
        <h1>Allele Binding Motifs and Distributions</h1>
    """

    # Adding images for each allele
    for allele in Alleles_to_test:
        html_content += f"<h2>{allele}</h2>\n"
        # Add images for binding motifs (assuming names are in a specific format)
        for length in range(8, 15):  # Lengths from 8 to 14
            image_filename = os.path.join(image_path, f"{allele}_PWM{length}.png")
            html_content += f"<img src='{image_filename}' alt='Binding Motif Length {length}' title='Binding Motif Length {length}'  width='{image_width}' height='{image_height}'/>\n"

        # Add image for peptide length distribution
        distribution_image_filename = os.path.join(image_path, f"{allele}_PLD.png")
        html_content += f"<img src='{distribution_image_filename}' alt='Peptide Length Distribution' title='Peptide Length Distribution'  width='{image_width}' height='{image_height}'/>\n"

    # End of the HTML file
    html_content += """
    </body>
    </html>
    """

    # Saving the HTML content to a file
    html_file_path = f'{output_dir}/Data_overview.html'
    with open(html_file_path, "w") as file:
        file.write(html_content)

    print(f'''
    HTML file created at {html_file_path}
    ''')
else:

    with open(file_output, 'w') as f:
        # Write the header comments first
        f.write('\n'.join(header_comments) + '\n')

        # Then write the CSV data from 'ligands'
        ligands.to_csv(f, sep='\t', index=False)

    print(f'''
    {file_output} is created
    ''' )

print('''
DONE
''')
