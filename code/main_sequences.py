from argparse import ArgumentParser
import numpy as np
import pandas as pd 
import warnings
import subprocess
import os
import shutil

warnings.simplefilter(action='ignore', category=Warning)
warnings.filterwarnings("ignore")

parser=ArgumentParser()
parser.add_argument("-i", "--InputFilePath",help="file path of the peptidess",default=None)
parser.add_argument("-s", "--SequenceFilePath",help="file path of the sequences")
parser.add_argument("-o", "--OutputFilePath",help="file path where to save results")
parser.add_argument("-l","--lib_path",help="Library directory")
parser.add_argument("-p","--score_peptides",help="to compute binding scores",default='0')
parser.add_argument("-m","--output_motifs",help="create the binding motifs",default='0')

args=parser.parse_args()


file_sequence=args.SequenceFilePath
if os.path.isfile(file_sequence) == False:
    print(f'''
    ERROR: No such file : {file_sequence} exists.
    ''')
    exit()



output_dir = args.OutputFilePath

if os.path.isfile(output_dir):
    print(f'''
    ERROR: A file named '{output_dir}' already exists.
        
     Please choose a different folder name to save multiple files in it.
    ''') 
    exit()


lib_path = args.lib_path
score_peptides = int(args.score_peptides)

if (args.InputFilePath is not None) and (score_peptides == 0):
    print(f'''
    Please add '-p 1' to the command line to compute the binding predictions for the peptide file '{args.InputFilePath}'.
    Or remove '-i {args.InputFilePath}' from the command line 
    ''')
    exit()
if score_peptides == 1:
    
    file_input=args.InputFilePath

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

        
output_motifs = int(args.output_motifs)


# Check if the output directory exists
if os.path.exists(output_dir):
    # Remove it if it exists
    shutil.rmtree(output_dir)

# Create the output directory
os.makedirs(output_dir, exist_ok=True)

# Construct the mafft command with the correct paths
# mafft_command = f"mafft --keeplength --add {file_sequence} {lib_path}/aligned_sequences.fasta > {output_dir}/final_alignment.fasta"
mafft_command = f"mafft --keeplength --add {file_sequence} {lib_path}/aligned_sequences.fasta > {output_dir}/final_alignment.fasta 2>/dev/null"

# Execute the mafft command
subprocess.run(mafft_command, shell=True)

def read_fasta_identifiers(file_path):
    identifiers = []
    with open(file_path, 'r') as file:
        for line in file:
            if line.startswith('>'):
                identifiers.append(line.strip())
    return identifiers

def filter_fasta_sequences(source_file, identifiers):
    filtered_sequences = []
    write_sequence = False
    with open(source_file, 'r') as infile:
        for line in infile:
            if line.startswith('>'):
                write_sequence = line.strip() in identifiers
            if write_sequence:
                filtered_sequences.append(line)
    return filtered_sequences

# Paths to your files
sequences_to_align_file = f'{file_sequence}'
final_alignment_file = f'{output_dir}/final_alignment.fasta'

# Extract identifiers
identifiers = read_fasta_identifiers(sequences_to_align_file)

# Filter sequences
filtered_sequences = filter_fasta_sequences(final_alignment_file, identifiers)

# Overwrite the final_alignment.fasta file with filtered sequences
with open(final_alignment_file, 'w') as outfile:
    outfile.writelines(filtered_sequences)

if output_motifs==1:
    import matplotlib.pyplot as plt
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
if score_peptides == 1 or output_motifs==1:
    from MixMHCpred import *
    from panPredictor import *

    def transform_allele(Allele_Pos_idx,Alleles,Alleles_seq,blosum,dictionary):
        Allele_pos = []
        new_alleles = []
        new_sequences = []
        for i in range(len(Alleles_seq)):
            allele_sequence="".join([Alleles_seq[i][z] for z in Allele_Pos_idx])
            
            if '-' in allele_sequence:
                print(f'gap `-` appears at the binding site sequence for {Alleles[i]}. No predictions can be done for this Allele')
                
            else:   
                new_alleles.append(Alleles[i])
                new_sequences.append(Alleles_seq[i])
                blos_seq = blosum_array(allele_sequence,dictionary,blosum)
                OHE_seq = one_hot_encode(allele_sequence)
                Allele_pos.append([x + y for x, y in zip(blos_seq, OHE_seq)])

        
        return Allele_pos,new_alleles,new_sequences
 
    # Initialize lists to hold allele names and sequences
    allele_names = []
    sequences = []


    # Read the FASTA file
    with open(final_alignment_file, 'r') as file:
        sequence = ''
        for line in file:
            line = line.strip()
            if line.startswith('>'):  # Header line
                if sequence:  # Save the previous sequence
                    sequences.append(sequence)
                    sequence = ''
                allele_name = line[1:]  # Remove '>'
                allele_names.append(allele_name)
            else:
                sequence += line
        sequences.append(sequence) 

    blosum = load_blosum(bl_path=f'{lib_path}/blosum62.txt')
    dictionary=load_dictionary(f'{lib_path}/blosum62.txt')
    blosum_t = np.transpose(blosum)
    AA_Human_freq = list(pd.read_csv(f'{lib_path}/proteome/AA_frequency_human_Methionine.csv')['Frequency'])
        

        
    Allele_Pos_idx = get_Allele_pos_idx(f'{lib_path}/Allele_pos')
    Allele_pos,allele_names,allele_sequences= transform_allele(Allele_Pos_idx,allele_names,sequences,blosum,dictionary)

    if len(allele_names)== 0:
        print('Only Alignment can be done')
        exit()

    L = [8,9,10,11,12,13,14]
    alphas = []
    PWMs_pred_dict = []
    Alleles_perRank_f = []

    for l in L:
        
        alphas.append([])

        PWMs_test = predict_model(lib_path,Allele_pos,PL=l)
        PWMS_pred_bl = Blosum_Corr_pred(blosum_t,PWMs_test,AAs=AA,N=5000,PL=l)

        PWMs_pred = Freq_Corr_spec(allele_names,AA_Human_freq,PWMS_pred_bl)
        if output_motifs==1:
            print(f'Creating Motifs for length {l}')
            for i in range(len(allele_names)):
                
                os.makedirs(f'{output_dir}/Motifs', exist_ok=True)
                safe_filename = allele_names[i].replace(":", "_")
                safe_filename = safe_filename.replace("*", "_")
                PWMs_pred[i][0].to_csv(f'{output_dir}/Motifs/{safe_filename}_PWM{l}.txt',sep='\t')

                plt.ioff()

                pwm_df = PWMs_pred[i][0].T
                pwm_df = pwm_df.div(pwm_df.sum(axis=1), axis=0)
                pssm_df = logomaker.transform_matrix(pwm_df, 
                                                    from_type='probability', 
                                                    to_type='information')

                fig, ax = plt.subplots(figsize=(4, 3))
                fig.set_facecolor('white'),ax.set_facecolor('white')
                ax.spines['top'].set_visible(False),ax.spines['right'].set_visible(False)
                ax.spines['left'].set_edgecolor('grey'),ax.spines['bottom'].set_edgecolor('grey')
                logo = logomaker.Logo(pssm_df,color_scheme=chemistry_colors,ax=ax)
                logo.ax.set_ylim([0,4])
                positions = range(1, len(pwm_df) + 1)
                logo.ax.set_xticks(positions),logo.ax.set_xticklabels(positions)
                logo.ax.grid(False)
                
                fig.savefig(f'{output_dir}/Motifs/{safe_filename}_PWM{l}.png', bbox_inches='tight',dpi=300)


        alphas[-1] = [[1.0] for t in range(len(allele_names))]
        PWMs_pred_dict.append(creation_PWM_dict_spec(PWMs_pred,AAs=AA,PL=l))

    L_dist_pred = L_predict_model(Allele_pos,lib_path)
    L_dist_pred = L_dist_test(allele_names,L_dist_pred)
    #### --------------------
    L_dist_pred.iloc[:, 1:] = L_dist_pred.iloc[:, 1:].div(L_dist_pred.iloc[:, 1:].sum(axis=1), axis=0)
    #### --------------------
    if output_motifs==1:
        L_dist_pred.to_csv(f'{output_dir}/Motifs/PLD_sequences.txt',sep='\t',index=False)

        for allele in allele_names:
            fig, ax = plt.subplots(figsize=(4, 3))
            fig.set_facecolor('white'),ax.set_facecolor('white')
            ax.plot(L, L_dist_pred[L_dist_pred['Allele']==allele].iloc[0,1:],color='#047FE5')
            ax.spines['top'].set_visible(False),ax.spines['right'].set_visible(False)
            ax.spines['left'].set_edgecolor('grey'),ax.spines['bottom'].set_edgecolor('grey')
            ax.set_ylim([0,1]),ax.set_xlabel('Peptide length')
            safe_filename = allele.replace(":", "_")
            safe_filename = safe_filename.replace("*", "_")
            fig.savefig(f'{output_dir}/Motifs/{safe_filename}_PLD.png', bbox_inches='tight',dpi=300)

def distance_to_training_special(lib_path,training_Alleles,seq):  

    Allele_Pos_idx = get_Allele_index(f'{lib_path}/Allele_pos')
    # print(training_Alleles)
    seq_data = pd.read_csv(f"{lib_path}/MHC_I_sequences.txt", delim_whitespace= True)

    Alleles_seq_training = [seq_data[seq_data['Allele']== allele]['Sequence'].iloc[0] for allele in training_Alleles]

    all_alleles = seq_data['Allele'].tolist()
    Alleles_seq = seq_data['Sequence'].tolist()

    Alleles_seq_predicted = seq


    dict_load=np.load(f'{lib_path}/blosum62_update.npy', allow_pickle=True)
    blosum62=dict_load.item()    
   


    training_sequences = ["".join([allele[z] for z in Allele_Pos_idx]) for allele in Alleles_seq_training]
    sequences_ = ["".join([allele[z] for z in Allele_Pos_idx]) for allele in Alleles_seq]
    predicted_sequences = ["".join([allele[z] for z in Allele_Pos_idx]) for allele in Alleles_seq_predicted]

    close_Alleles = []
    sim_scores = []

    Alleles_database = []
    scores_database = []

    for i in range(len(predicted_sequences)):

        Alleles_sim, score = sim_sequence(training_Alleles,training_sequences,predicted_sequences[i],blosum62)
        close_Alleles.append(Alleles_sim)
        sim_scores.append(score)

        Alleles_sim, score = sim_sequence(all_alleles,sequences_,predicted_sequences[i],blosum62)
        Alleles_database.append(Alleles_sim)
        scores_database.append(score)

    return close_Alleles,sim_scores,Alleles_database,scores_database



if output_motifs==1:
    alleles_list = pd.read_csv(f'{lib_path}/alleles_list.txt',sep='\t')
    Alleles = alleles_list['Allele'].tolist()
    closest_Alleles,distance_scores,Alleles_database,scores_database = distance_to_training_special(lib_path,Alleles,allele_sequences)
    # Alleles_quality_info = [f'{closest_Alleles[i]} ({np.round(distance_scores[i],4)})' for i in range(len(closest_Alleles))]

    # Path where your images are stored (change this to your actual path)
    image_path = f'./Motifs'
    # Image display size settings
    image_width = 200 
    image_height = 200
    # Start of the HTML file
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
    for xx in range(len(allele_names)) :
        html_content += f"<h2>{allele_names[xx]}</h2>\n"
        html_content += f"<h2>closest Allele = {closest_Alleles[xx]}, distance = {np.round(distance_scores[xx],4)}</h2>\n"
        html_content += f"<h2>closest Allele from database = {Alleles_database[xx]}, distance = {np.round(scores_database[xx],4)}</h2>\n"

        filename = allele_names[xx].replace(":", "_")
        filename = filename.replace("*", "_")
        # Add images for binding motifs (assuming names are in a specific format)
        for length in range(8, 15):  # Lengths from 8 to 14
            image_filename = os.path.join(image_path, f"{filename}_PWM{length}.png")
            html_content += f"<img src='{image_filename}' alt='Binding Motif Length {length}' title='Binding Motif Length {length}'  width='{image_width}' height='{image_height}'/>\n"

        # Add image for peptide length distribution
        distribution_image_filename = os.path.join(image_path, f"{filename}_PLD.png")
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

    print(f"HTML file created at {html_file_path}")



if score_peptides == 1:
    
    # file_input=args.InputFilePath

    # # Check input file
    # if os.path.isfile(file_input) == False:
    #     print(f'''
    #     ERROR: No such file : {file_input} exists.
    #     ''')
    #     exit()
    # else:
    #     with open(file_input, 'r') as f:
    #         peptide = [line.strip() for line in f if not line.startswith(">") and line.strip()]

    # if not peptide:
    #     print("Empty peptide file")
    #     exit()

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

    print('Binding Predictions ...')
    Ligands = pd.DataFrame({'Peptide':peptide})
    del peptide
    Ligands['length'] = Ligands.Peptide.str.len()
    Ligands_L = np.unique(Ligands['length'])


    prot_scores =  PWMs_proteome_scores_Length_spec(lib_path,allele_names,PWMs_pred_dict,L,alphas)
    scores = prot_scores

    bias = []
    standard_dev = []

    standard_dev = standard_dev + std_deviation(allele_names,scores,L)


    bias = bias + compute_bias(allele_names,L_dist_pred,scores,L)
    PWMs_prot_scores_corr = prot_scores_correction(allele_names,scores,bias,standard_dev,L)
    Alleles_perRank_f =Alleles_perRank_f + perRank_function(PWMs_prot_scores_corr)

    ligands = Ligands_scores_spec(Ligands,Ligands_L,allele_names,PWMs_pred_dict,alphas,bias,standard_dev,Alleles_perRank_f)
    columns_df = ['Peptide','Score_bestAllele','BestAllele','%Rank_bestAllele']
    for x in allele_names:
        columns_df+=[f'Score_{x}',f'%Rank_{x}']
    ligands = ligands[columns_df]

    ligands = ligands.round(6)
    
    header_comments = [
            "####################",
            "# Output from MixMHCpred (v3.0)",
            f"# Predictions for Alleles : {', '.join(allele_names)}",
            f"# Input file: {file_input}",
            "# MixMHCpred is freely available for academic users.",
            "# Private companies should contact Nadette Bulgin (nbulgin@lcr.org) at the Ludwig Institute for Cancer Research Ltd for commercial licenses.",
            "#",
            "# To cite MixMHCpred3.0, please refer to:",
            "# Gfeller et al. Improved predictions of antigen presentation and TCR recognition with MixMHCpred2.2 and PRIME2.0 reveal potent SARS-CoV-2 CD8+ T-cell epitopes , Cell Systems (2023).",
            "# Tadros et al., Predicting MHC-I ligands across alleles and species: How far can we go?, BioRxiv (2024).",
            "####################"
        ]

    with open(f'{output_dir}/Binding_predictions.txt', 'w') as f:
        f.write('\n'.join(header_comments))
        f.write('\n')

    with open(f'{output_dir}/Binding_predictions.txt', 'a') as f:
        ligands.to_csv(f, sep='\t', index=False)





print('''
DONE
''')
