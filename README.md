
# MixMHCpred3.0 Usage Guide

## Overview

MixMHCpred3.0 is a pan-Allele predictor designed for predicting peptide bindings to MHC alleles and can also perform sequence alignment and binding motif plotting. 

The product is provided free of charge, and, therefore, on an "as is" basis, without warranty of any kind.

FOR-PROFIT USERS:
If you plan to use MixMHCpred (version 3.0) or any data provided with the script in any for-profit
application, you are required to obtain a separate license.
To do so, please contact Nadette Bulgin (nbulgin@lcr.org) at the Ludwig Institute for Cancer Research Ltd.

If you use MixMHCpred3.0 in a publication, please cite:
Tadros et al. BioRxiv (2024).

For the older version (MixMHCpred2.2), please cite:
Gfeller et al. Improved predictions of antigen presentation and TCR recognition with MixMHCpred2.2 and PRIME2.0 reveal potent SARS-CoV-2 CD8+ T-cell epitopes , Cell Systems (2023).

For scientific questions, please contact David Gfeller (David.Gfeller@unil.ch)

Copyright (2024) David Gfeller


## NEW FEATURES OF VERSION 3.0

- A pan-Allele predictor (Use of neural networks)
- Expanded training set of MHC alleles of human and other species 
- Sequence Alignment
- HTML file to visualize the binding motifs and peptide length distributions (Optional)

## Prerequisites

- Python 3.x installed (available here: https://www.python.org/downloads/)
- Bash shell (default on macOS and Linux)
- For Windows users:
    - MixMHCpred is called through a bash script. So you need to have bash available on the Windows computers. If you don’t have one already, you could install git-bash for example, available here: https://gitforwindows.org/.


## Installation

1. Download the MixMHCpred3.0 repository to your local machine.
2. Ensure that Python 3.x is correctly installed and accessible from your terminal or command prompt.
3. Verify that the script `MixMHCpred` is executable:
    ```
    chmod +x MixMHCpred
    ```
    - For Mac:
    Depending on your security setup, you may have to manually allow MixMHCpred executable to run (Systems Preferences -> Security & Privacy -> General)
4. **Install Required Packages**: Run the `install_packages` executable to install the necessary Python packages and MAFFT. This step ensures that all dependencies are correctly set up for MixMHCpred3.0 to function.
    ```
    chmod +x install_packages
    ./install_packages
    ```
    or manually install the required python packages (available in `code/setup_pythonLibrary.txt`) and MAFFT alignment tool (available in https://mafft.cbrc.jp/alignment/software/)
5. **Python Virtual Environment** (Recommended): It's often recommended to create and use a Python virtual environment for running MixMHCpred3.0. This helps in managing dependencies and avoiding conflicts with other Python projects. You can create a virtual environment by running:
    ```
    python3 -m venv mixmhcpred_env
    ```
   Activate the virtual environment with:
    ```
    source mixmhcpred_env/bin/activate
    ```
    Install Required Packages
    ```
    chmod +x install_packages
    ./install_packages
    ```
6. (Optional) To run MixMHCpred from anywhere on your computer, add it in your path or simply run:
    ```
    chmod +x setup_path
    ./setup_path
    ```
    in MixMHCpred3.0 directory

## Usage

The script is invoked from the command line with various parameters that control its operation. Here is the basic syntax:

```
./MixMHCpred --help
```
 

### Options

- `-i, --input`: Absolute or relative path to the input file (FASTA format or list of peptides).
- `-o, --output`: Name of the output file or directory.
- `-a, --alleles`: List of MHC alleles separated by commas (e.g., `A0301,A2402,B1501,B3906,C0303,C0702`).
- `-s, --Sequence`: Name of the FASTA file of sequences to align.
- `-p, --peptides_scoring`: Enable (1) or disable (0) binding predictions. Default is 1 for peptides scoring, 0 for sequence alignment.
- `-m, --output_motifs`: Enable (1) or disable (0) plotting of logos and creation of an HTML file for motifs. Default is 0.

### Important Notes

- The script requires either `-a/--alleles` or `-s/--Sequence` to be specified.
- If `-s/--Sequence` is used, `-a/--alleles` should not be included, and vice versa.
- Ensure there are no spaces in the path to the MixMHCpred directory, as they are not supported.

## Examples

Predicting peptide binding to specific MHC alleles and plotting motifs:

```
./MixMHCpred -i input/test.txt -o output/test_out.txt -a A0101,A0301
```
```
./MixMHCpred -i input/test.txt -o output/test_out -a A0101,A0301 -m 1 
```

Aligning sequences from a FASTA file, plotting motifs and predicting peptide binding:

```
./MixMHCpred -s input/To_align_sequences.fasta -o output/sequences_predictions 
```
```
./MixMHCpred -s input/To_align_sequences.fasta -o output/sequences_predictions -m 1
```
```
./MixMHCpred -s input/To_align_sequences.fasta -o output/sequences_predictions -m 1 -i input/test.txt -p 1
```


## Troubleshooting

- Ensure all file paths are correctly specified and exist.
- Make sure there are no spaces in the paths.

## Support

For issues, questions, or contributions, please refer to the project's GitHub repository or contact the maintainer at [david.gfeller@unil.ch].
