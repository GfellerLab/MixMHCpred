/***************
 *Written by David Gfeller

 *The product is provided free of charge for academic users and, therefore, on an "as is" basis, without warranty of any kind.

 *For any question or commercial use, please ask david.gfeller@unil.ch

 *Copyright (2016) David Gfeller
 ****************/

#include <stdio.h>      /* printf, fgets */
#include <stdlib.h>     /* atoi */
#include <string.h>
#include <math.h> 
#include <iostream>
#include <fstream>
#include <algorithm>

using namespace std;

struct mycomparison
{
    bool operator() (double* lhs, double* rhs) {return (*lhs) > (*rhs);}
};

void load_peptides();

void load_pwm_binary();

void make_pred();

void comp_Pval();

int N;
int **peptides;
int *lg;
int Np;  // Number of peptide
char *output_file;
char *lib_dir;
char *input_file;
char *input_file_original;
char **alleles;
int rd;  //Random number for temporary output
double ***pwm9;
double ***pwm10;
double *P_val_bin;
int NP_val; // Number of P_values used.
double *P9;
double *P10;


int nh;
int *alleles_pres9;
int *alleles_pres10;



int main(int argc, char ** argv){

        
   
    output_file = new char[4096];
    strcpy(output_file, argv[1]);

    lib_dir = new char[4096];
    strcpy(lib_dir, argv[2]);

    rd=atoi(argv[3]);

    input_file_original = new char[4096];
    strcpy(input_file_original, argv[4]);
    
    input_file=new char[4096];
    sprintf(input_file, "%s/../tmp/%d/input.txt", lib_dir, rd);
    
    nh=atoi(argv[5]);
    
    alleles=new char*[nh];
    for(int i=0; i<nh; i++){
	alleles[i]=new char[4096];
	strcpy(alleles[i], argv[6+i]);
    }
    alleles_pres9=new int[nh];
    for(int i=0; i<nh; i++){
	alleles_pres9[i]=atoi(argv[6+nh+i]);
    }
    alleles_pres10=new int[nh];
    for(int i=0; i<nh; i++){
	alleles_pres10[i]=atoi(argv[6+2*nh+i]);
    }

    
    load_peptides();

    load_pwm_binary();

    comp_Pval();
    
    make_pred();
    
    return(0);
    
}

void comp_Pval(){

    fstream afile;
    
    int Nr9=100000;
    int Nr10=100000;

    char *letter;
    letter=new char[N+1];
    strcpy(letter, "ACDEFGHIKLMNPQRSTVWY");
    
    NP_val=117;
    P_val_bin=new double[NP_val];
    P9=new double[NP_val];
    P10=new double[NP_val];

    //Define the intervals
    for(int i=0; i<9; i++){
	P_val_bin[i]=0.0001+0.0001*i;
    }
    for(int i=9; i<18; i++){
	P_val_bin[i]=0.001+0.001*(i-9);
    }
    for(int i=18; i<NP_val; i++){
	P_val_bin[i]=0.01+0.01*(i-18);
    }
    
    int **rpep9;
    rpep9=new int*[Nr9];
    
    //Open the random peptides.
    char *rfile;
    rfile=new char[4096];

    sprintf(rfile, "%s/rand9.txt", lib_dir);
    afile.open(rfile, ios::in);
    for(int n=0; n<Nr9; n++){
        rpep9[n]=new int[9];
        for(int i=0; i<9; i++){
            afile>>rpep9[n][i];
        }
    }
    afile.close();
    
    int **rpep10;
    rpep10=new int*[Nr10];
    
    sprintf(rfile, "%s/rand10.txt", lib_dir);
    afile.open(rfile, ios::in);
    for(int n=0; n<Nr10; n++){
        rpep10[n]=new int[10];
        for(int i=0; i<10; i++){
            afile>>rpep10[n][i];
        }
    }
    afile.close();

    double *score9; score9=new double[Nr9];
    double *score10; score10=new double[Nr10];
    double tscore;
    int aa;
    int pos;
    
    //Compute the scores
    for(int n=0; n<Nr9; n++){
   	score9[n]=-1000;
     	for(int h=0; h<nh; h++){
            if(alleles_pres9[h]==1){
		tscore=0;
		for(int i=0; i<9; i++){
                    aa=rpep9[n][i];
                    tscore=tscore+log(pwm9[h][aa][i])/9;
		}
		if(tscore>score9[n]){
                    score9[n]=tscore;
		    pos=h;
		}
            }
        }	   
    }
        
    for(int n=0; n<Nr10; n++){
     	score10[n]=-1000;
	for(int h=0; h<nh; h++){
            if(alleles_pres10[h]==1){
		tscore=0;
		for(int i=0; i<10; i++){
                    aa=rpep10[n][i];
                    tscore=tscore+log(pwm10[h][aa][i])/10.0;
		}
		if(tscore>score10[n]){
		    score10[n]=tscore;
		}
            }
        }
	
    }
    
    //Sort the scores
    std::sort(score9, score9 + Nr9);
    std::sort(score10, score10 + Nr10);

    
    //Compute the ranks
    for(int n=0; n<NP_val; n++){
        P9[n]=score9[Nr9-1-int(Nr9*P_val_bin[n])];
        P10[n]=score10[Nr10-1-int(Nr10*P_val_bin[n])];
    }

    //for(int n=0; n<NP_val; n++){
    //	cout<<P_val_bin[n]<<"\t"<<P9[n]<<endl;
    //}
    
}


void load_peptides(){

     N=20;
    
    fstream afile;
    afile.open(input_file, ios::in);
    afile>>Np;

    peptides=new int*[Np];
    lg=new int[Np];
   
    
    for(int i=0; i<Np; i++){
	afile>>lg[i];
	peptides[i]=new int[lg[i]];
	for(int p=0; p<lg[i]; p++){
	    afile>>peptides[i][p];
	}
    }
    afile.close();
    
}



void load_pwm_binary(){

    pwm9=new double**[nh];
    for(int h=0; h<nh; h++){
	pwm9[h]=new double*[N];
	for(int n=0; n<N; n++){
	    pwm9[h][n]=new double[9];
	}
    }
    pwm10=new double**[nh];
    for(int h=0; h<nh; h++){
	pwm10[h]=new double*[N];
	for(int n=0; n<N; n++){
	    pwm10[h][n]=new double[10];
	}
    }
    
    char *pwm_file;
    pwm_file=new char[4096];
    sprintf(pwm_file, "%s/pwm/pwm.bin", lib_dir);

    
    int size;
    int tmp;
    double dtmp;
    char *al;
    al=new char[10];
    
    int N9, N10;
    double **tpwm9, **tpwm10;
    
    FILE *read_file = fopen(pwm_file, "rb");
    fread(&N9, 4, 1, read_file);
    fread(&N10, 4, 1, read_file);
    //cout<<N9<<" "<<N10<<endl;

    tpwm9=new double*[N];
    for(int i=0; i<N; i++){
	tpwm9[i]=new double[9];
    }
    tpwm10=new double*[N];
    for(int i=0; i<N; i++){
	tpwm10[i]=new double[10];
    }
    
    for(int n=0; n<N9; n++){
	fread(al, 8, 1, read_file);  // Read the allele
	fread(&tmp, sizeof(int), 1, read_file); // Read a separator
	fread(&size, sizeof(int), 1, read_file); // Read the number of positions
	fread(&dtmp, sizeof(double), 1, read_file); // Read a separator
	for(int i=0; i<N; i++){
	    fread(tpwm9[i], sizeof(double), 9, read_file);
	}
	for(int h=0; h<nh; h++){
	    if(strcmp(alleles[h], al) == 0){
		for(int i=0; i<N; i++){
		    for(int j=0; j<9; j++){
			pwm9[h][i][j]=tpwm9[i][j];
		    }
		}
	    }
	    
	}
    }
    
    for(int n=0; n<N10; n++){
	fread(al, 8, 1, read_file);
	fread(&tmp, sizeof(int), 1, read_file);
	fread(&size, sizeof(int), 1, read_file);
	fread(&dtmp, sizeof(double), 1, read_file);
	for(int i=0; i<N; i++){
	    fread(tpwm10[i], sizeof(double), 10, read_file);
	}
	for(int h=0; h<nh; h++){
	    if(strcmp(alleles[h], al) == 0){
		for(int i=0; i<N; i++){
		    for(int j=0; j<10; j++){
			pwm10[h][i][j]=tpwm10[i][j];
		    }
		}
	    }
	}
    
    }
    fclose(read_file);
    
}

void make_pred(){


    char *letter;
    letter=new char[N+1];
    strcpy(letter, "ACDEFGHIKLMNPQRSTVWY");
    
    double **score;
    score=new double*[nh];
    for(int h=0; h<nh; h++){
	score[h]=new double[Np];
    }
   
    double *max_score;
    max_score=new double[Np];

    int *max_pos;
    max_pos=new int[Np];

    //Compute the scores
    int aa;
    for(int i=0; i<Np; i++){
	max_score[i]=-10000;
	max_pos[i]=-1;
	for(int h=0; h<nh; h++){
	    if(lg[i]==9){
		if(alleles_pres9[h]==1){
		    score[h][i]=0;
		    for(int p=0; p<lg[i]; p++){
			aa=peptides[i][p];
			score[h][i]=score[h][i]+log(pwm9[h][aa][p]);
		    }
		    score[h][i]=1.0*score[h][i]/lg[i];
		    if(score[h][i]>max_score[i]){
			max_score[i]=score[h][i];
			max_pos[i]=h;
		    }
		}
	    }
	    if(lg[i]==10){
		if(alleles_pres10[h]==1){
		    score[h][i]=0;
		    for(int p=0; p<lg[i]; p++){
			aa=peptides[i][p];
			score[h][i]=score[h][i]+log(pwm10[h][aa][p]);
		    }
		    score[h][i]=1.0*score[h][i]/lg[i];
		    if(score[h][i]>max_score[i]){
			max_score[i]=score[h][i];
			max_pos[i]=h;
		    }
		}
	    }
	}
    }
     
    //Compute the ranks
    
    double *rank_all;
    rank_all=new double[Np];
    for(int i=0; i<Np; ++i)
    {
        rank_all[i] = max_score[i];
    }
    double *arrayofpointers[Np];
    for(int i=0; i<Np; ++i)
    {
        arrayofpointers[i]=rank_all + i;
    }

    std::sort(arrayofpointers, arrayofpointers + Np, mycomparison());

    double temp2;
    double temp1=*arrayofpointers[Np-1];
    *arrayofpointers[Np-1]=Np;

    for(int i=Np-2; i>=0 ; i--)
    {
	temp2=*arrayofpointers[i];
	if(*arrayofpointers[i]==temp1){
	    *arrayofpointers[i]=*arrayofpointers[i+1];
	} else {
	    temp1=*arrayofpointers[i];
	    *arrayofpointers[i] = i + 1;
	}
    }

    //Compute the rank for each allele

    double **rank;
    rank=new double*[nh];
    for(int h=0; h<nh; h++){
	
	rank[h]=new double[Np];
	for(int i=0; i<Np; ++i)
	{
	    rank[h][i] = score[h][i];
	}
	double *arrayofpointers[Np];
	for(int i=0; i<Np; ++i)
	{
	    arrayofpointers[i]=rank[h] + i;
	}
	
	std::sort(arrayofpointers, arrayofpointers + Np, mycomparison());

	double temp2;
	double temp1=*arrayofpointers[Np-1];
	*arrayofpointers[Np-1]=Np;
	
	for(int i=Np-2; i>=0 ; i--)
	{
	    temp2=*arrayofpointers[i];
	    if(*arrayofpointers[i]==temp1){
		*arrayofpointers[i]=*arrayofpointers[i+1];
	    } else {
		temp1=*arrayofpointers[i];
		*arrayofpointers[i] = i + 1;
	    }
	}
    }
    
   
    //Print the output
    
    FILE *pFile;
    pFile=fopen(output_file,"w");
    
    
    fprintf (pFile, "####################\n");
    fprintf (pFile, "# Output from MixMHCpred\n");
    fprintf (pFile, "# Alleles: %s",alleles[0]); for(int h=1; h<nh; h++){fprintf (pFile, ", %s", alleles[h]);} fprintf (pFile, "\n");
    fprintf (pFile, "# Input file: %s\n", input_file_original);
    fprintf (pFile, "# MixMHCpred is freely available for academic users.\n");
    fprintf (pFile, "# Private companies should contact eauffarth@licr.org or lfoit@licr.org at the Ludwig Institute for Cancer Research Ltd for commercial licenses.\n");
    fprintf (pFile, "#\n# To cite MixMHCpred, please refer to:\n");
    fprintf (pFile, "# Bassani-Sternberg et al. Deciphering HLA-I motifs across HLA peptidomes improves neo-antigen predictions and identifies allostery regulating HLA specificity (2017)\n");
    fprintf (pFile, "####################\n");

    fprintf (pFile, "Peptide\t");
    fprintf (pFile, "Max_score\tMax_allele\tRank\tP_val");
    for(int h=0; h<nh; h++){
	fprintf (pFile, "\t%s\tRank", alleles[h]);
    }
    fprintf (pFile, "\n");
    
    int cond;
    double pval;
    
    for(int i=0; i<Np; i++){
	for(int p=0; p<lg[i]; p++){
	    fprintf (pFile, "%c", letter[peptides[i][p]]);
	}
	
	if(max_score[i] > -10000){
	    fprintf (pFile, "\t%.6f\t",max_score[i]);
	    fprintf (pFile, "%s\t%.0f", alleles[max_pos[i]], rank_all[i]);
	    cond=0;
	    pval=1.0;
	    for(int j=0; j<NP_val && cond==0; j++){
		if(lg[i]==9){
		    if(max_score[i] > P9[j]){
			cond=1;
			pval=P_val_bin[j];
		    }
		}
		if(lg[i]==10){
		    if(max_score[i] > P10[j]){
			cond=1;
			pval=P_val_bin[j];
		    }
		}
	    }
	    fprintf (pFile, "\t%g", pval);
	} else {
	    fprintf (pFile, "\tNA\tNA\tNA\t");
	}
	for(int h=0; h<nh; h++){
	    if(lg[i]==9){
		if(alleles_pres9[h]==1){
		    fprintf (pFile, "\t%.6f\t%.0f", score[h][i], rank[h][i]);
		} else {
		    fprintf (pFile, "\tNA\tNA");
		}
	    }
	    if(lg[i]==10){
		if(alleles_pres10[h]==1){
		    fprintf (pFile, "\t%.6f\t%.0f", score[h][i], rank[h][i]);
		} else {
		    fprintf (pFile, "\tNA\tNA");
		}
	    }
	}
	fprintf (pFile, "\n");
    }
    fclose (pFile);
    
}



