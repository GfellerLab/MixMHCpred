/***************
 *Written by David Gfeller

 *The product is provided free of charge for academic users and, therefore, on an "as is" basis, without warranty of any kind.

 *For any question or commercial use, please ask david.gfeller@unil.ch

 *Copyright (2018) David Gfeller
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

void load_pwm();

void make_pred();

void comp_Pval();

int N;
int **peptides;
int *lg;
int Np;  // Number of peptide
int cys;
char *output_file;
char *lib_dir;
char *input_file;
char *input_file_original;
char **alleles;
int rd;  //Random number for temporary output
double *****pwm;
double ***w;
double *P_val_bin;
int NP_val; // Number of P_values used.
double *P;
double **P_allele;

int nh;
int *alleles_pres;
int **Nmotifs;
double **shifts;
int **nc;

int Cmax;
int Lmax; //Maximum length of peptides + 1
int Lmin;

int main(int argc, char ** argv){

    int inc=7;
    Lmax=15;
    Lmin=8;
    
    output_file = new char[4096];
    strcpy(output_file, argv[1]);

    lib_dir = new char[4096];
    strcpy(lib_dir, argv[2]);

    rd=atoi(argv[3]);

    input_file_original = new char[4096];
    strcpy(input_file_original, argv[4]);
    
    input_file=new char[4096];
    sprintf(input_file, "%s/../tmp/%d/input.txt", lib_dir, rd);

    cys=atoi(argv[5]);
    
    nh=atoi(argv[6]);

    //cout<<nh<<" HLA alleles"<<endl;
    if(cys==0){
	cout<<"Cystein containing peptides are disregarded\n";
    }
    
    alleles=new char*[nh];
    for(int i=0; i<nh; i++){
	alleles[i]=new char[4096];
	strcpy(alleles[i], argv[inc+i]);
    }
    alleles_pres=new int[nh];
    for(int i=0; i<nh; i++){
	alleles_pres[i]=atoi(argv[inc+nh+i]);
    }

    Nmotifs=new int*[Lmax];
    for(int l=Lmin; l<Lmax; l++){
	Nmotifs[l]=new int[nh];
	for(int i=0; i<nh; i++){
	    Nmotifs[l][i]=atoi(argv[inc+(l-Lmin+2)*nh+i]);
	}
    }

    shifts=new double*[Lmax];
    for(int l=Lmin; l<Lmax; l++){
	shifts[l]=new double[nh];
	for(int i=0; i<nh; i++){
	    shifts[l][i]=atof(argv[inc+(l-Lmin+2+Lmax-Lmin)*nh+i]);
	}
    }
    
    load_peptides();

    load_pwm();
    
    comp_Pval();
    
    make_pred();
    
    return(0);
    
}


void comp_Pval(){

    fstream afile;
    
    int Nr=100000;
       
    N=20;
    char *letter;
    letter=new char[N+1];
    strcpy(letter, "ACDEFGHIKLMNPQRSTVWY");
    
    NP_val=117;
    
    P=new double[NP_val];
    P_allele=new double*[nh];
    for(int h=0; h<nh; h++){
	P_allele[h]=new double[NP_val];
    }
    
    //Define the intervals
    P_val_bin=new double[NP_val];
    for(int i=0; i<9; i++){
	P_val_bin[i]=0.0001+0.0001*i;
    }
    for(int i=9; i<18; i++){
	P_val_bin[i]=0.001+0.001*(i-9);
    }
    for(int i=18; i<NP_val; i++){
	P_val_bin[i]=0.01+0.01*(i-18);
    }
    
    int ***rpep;
    rpep=new int**[Lmax];
    char *rfile;
    rfile=new char[4096];

    for(int l=Lmin; l<Lmax; l++){
	rpep[l]=new int*[Nr];
	//Open the random peptides.
	sprintf(rfile, "%s/rand%d.txt", lib_dir, l);
	afile.open(rfile, ios::in);
	for(int n=0; n<Nr; n++){
	    rpep[l][n]=new int[l];
	    for(int i=0; i<l; i++){
		afile>>rpep[l][n][i];
	    }
	}
	afile.close();
    }
    
    double *score;
    score=new double[(Lmax-Lmin)*Nr];

    double **score_allele;
    score_allele=new double*[nh];
    for(int i=0; i<nh; i++){
	score_allele[i]=new double[(Lmax-Lmin)*Nr];
    }
    
    double tscore;
    double ttscore;
    int aa;
    int pos;
    int t=0;
    
    //Compute the scores
    for(int l=Lmin; l<Lmax; l++){
	for(int n=0; n<Nr; n++){
	    score[t]=-1000;
	    if(cys==0){
		for(int i=0; i<l; i++){
		    if(rpep[l][n][i]==1){
			score[t]=-100;
		    }
		}
	    }
	    if(score[t]==-1000){
		for(int h=0; h<nh; h++){
		    if(alleles_pres[h]==1){
			tscore=0;
			for(int c=0; c<Nmotifs[l][h]; c++){
			    ttscore=1;
			    for(int i=0; i<l; i++){
				aa=rpep[l][n][i];
				ttscore=ttscore*pwm[l][h][c][aa][i];
			    }
			    tscore=tscore+w[l][h][c]*ttscore;
			}
			
			tscore=log(tscore)/(l*1.0);
			tscore=tscore-shifts[l][h];

			score_allele[h][t]=tscore;
			
			if(tscore>score[t]){
			    score[t]=tscore;
			    pos=h;
			}
		    }
		}
	    }
	    t++;
	}
    }
    //Sort the scores
    std::sort(score, score + Nr*(Lmax-Lmin));
    for(int i=0; i<nh; i++){
	std::sort(score_allele[i], score_allele[i] + Nr*(Lmax-Lmin));
    }
    
    //Compute the ranks
    for(int n=0; n<NP_val; n++){
	P[n]=score[Nr*(Lmax-Lmin)-1-int(Nr*(Lmax-Lmin)*P_val_bin[n])];
    }

    for(int i=0; i<nh; i++){
	for(int n=0; n<NP_val; n++){
	    P_allele[i][n]=score_allele[i][Nr*(Lmax-Lmin)-1-int(Nr*(Lmax-Lmin)*P_val_bin[n])];
	}
    }
        
}


void load_peptides(){

    N=20;
    Cmax=10;
    
    fstream afile;
    afile.open(input_file, ios::in);
    afile>>Np;

    //cout<<"A "<<Np<<endl;
    
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

void load_pwm(){

    pwm=new double****[Lmax];
    for(int l=Lmin; l<Lmax; l++){
	pwm[l]=new double***[nh];
	for(int h=0; h<nh; h++){
	    pwm[l][h]=new double**[Cmax];
	    for(int c=0; c<Cmax; c++){
		pwm[l][h][c]=new double*[N];
		for(int n=0; n<N; n++){
		    pwm[l][h][c][n]=new double[l];
		}
	    }
	}

    }
    w=new double**[Lmax];
    for(int l=Lmin; l<Lmax; l++){
	w[l]=new double*[nh];
	for(int h=0; h<nh; h++){
	    w[l][h]=new double[Cmax];
	}
    }

    char *pwm_file;
    ifstream myfile;
    string line;
    int tmp;
    double t;
    
    pwm_file=new char[4096];
   
    for(int l=Lmin; l<Lmax; l++){
	
    
	for(int h=0; h<nh; h++){
	    if(alleles_pres[h]==1){
		//Get the PWMs for the different motifs
		
		
		for(int c=0; c<Nmotifs[l][h]; c++){
		    sprintf(pwm_file, "%s/pwm/class1_%d/%s_%d.txt", lib_dir, l, alleles[h], c+1);
		    myfile.open(pwm_file);
		    for(int i=0; i<5; i++){
			getline (myfile,line);
		    }
		    myfile >> w[l][h][c];
		    myfile >> tmp;
		    myfile >> tmp;
		    for(int i=0; i<N; i++){
			myfile >> line;
			for(int j=0; j<l; j++){
			    myfile >> pwm[l][h][c][i][j];
			}
		    }
		    myfile.close();
		}
	    }
	    else{
		Nmotifs[l][h]=0;
	    }
	}
    
	for(int h=0; h<nh; h++){
	    t=0;
	    for(int c=0; c<Nmotifs[l][h]; c++){
		t=t+w[l][h][c];
	    }
	    for(int c=0; c<Nmotifs[l][h]; c++){
		w[l][h][c]=w[l][h][c]/t;
	    }
	}
    }
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
    double tscore;
    
    double *max_score;
    max_score=new double[Np];

    int *max_pos;
    max_pos=new int[Np];

    double ttscore;
    
    //Compute the scores
    int aa;
    for(int i=0; i<Np; i++){
	
	max_score[i]=-1000;
	max_pos[i]=-1;

	if(cys==0){
	    for(int p=0; p<lg[i]; p++){
		if(peptides[i][p]==1){
		    max_score[i]=-100;
		    max_pos[i]=0;
		    for(int h=0; h<nh; h++){
			score[h][i]=-100;
		    }
		}
	    }
	}
	
	if(max_score[i]==-1000){  // Either Cys are not excluded or the peptide does not contain a Cys
	    for(int h=0; h<nh; h++){
		if(alleles_pres[h]==1){
		    score[h][i]=0;
		    //cout<<Nmotifs[lg[i]][h]<<endl;
		    for(int c=0; c<Nmotifs[lg[i]][h]; c++){
			tscore=1;
			for(int p=0; p<lg[i]; p++){
			    aa=peptides[i][p];
			    tscore=tscore*(pwm[lg[i]][h][c][aa][p]);
			}
			score[h][i]=score[h][i]+w[lg[i]][h][c]*tscore;
		    }
		    //cout<<score[h][i]<<endl;
		    score[h][i]=log(score[h][i])/lg[i];
		    score[h][i]=score[h][i]-shifts[lg[i]][h];

		    if(score[h][i]>max_score[i]){
			max_score[i]=score[h][i];
			max_pos[i]=h;
		    }
		}
	    }
	}
    }
     
    //Compute the ranks
    
    /* double *rank_all;
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
    
    */
    //Print the output
    
    FILE *pFile;
    pFile=fopen(output_file,"w");
    
    
    fprintf (pFile, "####################\n");
    fprintf (pFile, "# Output from MixMHCpred (v2.0.2)\n");
    fprintf (pFile, "# Alleles: %s",alleles[0]); for(int h=1; h<nh; h++){fprintf (pFile, ", %s", alleles[h]);} fprintf (pFile, "\n");
    fprintf (pFile, "# Input file: %s\n", input_file_original);
    fprintf (pFile, "# MixMHCpred is freely available for academic users.\n");
    fprintf (pFile, "# Private companies should contact eauffarth@licr.org or lfoit@licr.org at the Ludwig Institute for Cancer Research Ltd for commercial licenses.\n");
    fprintf (pFile, "#\n# To cite MixMHCpred2.0.2, please refer to:\n");
    fprintf (pFile, "# Bassani-Sternberg et al. Deciphering HLA-I motifs across HLA peptidomes improves neo-antigen predictions and identifies allostery regulating HLA specificity, PLoS Comp Bio (2017).\n");
    fprintf (pFile, "# Gfeller et al. The length distribution and multiple specificity of naturally presented HLA-I ligands, J Immunol (2018).\n");
    fprintf (pFile, "####################\n");

    fprintf (pFile, "Peptide\t");
    //fprintf (pFile, "Max_score\tMax_allele\tRanking\tP_val");
    fprintf (pFile, "Score_bestAllele\tBestAllele\t%%Rank_bestAllele");
    for(int h=0; h<nh; h++){
	//fprintf (pFile, "\t%s\tRanking\tP_val", alleles[h]);
 	fprintf (pFile, "\tScore_%s\t%%Rank_%s", alleles[h], alleles[h]);
    }
    fprintf (pFile, "\n");
    
    int cond;
    double pval;

    for(int i=0; i<Np; i++){
	for(int p=0; p<lg[i]; p++){
	    fprintf (pFile, "%c", letter[peptides[i][p]]);
	}
	
	if(max_score[i] > -10000){
	    fprintf (pFile, "\t%.6f\t", max_score[i]);
	    //fprintf (pFile, "%s\t%.0f", alleles[max_pos[i]], rank_all[i]);
	    fprintf (pFile, "%s", alleles[max_pos[i]]);
	    cond=0;
	    pval=1.0;
	    for(int j=0; j<NP_val && cond==0; j++){
		if(max_score[i] > P[j]){
		    cond=1;
		    pval=P_val_bin[j];
		}
	    }
	    fprintf (pFile, "\t%g", 100*pval);
	} else {
	    fprintf (pFile, "\tNA\tNA\tNA\t");
	}
	for(int h=0; h<nh; h++){
	    if(alleles_pres[h]==1){
		//fprintf (pFile, "\t%.6f\t%.0f", score[h][i], rank[h][i]);
		fprintf (pFile, "\t%.6f", score[h][i]);
		cond=0;
		pval=1.0;
		for(int j=0; j<NP_val && cond==0; j++){
		    if(score[h][i] > P_allele[h][j]){
			cond=1;
			pval=P_val_bin[j];
		    }
		}
		fprintf (pFile, "\t%g", 100*pval);
	    
	    } else {
		fprintf (pFile, "\tNA\tNA");
	    }
	}
	fprintf (pFile, "\n");
    }
    fclose (pFile);
}
