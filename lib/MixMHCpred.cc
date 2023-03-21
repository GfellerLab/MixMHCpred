/***************
 *Written by David Gfeller

 *The product is provided free of charge for academic users and, therefore, on an "as is" basis, without warranty of any kind.

 *For any question or commercial use, please ask david.gfeller@unil.ch

 *Copyright (2022) David Gfeller
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

void load_Pval();

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
double *****pwm;
double ***w;

int nh;
char **alleles_map;
int **Nmotifs;
double **shifts;
double **shifts_sd;
int **nc;
double **rank_val;
double **rank_thr;
int Nthr;
double **min_rank;

int Cmax;
int Lmax; //Maximum length of peptides + 1
int Lmin;

int main(int argc, char ** argv){

    int inc=6;
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
    sprintf(input_file, "%s/../temp/%d/input.txt", lib_dir, rd);

    nh=atoi(argv[5]);

    alleles=new char*[nh];
    for(int i=0; i<nh; i++){
	alleles[i]=new char[4096];
	strcpy(alleles[i], argv[inc+i]);
    }

    alleles_map=new char*[nh];
    for(int i=0; i<nh; i++){
	alleles_map[i]=new char[4096];
	strcpy(alleles_map[i], argv[inc+nh+i]);
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

    shifts_sd=new double*[Lmax];
    for(int l=Lmin; l<Lmax; l++){
	shifts_sd[l]=new double[nh];
	for(int i=0; i<nh; i++){
	  shifts_sd[l][i]=atof(argv[inc+(l-Lmin+2+2*(Lmax-Lmin))*nh+i]);
	}
    }

    //Read the number of thresholds
    char file [4096];
    string nada="";
    sprintf(file, "%s/PerRank/A0201.txt",  lib_dir);
    string line;
    ifstream myfile (file);
    if (myfile.is_open()) {
	Nthr=0;
	while (! myfile.eof() ) {
	    getline (myfile,line);
	    if(line.compare(nada) != 0) {
		Nthr++;
	    }
	}
	myfile.close();
    } else {
	cout << "Unable to open file: " << file << endl;
	exit(2);
    }

    rank_val=new double*[nh];
    for(int i=0; i<nh; i++){
	rank_val[i]=new double[Nthr];
    }

    rank_thr=new double*[nh];
    for(int i=0; i<nh; i++){
	rank_thr[i]=new double[Nthr];
    }

    load_peptides();

    load_pwm();

    load_Pval();
    //comp_Pval();

    make_pred();

    return(0);

}

//Load the %rank corresponding to different scores.
void load_Pval(){

    fstream afile;

    char *rfile;
    rfile=new char[4096];

    for(int i=0; i<nh; i++){
	sprintf(rfile, "%s/PerRank/%s.txt", lib_dir, alleles_map[i]);
	afile.open(rfile, ios::in);
	for(int j=0; j<Nthr; j++){
	    afile>>rank_val[i][j];
	    rank_val[i][j]=log(rank_val[i][j]);
	    afile>>rank_thr[i][j];
	}
	afile.close();
    }
}



void load_peptides(){

    N=20;
    Cmax=10;

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

	    for(int c=0; c<Nmotifs[l][h]; c++){
		sprintf(pwm_file, "%s/pwm/class1_%d/%s_%d.txt", lib_dir, l, alleles_map[h], c+1);
		myfile.open(pwm_file);
		for(int i=0; i<5; i++){
		    getline (myfile,line);
		}
		myfile >> w[l][h][c];
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
    double **rank;
    rank=new double*[nh];
    for(int h=0; h<nh; h++){
	rank[h]=new double[Np];
    }

    int best_pos;
    double best_rank;
    int cond;

    double tscore;
    double min_score=-100;
    double ttscore;

    //Compute the scores
    int aa;
    for(int i=0; i<Np; i++){

	for(int h=0; h<nh; h++){
	    score[h][i]=0;


	    for(int c=0; c<Nmotifs[lg[i]][h]; c++){
		tscore=1;

		for(int p=0; p<lg[i]; p++){
		    aa=peptides[i][p];
		    tscore=tscore*(pwm[lg[i]][h][c][aa][p]);
		}
		score[h][i]=score[h][i]+w[lg[i]][h][c]*tscore;
	    }
	    score[h][i]=log(score[h][i])/lg[i];
	    score[h][i]=(score[h][i]-shifts[lg[i]][h])/shifts_sd[lg[i]][h];

	    //Now compute the ranks
	    if(score[h][i]>=rank_thr[h][0]){
		rank[h][i]=exp(rank_val[h][0]);
	    } else if(score[h][i]<=rank_thr[h][Nthr-1]){
		rank[h][i]=100;
	    } else  {
		//Do the extrapolation between bins of ranks
		cond=1;
		for(int j=1; j<Nthr && cond==1; j++){
		    if(score[h][i]>=rank_thr[h][j]){
			rank[h][i]=exp(rank_val[h][j-1]+(rank_val[h][j]-rank_val[h][j-1])/(rank_thr[h][j]-rank_thr[h][j-1])*(score[h][i]-rank_thr[h][j-1]));
			cond=0;
		    }
		}
	    }
	}
    }


    //Print the output

    FILE *pFile;
    pFile=fopen(output_file,"w");


    fprintf (pFile, "####################\n");
    fprintf (pFile, "# Output from MixMHCpred (v2.2)\n");
    fprintf (pFile, "# Alleles: %s",alleles[0]); for(int h=1; h<nh; h++){fprintf (pFile, ", %s", alleles[h]);} fprintf (pFile, "\n");
    fprintf (pFile, "# Input file: %s\n", input_file_original);
    fprintf (pFile, "# MixMHCpred is freely available for academic users.\n");
    fprintf (pFile, "# Private companies should contact eauffarth@licr.org or lfoit@licr.org at the Ludwig Institute for Cancer Research Ltd for commercial licenses.\n");
    fprintf (pFile, "#\n# To cite MixMHCpred2.2, please refer to:\n");
    fprintf (pFile, "# Gfeller et al. Improved predictions of antigen presentation and TCR recognition with MixMHCpred2.2 and PRIME2.0 reveal potent SARS-CoV-2 CD8+ T-cell epitopes , Cell Systems (2023).\n");
    fprintf (pFile, "# \n");
    fprintf (pFile, "####################\n");

    fprintf (pFile, "Peptide\t");
    fprintf (pFile, "Score_bestAllele\tBestAllele\t%%Rank_bestAllele");
    for(int h=0; h<nh; h++){
 	fprintf (pFile, "\tScore_%s\t%%Rank_%s", alleles[h], alleles[h]);
    }
    fprintf (pFile, "\n");

    double pval;

    for(int i=0; i<Np; i++){
	for(int p=0; p<lg[i]; p++){
	    fprintf (pFile, "%c", letter[peptides[i][p]]);
	}

	//Find the allele with the best %rank
	best_pos=-1;
	best_rank=101;
	for(int h=0; h<nh; h++){
	    if(best_rank>rank[h][i] && score[h][i] > min_score){
		best_rank=rank[h][i];
		best_pos=h;
	    }
	}

	if(best_pos > -1){
	    fprintf (pFile, "\t%.6f\t", score[best_pos][i]);
	    fprintf (pFile, "%s", alleles[best_pos]);
	    fprintf (pFile, "\t%g", rank[best_pos][i]);
	    for(int h=0; h<nh; h++){
		fprintf (pFile, "\t%.6f", score[h][i]);
		fprintf (pFile, "\t%g", rank[h][i]);
	    }
	} else {
	    fprintf (pFile, "\tNA\tNA\tNA\t");
	    for(int h=0; h<nh; h++){
		fprintf (pFile, "\tNA\tNA");
	    }
	}

	fprintf (pFile, "\n");
    }
    fclose (pFile);
}
