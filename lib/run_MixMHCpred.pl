
################
# Written by David Gfeller.

# The product is provided free of charge for academic users and, therefore, on an "as is" basis, without warranty of any kind.

# For any question or commercial use, please ask david.gfeller@unil.ch.

# Copyright (2016) David Gfeller.
################

use strict;
use Getopt::Long;

###########
# Run the MixMHCpred script
# Check input
###########


my ($alleles, $input, $dir, $output_file, $lib_dir);

GetOptions ("alleles=s" => \$alleles,    # allele
	    "input=s"   => \$input,      # input file
	    "dir=s"   => \$dir,      # current path
	    "output=s"   => \$output_file,      # output file
	    "lib=s"   => \$lib_dir)      # MixMHCpred package
    or die("Error in command line arguments\n");


##########
# Check allele name
##########


my @allele_input=split(",", $alleles);
my @allele_list=();
my $h;
my $nh=scalar @allele_input;


foreach $h (@allele_input){
    
    if(substr($h, 0, 4) eq "HLA-"){
	$h=substr($h, 4, length($h)-4);
    }
    if(substr($h, 1, 1) eq "*"){
	$h=substr($h, 0, 1).substr($h, 2, length($h)-2);
    }
    if(substr($h, 3, 1) eq ":"){
	$h=substr($h, 0, 3).substr($h, 4, length($h)-3);
    }
    push @allele_list, $h;
}


###########
# Check input file
###########

my @peptide=();
my @peptide_num=([]);
my @seq=();
my $l;
my $p;
my $s;

open IN, "$input", or die;
while($l=<IN>){
    $l =~ s/\r?\n$//;
    chomp($l);
    if(substr($l, 0, 1) ne ">" && $l ne ""){
	push @peptide, $l;
    }
}
close IN;

my %map=("A", 0, "C", 1, "D", 2, "E", 3, "F", 4, "G", 5, "H", 6, "I", 7, "K", 8, "L", 9, "M", 10, "N", 11, "P", 12, "Q", 13, "R", 14, "S", 15, "T", 16, "V", 17, "W", 18, "Y", 19);
my ($ct1, $ct2, $Lp);
my @lg=();
my @lg_pres=(); $lg_pres[9]=0; $lg_pres[10]=0;

$ct1=0;
foreach $p (@peptide){
    @seq=split('', $p);
    $ct2=0;
    foreach $s (@seq){
	if(!exists $map{$s}){
	    print "Error in input sequence - unknown amino acids $s: $p\n";
	    exit;
	} else {
	    $peptide_num[$ct1][$ct2]=$map{$s};
	}
	$ct2++;
    }
    if($ct2 != 9 && $ct2 != 10){
	print "Non-supported peptide length: $p\t$ct2 amino acids\n";
	exit;
    } else {
	push @lg, $ct2;
	$lg_pres[$ct2]=1;
    }
    $ct1++;
}

#print "9:\t$lg_pres[9]\n10:\t$lg_pres[10]\n";


#######################################
#Make sure the alleles are characterized
#######################################

my %cond_allele=();
my %Nmotifh=();
my %Ctermh=();
my %Ntermh=();
my @a=();
my @allele_list_pres9=(); for(my $i=0; $i<$nh; $i++){$allele_list_pres9[$i]=0;}
my @allele_list_pres10=(); for(my $i=0; $i<$nh; $i++){$allele_list_pres10[$i]=0;}
my @Nmotif9=(); for(my $i=0; $i<$nh; $i++){$Nmotif9[$i]=0;}
my @Nmotif10=(); for(my $i=0; $i<$nh; $i++){$Nmotif10[$i]=0;}
my @Cterm9=(); for(my $i=0; $i<$nh; $i++){$Cterm9[$i]=0;}
my @Cterm10=(); for(my $i=0; $i<$nh; $i++){$Cterm10[$i]=0;}
my @Nterm9=(); for(my $i=0; $i<$nh; $i++){$Nterm9[$i]=0;}
my @Nterm10=(); for(my $i=0; $i<$nh; $i++){$Nterm10[$i]=0;}
my $i;
my $t9;
my$t10;


if($lg_pres[9]==1){
   foreach $h (@allele_list){
	$cond_allele{$h}=0;
	$Nmotifh{$h}=0;
	$Ctermh{$h}=0;
	$Ntermh{$h}=0;
    }
    open IN, "$lib_dir/allele_list.txt", or die;
    $l=<IN>;
    while($l=<IN>){
	$l =~ s/\r?\n$//;
	chomp($l);
	@a=split("\t", $l);
	foreach $h (@allele_list){
	    if($a[1] eq $h && $a[0] eq 9){
		$cond_allele{$h}=1;
		$Nmotifh{$h}=$a[3];
		$Ctermh{$h}=$a[4];
		$Ntermh{$h}=$a[5];
	    }
	}
    }
    close IN;
    $t9=0;
    $i=0;

    foreach $h (@allele_list){
	#print "$h\n";
	if($cond_allele{$h}==0){
	    print "Uncharacterized allele for 9-mers: $h\n";
	    #exit();
	    $t9++;
	    $allele_list_pres9[$i]=0;
	    $Nmotif9[$i]=0;
	    $Cterm9[$i]=0;
	    $Nterm9[$i]=0;
	} else {
	    $allele_list_pres9[$i]=1;
	    $Nmotif9[$i]=$Nmotifh{$h};
	    $Cterm9[$i]=$Ctermh{$h};
	    $Nterm9[$i]=$Ntermh{$h};
	}
	$i++
    }
    
    if($t9==$nh){
	#print "None of your alleles are characterized in MixMHCpred for 9-mers...\n";
    }
}


if($lg_pres[10]==1){
    foreach $h (@allele_list){
	$cond_allele{$h}=0;
	$Nmotifh{$h}=0;
	$Ctermh{$h}=0;
	$Ntermh{$h}=0;
    }
    open IN, "$lib_dir/allele_list.txt", or die;
    $l=<IN>;
    while($l=<IN>){
	$l =~ s/\r?\n$//;
	chomp($l);
	@a=split("\t", $l);
	foreach $h (@allele_list){
	    if($a[1] eq $h && $a[0] eq 10){
		$cond_allele{$h}=1;
		$Nmotifh{$h}=$a[3];
		$Ctermh{$h}=$a[4];
		$Ntermh{$h}=$a[5];
	    }
	}
    }
    $t10=0;
    $i=0;
   
    foreach $h (@allele_list){
	if($cond_allele{$h}==0){
	    print "Uncharacterized allele for 10-mers: $h\n";
	    $t10++;
	    $allele_list_pres10[$i]=0;
	    $Nmotif10[$i]=0;
	    $Cterm10[$i]=0;
	    $Nterm10[$i]=0;
	} else {
	    $allele_list_pres10[$i]=1;
	    $Nmotif10[$i]=$Nmotifh{$h};
	    $Cterm10[$i]=$Ctermh{$h};
	    $Nterm10[$i]=$Ntermh{$h};
	}
	$i++
    }
       
    if($t10==$nh){
	#print "None of your alleles are characterized in MixMHCpred for 10-mers...\n";
    }
}



#Stop if only 9-mers and no predictions available
if( $t9==$nh && $lg_pres[10]==0 ){
    print "Sorry, no predictions available for your alleles...\n";
    exit();
}
#Stop if only 10-mers and no predictions available for all alleles
if($t10==$nh && $lg_pres[9]==0){
    print "Sorry, no predictions available for your alleles with 10-mers...\n";
    exit();
}
#Stop if no predictions available
if($t9==$nh && $t10==$nh){
    print "Sorry, no predictions available for your alleles...\n";
    exit();
}

#############
# Print input peptides
#############

my $rd=int(rand(1000000));
#$rd=101;

system("mkdir -p $lib_dir/../tmp/$rd");
open OUT, ">$lib_dir/../tmp/$rd/input.txt";
print OUT "$ct1\n";
for(my $i=0; $i<$ct1; $i++){
    print OUT "$lg[$i] $peptide_num[$i][0]";
    for(my $j=1; $j<$lg[$i]; $j++){
	print OUT " $peptide_num[$i][$j]";
    }
    print OUT "\n";
}
close OUT;

############
# Run the predictions
############

#die;

#print "$output_file $lib_dir $rd $input $nh @allele_list @allele_list_pres9 @allele_list_pres10 @Nmotif9 @Nmotif10 @Cterm9 @Cterm10 @Nterm9 @Nterm10\n";

system("$lib_dir/MixMHCpred.x $output_file $lib_dir $rd $input $nh @allele_list @allele_list_pres9 @allele_list_pres10 @Nmotif9 @Nmotif10 @Cterm9 @Cterm10 @Nterm9 @Nterm10");

if(-d "$lib_dir/../tmp/$rd"){
    system("rm -fr $lib_dir/../tmp/$rd/");
}
