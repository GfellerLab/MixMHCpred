
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


my ($alleles, $input, $dir, $output_file, $lib_dir, $cys);

GetOptions ("alleles=s" => \$alleles,    # allele
	    "input=s"   => \$input,      # input file
	    "dir=s"   => \$dir,      # current path
	    "output=s"   => \$output_file,      # output file
	    "lib=s"   => \$lib_dir,      # MixMHCpred package
 	    "cys=i"   => \$cys)      # Do not include Cystein containing peptides
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

open IN, "$input", or die "No input\n";
while($l=<IN>){
    $l =~ s/\r?\n$//;
    chomp($l);
    if(substr($l, 0, 1) ne ">" && $l ne ""){
	push @peptide, $l;
    }
}
close IN;


if(scalar @peptide == 0){
    print "Empty peptide file\n";
    exit;
}

my $Lmin=8;
my $Lmax=14;

my %map=("A", 0, "C", 1, "D", 2, "E", 3, "F", 4, "G", 5, "H", 6, "I", 7, "K", 8, "L", 9, "M", 10, "N", 11, "P", 12, "Q", 13, "R", 14, "S", 15, "T", 16, "V", 17, "W", 18, "Y", 19);
my ($ct1, $ct2, $Lp);
my @lg=();
my @lg_pres=(); # Keep track whether peptides of this lenght are found in input
for(my $le=$Lmin; $le<=$Lmax; $le++){
    $lg_pres[$le]=0;
}

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
    if($ct2 < $Lmin || $ct2 > $Lmax){
	print "Incompatible peptide length: $p\t$ct2. \nOnly peptides of length $Lmin-$Lmax are supported\n";
	exit;
    } else {
	push @lg, $ct2;
	$lg_pres[$ct2]=1;
    }
    $ct1++;
}


#######################################
# Load information about the allelels
#######################################



my @a=();
my @Nmotif=([]); for(my $le=$Lmin; $le<=$Lmax; $le++){ for(my $i=0; $i<$nh; $i++){$Nmotif[$le][$i]=0;}}
my @shifts=([]); for(my $le=$Lmin; $le<=$Lmax; $le++){ for(my $i=0; $i<$nh; $i++){$shifts[$le][$i]=0;}}
my @allele_pres=(); for(my $i=0; $i<$nh; $i++){ $allele_pres[$i]=0; }

my $t;
my $cond=0;

open IN, "$lib_dir/allele_list.txt", or die;
$l=<IN>;
while($l=<IN>){
    $l =~ s/\r?\n$//;
    chomp($l);
    @a=split("\t", $l);
    $t=0;
    foreach $h (@allele_list){
	if($a[1] eq $h){
	    $Nmotif[$a[0]][$t]=$a[2];
	    $allele_pres[$t]=1;
	    $cond=1;
	}
	$t++;
    }
}
close IN;

open IN, "$lib_dir/shifts.txt", or die;
while($l=<IN>){
    $l =~ s/\r?\n$//;
    chomp($l);
    @a=split("\t", $l);
    $t=0;
    foreach $h (@allele_list){
	if($a[1] eq $h){
	    $shifts[$a[0]][$t]=$a[2];
	    if(substr($h, 0, 1) eq "C"){
		#$shifts[$a[0]][$t]=$a[2]+0.3;
	    }
	    #$shifts[$a[0]][$t]=0;   #Set all the shifts to 0.
	}
	$t++;
    }
}

######################################
# Stop if no predictions are available
######################################


if( $cond==0 ){
    print "Predictions are not available in MixMHCpred for any of your alleles...\n";
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

my @shifts_all=();
my @Nmotif_all=();

for(my $l=$Lmin; $l<=$Lmax; $l++){
    push @Nmotif_all, @{$Nmotif[$l]};
    push @shifts_all, @{$shifts[$l]};
}

#print "$output_file $lib_dir $rd $input $cys $nh @allele_list @allele_pres @Nmotif_all @shifts_all\n";
#die;


system("$lib_dir/MixMHCpred.x $output_file $lib_dir $rd $input $cys $nh @allele_list @allele_pres @Nmotif_all @shifts_all");

if(-d "$lib_dir/../tmp/$rd"){
    system("rm -fr $lib_dir/../tmp/$rd/");
}
