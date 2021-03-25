#! /usr/bin/perl -w
use strict;
use Statistics::Distributions qw(chisqrprob);
if(@ARGV != 0){
	print "Error: perl ~.pl.\n\n";
	exit;
}

##
my $folder="Example_data";
my $p=1e-5;


## iGPS file
open IGPS,"$folder/All.igps.txt" or die $!;
	#[0]ID	[1]Position	[2]Code	Peptide	M. ID	Gene Name	Kinase ID	[7]Kinase Name	Interaction	[9]Predictor	Score	Cutoff
my %pk_list;	# list of pk
my %site;	# site -> pk list
while(<IGPS>){
	chomp;
	if($_=~/Interaction information/){
		last;
	}elsif($_=~/[#@]/ or $_ eq ''){
		next;
	}else{
		my @a=split(/\t/,$_);
		my $psite="$a[0]-$a[1]-$a[2]";	# ID-Position-Code
		my $kinase="$a[7]-$a[9]";	# Kinase Name-Predictor
		$pk_list{$kinase}="$a[7]\t$a[6];$a[9]";	# Kinase Name	Kinase ID;[9]Predictor
		if(exists $site{$psite}){
			$site{$psite}="$site{$psite};$kinase";
		}else{
			$site{$psite}=$kinase;
		}
	}
}
close IGPS;
# kinase list
my $i = keys(%pk_list);
print "Kinase: $i\n";
open OUT,"> Example_data/kinase-list.txt" or die $!;
foreach my $k(sort keys %pk_list){
	print OUT "$k\t$pk_list{$k}\n";
}
close OUT;
# ssKSRs site-pk
$i = keys(%site);
print "pSite: $i\n";
open OUT, "> Example_data/ssKSRs-site.txt" or die $!;
foreach my $k(sort keys %site){
	print OUT "$k\t$site{$k}\n";
}
close OUT;

## Files
my @file=("DMSO-1.txt", "DMSO-2.txt", "DMSO-3.txt", "DOX-1.txt", "DOX-2.txt", "DOX-3.txt");
my @name=("DMSO-1", "DMSO-2", "DMSO-3", "DOX-1", "DOX-2", "DOX-3");
#
my $head="Kinase\tPredictor";
my @pk_sum=(0) x ($#file+1);
print "@pk_sum\n";

##
foreach my $k(0..$#file){
	$head=$head."\t$name[$k]";
	#
	open IN,"$folder\\$file[$k]" or die $!;
	my %tmp;
	my $all=0;	#
	my $pre=0;	#
	while(<IN>){
		$all++;
		chomp;
		my @a=split(/\t/,$_);
		my $substrate="$a[0]-$a[1]-$a[2]";	# ID-Position-Code
		if(exists $site{$substrate}){
			$pk_sum[$k]=$pk_sum[$k]+$a[4];	# 
			$pre++;
			my @b=split(/;/,$site{$substrate});
			foreach my $q(0..$#b){
				if(exists $tmp{$b[$q]}){
					$tmp{$b[$q]}=$tmp{$b[$q]}+$a[4];
				}else{
					$tmp{$b[$q]}=$a[4];
				}
			}
		}
	}
	print "$k\t$folder\\$file[$k]\t$name[$k]\t$all\t$pre\t".($pre/$all)."\n";
	close IN;
	#
	foreach my $q(sort keys %pk_list){
		if(exists $tmp{$q}){
			$pk_list{$q}=$pk_list{$q}."\t".$tmp{$q};
			#$pk_sum[$k]=$pk_sum[$k]+$tmp{$q};
		}else{
			$pk_list{$q}=$pk_list{$q}."\t"."NA";
			#print "$q\n";
		}
	}
}
print "@pk_sum\n";

## 
my @vs=("3","4","5");
open VIEW,"> $folder/kinase-view-S.txt" or die $!;
print VIEW "Kinase\tID\tPredictor";
my @outfile;
$i=0;
foreach my $k(0..$#vs){
	my @a=split(/,/,$vs[$k]);
	foreach my $q(0..$#a){
		my $filename="$name[$k]-$name[$a[$q]].Intensity";
		print VIEW "\t$name[$k]-$name[$a[$q]]-E\t$name[$k]-$name[$a[$q]]-p";
		$outfile[$i]=$filename;
		$i++;
	}
}
print VIEW "\n";

##
# open OUT,"> Kinase.txt" or die $!;
# print OUT "$head\n";
if(-d $folder){
	unlink glob "$folder\\*.Intensity";
}else{
	mkdir $folder;
}
my %up;
my %down;
foreach my $k(sort keys %pk_list){
	# print OUT "$pk_list{$k}\n";
	my @a=split(/\t/,$pk_list{$k});
	$k=~s/-/\t/;
	$a[1]=~s/;/\t/;
	print VIEW "$a[0]\t$a[1]";
	foreach my $q(0..$#vs){
		my @b=split(/,/,$vs[$q]);
		if($a[$q+2] eq 'NA'){
			$a[$q+2]=0;
			foreach my $i(0..$#b){
				if($a[$b[$i]+2] eq 'NA'){
					#$a[$b[$i]+2]=0;
					print VIEW "\tNA\tNA";
					next;
				}else{
					my $s=&yates($a[$q+2],($pk_sum[$q]-$a[$q+2]),$a[$b[$i]+2],($pk_sum[$b[$i]]-$a[$b[$i]+2]) );
					#my $e=$a[$b[$i]+2]/$pk_sum[$b[$i]]/($a[$q+2]/$pk_sum[$q]);
					my $e='Inf';
					my $pvalue=&chisqrprob(1,$s)+0;
					my $filename="$name[$q]-$name[$b[$i]].Intensity";
					print VIEW "\t$e\t$pvalue";
					#
					if($pvalue < $p){
						my $id="$filename;$pvalue;$k";
						$up{$id}="$k\t$a[$q+2]\t".($pk_sum[$q]-$a[$q+2])."\t$a[$b[$i]+2]\t".($pk_sum[$b[$i]]-$a[$b[$i]+2])."\t$e\t$s\t$pvalue";
					}
				}
			}
		}else{
			foreach my $i(0..$#b){
				if($a[$b[$i]] eq 'NA'){
					$a[$b[$i]+2]=0;
				}
				my $s=&yates($a[$q+2],($pk_sum[$q]-$a[$q+2]),$a[$b[$i]+2],($pk_sum[$b[$i]]-$a[$b[$i]+2]) );
				my $pvalue=&chisqrprob(1,$s)+0;
				my $e=$a[$b[$i]+2]/$pk_sum[$b[$i]]/($a[$q+2]/$pk_sum[$q]);
				my $filename="$name[$q]-$name[$b[$i]].Intensity";
				print VIEW "\t$e\t$pvalue";
				#
				if($pvalue < $p){
					my $id="$filename;$pvalue;$k";
					if($e>=1){
						$up{$id}="$k\t$a[$q+2]\t".($pk_sum[$q]-$a[$q+2])."\t$a[$b[$i]+2]\t".($pk_sum[$b[$i]]-$a[$b[$i]+2])."\t$e\t$s\t$pvalue";
					}else{
						$down{$id}="$k\t$a[$q+2]\t".($pk_sum[$q]-$a[$q+2])."\t$a[$b[$i]+2]\t".($pk_sum[$b[$i]]-$a[$b[$i]+2])."\t$e\t$s\t$pvalue";
					}
				}
			}
		}
	}
	print VIEW "\n";
} 
# close OUT;
close VIEW;

## 
foreach my $k(sort by_file_p keys %up){
	my @a=split(/;/,$k);
	open OUT,">> $folder/$a[0]" or die $!;
	print OUT "$up{$k}\n";
	close OUT;
}
foreach my $k(0..$#outfile){
	open OUT,">> $folder/$outfile[$k]" or die $!;
	print OUT "\n";
	close OUT;
}
foreach my $k(sort by_file_p keys %down){
	my @a=split(/;/,$k);
	open OUT,">> $folder/$a[0]" or die $!;
	print OUT "$down{$k}\n";
	close OUT;
}



##
@file= glob "$folder/*.Intensity";
foreach my $k(0..$#file){
	print "$file[$k]\n";
	open IN,"$file[$k]" or die $!;
	$file[$k]=$file[$k].".kinase";
	open OUT,"> $file[$k]" or die $!;
	#
	while(<IN>){
		chomp;
		my @a=split(/\t/,$_);
		if($#a>=8){
			if($a[8]<$p){
				print OUT "$a[0]\t$a[8]\n";
			}
		}
	}
	close IN;
	close OUT;
}



## Sub
sub yates{
	my ($a,$b,$c,$d)=@_;
	my $n=$a+$b+$c+$d;
	#print "A:$a\t$b\t$c\t$d\t$n\n";
	my $s=0;
	if(($a+$b)*($a+$c)/$n<5 or ($a+$b)*($b+$d)/$n<5 or ($c+$d)*($a+$c)/$n<5 or($c+$d)*($b+$d)/$n<5 ){
		$s=$n*((abs($a*$d-$b*$c))-$n/2)**2/(($a+$b)*($c+$d)*($a+$c)*($b+$d));
	}else{
		$s=$n*($a*$d-$b*$c)**2/(($a+$b)*($c+$d)*($a+$c)*($b+$d));
	}
	#print "$s\n";
	return $s;
}
sub by_file_p{
	my @a=split(/;/,$a);
	my @b=split(/;/,$b);
	if($a[0] eq $b[0]){
		return $a[1] <=> $b[1];
	}else{
		return $a[0] cmp $b[0];
	}
}

