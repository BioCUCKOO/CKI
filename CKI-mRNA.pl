#! /usr/bin/perl -w
use strict;
# perl ~.pl
if(@ARGV != 0){
	print "Error: perl ~.pl.\n\n";
	exit;
}

##
my $p=0.01;
my @file= glob "Example_data/mRNA-*.txt";

##
my %view;
foreach my $k(0..$#file){
	open IN,"$file[$k]" or die $!;
	my %uniq;
	while(<IN>){
		chomp;
		my @a=split(/\t/,$_);
		if(exists $uniq{$a[0]}){
			next;
		}else{
			$view{$a[0]}=$a[1];
			$uniq{$a[0]}=1;
		}
	}
	close IN;
}

##
open VIEW,"> View-RNA.view" or die $!;
print VIEW "ID\tName";
foreach my $k(0..$#file){
	print "$file[$k]\n";
	open IN,"$file[$k]" or die $!;
	$file[$k]=~s/.txt//;
	print VIEW "\t$file[$k]-log2(fold_change)\t$file[$k]-p_value";
	my $file=$file[$k].".mRNA.kinase";
	open OUT,"> $file" or die $!;
	#
	my %uniq;
	while(<IN>){
		chomp;
		my @a=split(/\t/,$_);
		if(exists $uniq{$a[0]}){
			next;
		}else{
			$view{$a[0]}=$view{$a[0]}."\t$a[1]\t$a[3]";
			$uniq{$a[0]}=1;
			#
			if($a[1]=~/,/){
				$a[1]=(split(/,/,$a[1]))[0];
			}
			if($a[3]<$p){
				print OUT "$a[0]\t$a[3]\n";
			}
		}
	}
	close IN;
	#
	foreach my $k(sort keys %view){
		if(exists $uniq{$k}){
			next;
		}else{
			$view{$k}=$view{$k}."\tNA\tNA";
		}
	}
}
print VIEW "\n";

##
foreach my $k(sort keys %view){
	my $view=$k;
	print VIEW "$view\t$view{$k}\n";
}