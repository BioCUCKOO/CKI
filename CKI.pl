#! /usr/bin/perl -w
use strict;
# perl ~.pl
if(@ARGV != 0){
	print "Error: perl ~.pl\n\n";
	exit;
}

my $folder="Example_data";
my @file=glob "$folder/*.kinase";

## Kinase Name
open NAME,"$folder/Kinase-Name-R.txt" or die $!;
my %name;
while(<NAME>){
	chomp;
	my @a=split(/\t/,$_);
	if($a[2] ne ''){
		$name{$a[0]}=$a[2];
	}
}
close NAME;

## Kinase family
open FAM,"$folder/kinase-view-S.txt" or die $!;
my %family;
while(<FAM>){
	chomp;
	my @a=split(/\t/,$_);
	$family{$a[0]}=$a[2];
}
close FAM;

##
my %kinase;
my %new;
foreach my $k(0..$#file){
	open IN,"$file[$k]" or die $!;
	my %uniq;
	while(<IN>){
		my @a=split(/\t/,$_);
		#
		
		if(exists $name{$a[0]}){
			if(exists $family{$a[0]}){
				$new{$name{$a[0]}}=$family{$a[0]};
			}
			$a[0]=$name{$a[0]};
		}
		if(exists $uniq{$a[0]}){
			next;
		}else{
			$kinase{$a[0]}++;
			$uniq{$a[0]}=1;
		}
	}
	close IN;
}

##
open OUT,"> $folder/CKI-kinase.txt" or die $!;
print OUT "Name\tFamily\tNumber";
foreach my $k(0..$#file){
	print "$file[$k]\n";
	open IN,"$file[$k]" or die $!;
	$file[$k]=~s/.kinase//;
	$file[$k]=~s/$folder\///;
	print OUT "\t$file[$k]";
	#
	my %uniq;
	while(<IN>){
		chomp;
		my @a=split(/\t/,$_);
		if(exists $name{$a[0]}){
			$a[0]=$name{$a[0]};
		}
		#
		if(exists $uniq{$a[0]}){
			next;
		}else{
			$kinase{$a[0]}=$kinase{$a[0]}."\t$a[1]";
			$uniq{$a[0]}=1;
		}
	}
	close IN;
	#
	foreach my $k(sort keys %kinase){
		if(exists $uniq{$k}){
			next;
		}else{
			$kinase{$k}=$kinase{$k}."\t";
		}
	}
}
print OUT "\n";

##
foreach my $k(sort keys %kinase){
	if(exists $new{$k}){
		print OUT "$k\t$new{$k}\t$kinase{$k}\n";
	}else{
		print OUT "$k\t\t$kinase{$k}\n";
	}
	
}