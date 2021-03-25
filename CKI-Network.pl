#! /usr/bin/perl -w
use strict;
use Statistics::Distributions qw(chisqrprob);
# perl ~.pl
if(@ARGV != 0){
	print "Error: perl ~.pl.\n\n";
	exit;
}

###
my @vs=("4,7","5,8","6,9");
my $p=0.00001;

##
open SSKSR,"Example_data/ssKSRs-site.txt" or die $!;
my %ssksr;
while(<SSKSR>){
	chomp;
	my @a=split(/\t/,$_);
	$ssksr{$a[0]}=$a[1];
}
close SSKSR;

#### 
foreach my $i(0..$#vs){
	my @a=split(/,/,$vs[$i]);
	my $bg=$a[0];
	my $fg=$a[1];
	##
	my %up;
	my %down;
	my $sum_up=0;
	my $sum_down=0;
	my $file;
	open VIEW, "Example_data/Site-Raw.txt" or die $!;
	while(<VIEW>){
		if($_=~/^Protein/){
			chomp;
			my @a=split(/\t/,$_);
			$file="$a[$bg]-$a[$fg].Network";
			next;
		}else{
			chomp;
			my @a=split(/\t/,$_);
			my $p="$a[0]-$a[1]-$a[2]";
			if($a[$bg]==0 or $a[$fg]==0){
				next;
			}elsif($a[$fg]/$a[$bg]>=1){
				if(exists $ssksr{$p}){
					my $ratio=$a[$fg]/$a[$bg];
					#
					my @b=split(/;/,$ssksr{$p});
					foreach my $k(0..$#b){
						$sum_up=$sum_up+$ratio;
						if(exists $up{$b[$k]}){
							$up{$b[$k]}=$up{$b[$k]}+$ratio;
						}else{
							$up{$b[$k]}=$ratio;
						}
					}
				}
			}else{
				if(exists $ssksr{$p}){
					my $ratio=$a[$bg]/$a[$fg];
					#
					my @b=split(/;/,$ssksr{$p});
					foreach my $k(0..$#b){
						$sum_down=$sum_down+$ratio;
						if(exists $down{$b[$k]}){
							$down{$b[$k]}=$down{$b[$k]}+$ratio;
						}else{
							$down{$b[$k]}=$ratio;
						}
					}
				}
			}
		}
	}
	close VIEW;
	##
	open OUT,"> Example_data/$file" or die $!;
	foreach my $k(sort keys %up){
		if(! exists $down{$k}){
			$down{$k}=0;
			my $e="Inf";
			my $s=&yates($down{$k},($sum_down-$down{$k}),$up{$k},($sum_up-$up{$k}) );
			my $pvalue=&chisqrprob(1,$s)+0;
			print OUT "$k\t$down{$k}\t".($sum_down-$down{$k})."\t$up{$k}\t".($sum_up-$up{$k})."\t$e\t$s\t$pvalue\n";
		}else{
			my $e=$up{$k}/$sum_up/($down{$k}/$sum_down);
			my $s=&yates($down{$k},($sum_down-$down{$k}),$up{$k},($sum_up-$up{$k}) );
			my $pvalue=&chisqrprob(1,$s)+0;
			print OUT "$k\t$down{$k}\t".($sum_down-$down{$k})."\t$up{$k}\t".($sum_up-$up{$k})."\t$e\t$s\t$pvalue\n";
		}
	}
}

###
my @file= glob "Example_data/*.Network";
foreach my $k(0..$#file){
	print "$file[$k]\n";
	open IN,"$file[$k]" or die $!;
	$file[$k]=~s/.txt//;
	$file[$k]=$file[$k].".kinase";
	open OUT,"> $file[$k]" or die $!;
	#
	while(<IN>){
		chomp;
		my @a=split(/\t/,$_);
		$a[0]=(split(/-/,$a[0]))[0];
		if($#a>=7){
			if($a[7]<$p){
				print OUT "$a[0]\t$a[7]\n";
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
