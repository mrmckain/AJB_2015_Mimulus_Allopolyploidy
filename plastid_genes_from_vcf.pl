#!/usr/bin/perl -w
use strict;

####Split VCF Chloroplast file Ignoring MT-CP regions and only within cp genes

####USAGE: ignore_regions plastid_regions vcf_file

my %plastid_regions;

my %regions_to_keep;

my %ignore;
open my $ignore_file, "<", $ARGV[0];
while(<$ignore_file>){
	chomp;
	my @tarray = split /\s+/;

	for (my $i=$tarray[0]; $i <= $tarray[1]; $i++){
		$ignore{$i}=1;
	}
}

open my $annotation_file, "<", $ARGV[1];
ANNO: while(<$annotation_file>){
		chomp;
		if(/\~\~\~/){
			next;
		}
		else{
			my @tarray = split /\s+/;
			for (my $i = $tarray[1]; $i <= $tarray[2]; $i++){
				if(exists $ignore{$i}){
					delete $plastid_regions{$tarray[0]};
					next ANNO;
				}
				elsif($i >= 127647){
					next ANNO;
				}
				else{
					$plastid_regions{$tarray[0]}{$i}=1;
				}
			}
		}
}

for my $gene (keys %plastid_regions){
	for my $pos (sort {$a<=>$b} keys %{$plastid_regions{$gene}}){
		$regions_to_keep{$pos}=$gene;
	}
}
my %all_genes;
my @labels;
open my $vcf_file, "<", $ARGV[2];
while(<$vcf_file>){
	if(/^\#\#/){
		next;
	}
	elsif(/CHROM/){
		@labels = split /\s+/;
	}

	else{
		my @tarray = split /\s+/;
		if(exists $regions_to_keep{$tarray[1]}){
			my @bases;
			push(@bases, $tarray[3]);
			if($tarray[4] =~ /,/){
				push(@bases, split(",", $tarray[4]));
			}
			else{
				if($tarray[4] !~ /\./){
					push(@bases, $tarray[4]);
				}
			}
			
			for (my $i = 9; $i <= ((scalar @tarray) - 1); $i++){
				$tarray[$i] =~ /(.)\/(.)/;
				if($1 ne "."){
					$all_genes{$regions_to_keep{$tarray[1]}}{$labels[$i]}.=$bases[$1];
				}
			}
		}
	}
	
}

for my $geneid (keys %all_genes){
	open my $tout, ">", $geneid . "_mimulus_cp_from_vcf.fsa";
	for my $sample (keys %{$all_genes{$geneid}}){
		print $tout ">$sample\n$all_genes{$geneid}{$sample}\n";
	}
}
