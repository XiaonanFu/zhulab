#!/usr/bin/perl -w
use strict;
use Getopt::Long;

# ------------------------------------------------------------------
# Usage: This code was design to remove PCR induced bias; edited by Xiaonan Fu
#        email: xiaonan.fu@gmail.com
# ------------------------------------------------------------------
my $usage = <<USAGE;

 perl step1_blastClapse.pl <parameters>
  -b 	<str>	Input file of blast result
  
USAGE


my $blast="";
my $annotation="";
my %tag;
my %uniqcode;
my %biotype;

GetOptions ("b=s" => \$blast
			);

if ($blast eq "") {
    print $usage,"\n";
    exit;
}


readBlast($blast);

output();



sub readBlast{
	my $file = shift;
	open(BLAST,"$file") || die;
	while(my $ln=<BLAST>){
		chomp $ln;
		my @content = split(/\t/,$ln);

		if($content[9]>$content[10]){
			($content[9],$content[10]) = ($content[10],$content[9]);
		}
		#minimum cut
		if(($content[10]-$content[9])<20){
			next;
		}
		$tag{$content[0]}=$ln;
		my $tagcode = " ";
		if($content[0]=~/\_/){
			$tagcode = [split(/\_/,$content[0])]->[1];
		}
		if(not exists $uniqcode{$content[1].$content[9].$content[10].$tagcode}){
			$uniqcode{$content[1].$content[9].$content[10].$tagcode} = $content[0];
		} else {
			$uniqcode{$content[1].$content[9].$content[10].$tagcode} .= '|'.$content[0];
		}
	}
}


sub output{
	foreach my $key(sort keys %uniqcode){
		if($uniqcode{$key}=~/\|/){
			my @taglist = split(/\|/,$uniqcode{$key});
			print $tag{$taglist[0]},"\t",join("\|",@taglist),"\n";
		} else {
			print $tag{$uniqcode{$key}},"\n";
		}
	}
}
