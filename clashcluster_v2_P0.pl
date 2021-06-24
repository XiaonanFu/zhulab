#!/usr/bin/perl -w 
use strict;
use Math::Random;
use List::Util qw(max);
use Math::CDF qw /ppois pchisq/;
use Getopt::Long;

# ------------------------------------------------------------------
# Usage: This code was design for identif the ago cluster form CLIPSEQ data
#        based on the FDR method; edited by Xiaonan Fu
#        email: xiaonan.fu@gmail.com
# ------------------------------------------------------------------
my $usage = <<USAGE;

 perl CLIP_cluster.pl <parameters>
  -e 	<str>	Input file of expression data
  -c	<str>	Input file of CLIP-SEQ or CLASH aligment file
  -s	<str>	the time point of sample[24hPBM,30hPE,3hPE,48hPBM,NBF]
  -t	<str>	Input file of transcript sequence file[FASTA]
  -b	<str>	Input file of transcript annotation file[tab delimited]
  

USAGE

my $fileExpression = "";
my $fileCLIP = "";
my $sam_timepoint = "";
my $fileTranscript = "";
my $fileAnnotation = "";
my %clip;

GetOptions ("e=s" => \$fileExpression,
			"c=s" => \$fileCLIP,
			"s=s" => \$sam_timepoint,
            "t=s" => \$fileTranscript,
			"b=s" => \$fileAnnotation);

if ($fileExpression eq "" or $fileCLIP eq "" or $fileTranscript eq "" or $fileAnnotation eq "" or $sam_timepoint eq "") {
    print $usage,"\n";
    exit;
}

#total uniq tags, this varis from sample to sample
my %z;
my %zgene;
#the expression of transcript
my %transExpression;

#the sequence of the transcript
my %transcript;

#the lenght of the transcript
my %translength;

#the gene description
my %geneDes;

#the genetic region
my %genePos;

#the gene types
my %geneType;

#samplelist
my %sampleList;



print "running start\t";printTime();

#loading the required files
print "loading the required files\n";
readTranscript($fileTranscript);
readAnnotation($fileAnnotation);
readExpression($fileExpression);

print "loading the CLIP files\n";
readCLIP("P01_clapseBlast.final","P01");
readCLIP("P02_clapseBlast.final","P02");
readCLIP("P03_clapseBlast.final","P03");

#readCLIP("301_clapseBlast.final","301");
#readCLIP("302_clapseBlast.final","302");
#readCLIP("303_clapseBlast.final","303");

#readCLIP("SF1_clapseBlast.final","SF1");
#readCLIP("SF2_clapseBlast.final","SF2");
#readCLIP("SF3_clapseBlast.final","SF3");

#readCLIP("24L_clapseBlast.final","24L");
#readCLIP("241_clapseBlast.final","241");
#readCLIP("242_clapseBlast.final","242");
#readCLIP("243_clapseBlast.final","243");

#readCLIP("481_clapseBlast.final","481");
#readCLIP("482_clapseBlast.final","482");
#readCLIP("483_clapseBlast.final","483");


#system("perl clusterCall.pl P01_clapseBlast.final >clusterP01.region");
#system("perl clusterCall.pl P02_clapseBlast.final >clusterP02.region");
#system("perl clusterCall.pl P03_clapseBlast.final >clusterP03.region");


#system("perl clusterCall.pl 301_clapseBlast.final >cluster301.region");
#system("perl clusterCall.pl 302_clapseBlast.final >cluster302.region");
#system("perl clusterCall.pl 303_clapseBlast.final >cluster303.region");


#system("perl clusterCall.pl SF1_clapseBlast.final >clusterSF1.region");
#system("perl clusterCall.pl SF2_clapseBlast.final >clusterSF2.region");
#system("perl clusterCall.pl SF3_clapseBlast.final >clusterSF3.region");

#system("perl clusterCall.pl 24L_clapseBlast.final >cluster24L.region");
#system("perl clusterCall.pl 241_clapseBlast.final >cluster241.region");
#system("perl clusterCall.pl 242_clapseBlast.final >cluster242.region");
#system("perl clusterCall.pl 243_clapseBlast.final >cluster243.region");


#system("perl clusterCall.pl 481_clapseBlast.final >cluster481.region");
#system("perl clusterCall.pl 482_clapseBlast.final >cluster482.region");
#system("perl clusterCall.pl 483_clapseBlast.final >cluster483.region");



my %pval;
readCluster("clusterP01.region","P01");
readCluster("clusterP02.region","P02");
readCluster("clusterP03.region","P03");

#readCluster("cluster301.region","301");
#readCluster("cluster302.region","302");
#readCluster("cluster303.region","303");

#readCluster("clusterSF1.region","SF1");
#readCluster("clusterSF2.region","SF2");
#readCluster("clusterSF3.region","SF3");

#readCluster("cluster241.region","241");
#readCluster("cluster242.region","242");
#readCluster("cluster243.region","243");
#readCluster("cluster24L.region","24L");

#readCluster("cluster481.region","481");
#readCluster("cluster482.region","482");
#readCluster("cluster483.region","483");

my %pval_combine;
my $totalk = keys %sampleList;
pvalCalculation();
print "peak calling by cubic spline analysis\n";
#peak calling by cubic spline analysis
system("perl clusterCall.pl $fileCLIP >cluster.region");
system("python peakcall.py cluster.region cluster.peak");


print "Filter peaks\n";
my $cutoff = 0.05;
readPeak("cluster.peak");
print "finished\n";printTime();




sub pvalCalculation{
	my $df = 2*$totalk;
	foreach my $key(sort keys %pval){
		foreach my $sky(sort {$a<=>$b} keys %{$pval{$key}}){
			my $sum = 0;
			foreach my $sample(keys %sampleList){
				if(exists $pval{$key}{$sky}{$sample}){
					if($pval{$key}{$sky}{$sample}==0){
						$pval{$key}{$sky}{$sample} = 0.000000001;
					}
					$sum+=log($pval{$key}{$sky}{$sample});
				} else {
					$sum+=log(0.5);
				}
			}
			$pval_combine{$key}{$sky}=1-pchisq(-2*$sum,$df);
		}
	}
}


sub readCluster{
	my $file = shift;
	my $type = shift;
	$sampleList{$type} = 1;
	open(CLUSTER,"$file") || die;
	while(my $ln=<CLUSTER>){
		chomp $ln;
		my @data = split(/\t/,$ln);
		my $keyid = $data[0];
		my $lamda = 1+$zgene{$keyid}{$type}/$translength{$keyid}*45;
		my ($pos,$tag) = split(/\|/,$data[1]);
		my @p = split(/\,/,$pos);
		my @t = split(/\,/,$tag);
		for(my $i=0;$i<=$#p;$i++){
			$pval{$keyid}{$p[$i]}{$type} = 1-ppois($t[$i],$lamda);
		}
	}
}

#####loading the cluster information
sub readPeak{
	my $file= shift;
	open(PEAK,"$file") || die;
	my %ppval;
	my %ppos;
	while(my $ln=<PEAK>){
		chomp $ln;
		my @dat = split(/\s+/,$ln);
		my $s = $dat[1] - 30;
		my $e = $dat[1] + 30;
		for(my $i=$s;$i<=$e;$i++){
			if(exists $pval_combine{$dat[0]}{$i} and $pval_combine{$dat[0]}{$i}<=$cutoff){
				$clip{$dat[0]."\t".$dat[1]} = $ln;
				if(not exists $ppval{$dat[0]."\t".$dat[1]}){
					$ppval{$dat[0]."\t".$dat[1]}=$pval_combine{$dat[0]}{$i};
					$ppos{$dat[0]."\t".$dat[1]} = $i;
				} else {
					$ppval{$dat[0]."\t".$dat[1]}.=','.$pval_combine{$dat[0]}{$i};
					$ppos{$dat[0]."\t".$dat[1]}.=','.$i;
				}
			}
		}
	}
	close(PEAK);
	open(FINAL,">finalclip.txt") || die;
	foreach my $key(keys %clip){
		print FINAL $clip{$key},"\t",$ppos{$key},"\t",$ppval{$key},"\n";
	}
	close(FINAL);
}

#read the expression of RNA-seq
sub readExpression{
	my $file = shift;
	my @sample;
	open(ARY, "$file") || die "can't open $!";
	while(my $line=<ARY>){
		my @data = split(/\t/,$line);
		if($.==1){
			@sample = @data;
		} else {
			for(my $i=1;$i<=$#sample;$i++){
				$transExpression{$data[0]}{$sample[$i]} = $data[$i];
			}
		}
	}
	close(ARY);
}


sub readCLIP{
	my $clipfile = shift;
	my $type = shift;
	open(CLIP, "$clipfile") || die "can't open $!";
	while(my $line=<CLIP>){
		my @data = split(/\t/,$line);
		$zgene{$data[1]}{$type}++;
	}
	close(CLIP);
}

#read in the transcript sequence
sub readTranscript{
	my $file = shift;
	open(TRAN, "$file") || die "can't open $!";
	my $id;
	while(my $line=<TRAN>){
		chomp $line;
		if($line=~/^\>/){
			$line=~s/\>//g;
			$id = $line;
		} else {
			$transcript{$id}=$line;
			$translength{$id} = length($line);
		}
	}
	close(TRAN);
}

#read in the transcript annotation
sub readAnnotation{
	my $file = shift;
	open(FF, "$file") || die "can't open $!";
	while(my $line=<FF>){
		chomp $line;
		my @data = split(/\t/,$line);
		$geneDes{$data[1]} = $data[6];
		$genePos{$data[1]} = join("\t",@data[3..5]);
		$geneType{$data[1]} = $data[2];
	}
	close(FF);
}


sub printTime{
    my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime();
    printf("%02d:%02d:%02d", $hour, $min, $sec);
    print "\n";
}
