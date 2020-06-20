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
  -e 	<str>	Input file of Gff
  -c	<str>	Input file of CLIP-SEQ or CLASH aligment file
  -s	<str>	the time point of sample[24hPBM,30hPE,3hPE,48hPBM,NBF]
  -t	<str>	Input file of transcript sequence file[FASTA]
  -b	<str>	Input file of transcript annotation file[tab delimited]
  

USAGE

my $fileGff = "";
my $fileCLIP = "";
my $sam_timepoint = "";
my $fileTranscript = "";
my $fileAnnotation = "";
my %clip;

GetOptions ("e=s" => \$fileGff,
			"c=s" => \$fileCLIP,
			"s=s" => \$sam_timepoint,
            "t=s" => \$fileTranscript,
			"b=s" => \$fileAnnotation);

if ($fileGff eq "" or $fileCLIP eq "" or $fileTranscript eq "" or $fileAnnotation eq "" or $sam_timepoint eq "") {
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

my %transcript2gene;
my %id2Name;
my %geneStrand;
my %chr2GR;

print "running start\t";printTime();

#loading the required files
print "loading the required files\n";
readAnnotation($fileAnnotation);
readTranscript($fileTranscript);
readgff($fileGff);

print "loading the CLIP files\n";
#readCLIP("C01_clapseBlast.final","C01");
#readCLIP("C02_clapseBlast.final","C02");
#readCLIP("C03_clapseBlast.final","C03");
readCLIP("C04_clapseBlast.final","C04");

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

#system("perl clusterCall.pl C01_clapseBlast.final >clusterC01.region");
#system("perl clusterCall.pl C02_clapseBlast.final >clusterC02.region");
#system("perl clusterCall.pl C03_clapseBlast.final >clusterC03.region");
system("perl clusterCall.pl C04_clapseBlast.final >clusterC04.region");


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
#readCluster("clusterC01.region","C01");
#readCluster("clusterC02.region","C02");
#readCluster("clusterC03.region","C03");
readCluster("clusterC04.region","C04");

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
	open(PPV,">pvalue.txt") || die;
	while(my $ln=<CLUSTER>){
		chomp $ln;
		my @data = split(/\t/,$ln);
		if($data[0] eq "rDNA") {next;}
		my ($pos,$tag) = split(/\|/,$data[1]);
		my @p = split(/\,/,$pos);
		my @t = split(/\,/,$tag);
		my $retrGe = mapChrGene($data[0],$p[0],$p[$#p]);
		my ($keyid,$tranid,$biotypeid) = split(/\t/,$retrGe);
		if($retrGe eq "null" or $biotypeid eq "tRNA" or $biotypeid eq "rRNA" or $biotypeid eq "snRNA" or $biotypeid eq "snoRNA"){
			next;
		}
		my $lamda = 1+$zgene{$keyid}{$type}/$translength{$keyid}*45;
		for(my $i=0;$i<=$#p;$i++){
			$pval{$data[0]}{$p[$i]}{$type} = 1-ppois($t[$i],$lamda);
			print PPV $data[0],"\t",$p[$i],"\t",$zgene{$keyid}{$type},"\t",1-ppois($t[$i],$lamda),"\n";
		}
	}
	close(CLUSTER);
	close(PPV);
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
		my ($chr,$midpeak) = split(/\t/,$key);
		my ($keyid,$tranid,$biotypeid) = split(/\t/,mapChrGene($chr,$midpeak-30,$midpeak+30));
		print FINAL $clip{$key},"\t",$keyid,"\t",$biotypeid,"\t",$ppos{$key},"\t",$ppval{$key},"\n";
	}
	close(FINAL);
}


sub readCLIP{
	my $clipfile = shift;
	my $type = shift;
	my $totalZ=0;
	open(CLIP, "$clipfile") || die "can't open $!";
	while(my $line=<CLIP>){
		my @data = split(/\t/,$line);
		if($data[0]=~/^dme\-/){next;}
		my ($a,$b,$c) = split(/\t/,mapChrGene($data[1],$data[9],$data[10]));
		my $geneid = $a;
		$zgene{$geneid}{$type}++;
		if ($c eq "CDS" or $c eq "3UTR" or $c eq "5UTR"){
			$totalZ++;
		}
	}
	close(CLIP);
	print $totalZ,"\n";exit;
}

#read in the transcript sequence
sub readTranscript{
	my $file = shift;
	open(TRAN, "$file") || die "can't open $!";
	my $id;
	my $gene = "";
	while(my $line=<TRAN>){
		chomp $line;
		if($line=~/^\>/){
			$line=~s/\>//g;
			$id = $line;
			$gene = $transcript2gene{$id};
		} else {
			$transcript{$id}=$line;
			if(exists $translength{$gene}){
				if($translength{$gene} < length($line)){
					$translength{$gene} = length($line);
				}
			} else {
				$translength{$gene} = length($line);
			}
		}
	}
	close(TRAN);
}

sub printTime{
    my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime();
    printf("%02d:%02d:%02d", $hour, $min, $sec);
    print "\n";
}


#read in all the annotation
sub readAnnotation{
	my $file = shift;
	open(FF, "$file") || die "can't open $!";
	while(my $line=<FF>){
		chomp $line;
		my @data = split(/\t/,$line);
		$geneDes{$data[1]} = $data[7];
		$transcript2gene{$data[1]} = $data[0];
		$genePos{$data[0]} = join("\t",@data[3..5]);
		$geneType{$data[0]} = $data[2];
		$id2Name{$data[1]} = $data[2];
	}
}


#read in all the chr_locus to gene annotation
sub readgff{
	my $gfffile = shift;
	open(GF,"$gfffile") || die;
	while(my $ln=<GF>){
		chomp $ln;
		my @dat = split(/\t/,$ln);
		if($dat[2] ne "CDS" and $dat[2] ne "3UTR" and $dat[2] ne "5UTR" and $dat[2] ne "rRNA" and $dat[2] ne "tRNA" and $dat[2] ne "snRNA" and $dat[2] ne "snoRNA"){
			next;
		}
		my $code = substr($dat[3],0,2);
		#gene_id "FBgn0031081"; gene_symbol "Nep3"; transcript_id "FBtr0307554";
		my @vvd = split(/[\",\;]/,$dat[8]);
		my $geneid = $vvd[1];
		my $trid = $vvd[7];
		$geneStrand{$trid} = $dat[6];
		#print $geneid."\t".$trid."\t".$dat[2],"\n";
		$chr2GR{$dat[0]}{$code}{$dat[3].'_'.$dat[4]}=$geneid."\t".$trid."\t".$dat[2];
	}
	close(GF);
}

sub mapChrGene{
	my $chr = shift;
	my $st = shift;
	my $ed = shift;
	my $ccd = substr($st,0,2);
	my $mapResult = "null";
	if(exists $chr2GR{$chr}{$ccd}){
		foreach my $key(keys %{$chr2GR{$chr}{$ccd}}){
			my ($a, $b) = split(/\_/,$key);
			if(overlap($st,$ed,$a,$b)>0){
				$mapResult = $chr2GR{$chr}{$ccd}{$key};
				last;
			}
		}
	}
	return($mapResult);
}

sub overlap{
	my ($a,$b,$c,$d)=@_;
	if($c>=$b){
		return($b-$c);
	} elsif($a>=$d){
		return($d-$a);
	} elsif($c>=$a and $c<=$b){
		if($d<=$b){
			return($d-$c);
		} else {
			return($b-$c);
		}
	} elsif($a>=$c and $a<=$d){
		if($b<=$d){
			return($b-$a);
		} else {
			return($d-$a);
		}
	}
}
