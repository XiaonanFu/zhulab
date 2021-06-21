#!/usr/bin/perl -w 
use strict;
use Getopt::Long;

# ------------------------------------------------------------------
# Usage: This code was design for identify miRNA-mRNA hybrids
#        from fastq sequencing files. 01-11-16
# Require: PEAR - Paired-End reAd mergeR
#		   (http://sco.h-its.org/exelixis/web/software/pear/)
#		   Blast - NCBI blast software
#		   (ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/)
#		   Flexbar - flexible barcode and adapter removal for sequencing platforms
#		   (http://sourceforge.net/projects/flexbar/)
#        email: xnfu@vt.edu
# ------------------------------------------------------------------

my $usage = <<USAGE;

 perl clashCallAngam.pl <parameters>
  -f 	<str>	fastq file of left sequences
  -r	<str>	fastq file of right sequences
  -a	<str>	annotation file

USAGE

my $fastqLeft = "";
my $fastqRight = "";
my $annotation = "";

GetOptions ("f=s" => \$fastqLeft,
			"r=s" => \$fastqRight,
			"a=s" => \$annotation
			);

if ($fastqLeft eq "" or $fastqRight eq "" or $annotation eq "") {
    print $usage,"\n";
    exit;
}

my $transcriptFile = "transupdateP44.fa";
my $chromosomeFile = "AgamChr.fa";
my $sampleName = [split(/\_/,$fastqLeft)]->[0];
my %fastaSample;
my %updateFa;
my %clash;
my %geneDes;
my %genePos;
my %geneType;
#=================================main==================================
my $datestring = localtime();
print "Start running,the time now is $datestring\n";

#step 1: Prepare DB for blast
#prepareDBforBlast();

#step 2: Clean up the sequenced reads from fastq files
print "fastq cleaning...\n";
#fastqClean();

#step 3: Mapping the reads to transcripts and genome using blast
print "Mapping reads...\n";
readsMapping();
#cleanPCRbias();

#step 4: Call the chimeric reads for interactions between miRNAs and mRNAs
#parameters: gap<=8
print "calling and output...\n";
chimericCall(8);

print "Annotate the result...\n";
readAnnotation($annotation);
resultOutput();

$datestring = localtime();
print "Running finished at $datestring\n";

#=================================function==============================
#prepare indexed database for blast
sub prepareDBforBlast{
	#system("makeblastdb -in $transcriptFile -input_type fasta -dbtype nucl -parse_seqids -out AgamP44Uv2");
	#system("makeblastdb -in $chromosomeFile -input_type fasta -dbtype nucl -parse_seqids -out AgamChr");
}

#clean up the raw fastq file including merging paired sequences;
#remove bad quality reads; remove the PCR duplicates
sub fastqClean{
	#-f forward fastq file, -r reverse fastq file, -n minimum length after merging, -j threads
	system("pear -f $fastqLeft -r $fastqRight -o $sampleName -n 15 -j 8");
	system("cat ${sampleName}.assembled.fastq ${sampleName}.unassembled.forward.fastq >${sampleName}.fastq");
	system("flexbar -t ${sampleName}.clip1 -f sanger -r ${sampleName}.fastq -as AGATCGGAAGAGCACACGTCT -ao 4 -u 3 -m 16 -n 6");
	system("fastx_collapser -Q 33 -i ${sampleName}.clip1.fastq -o ${sampleName}.fa");
	system("flexbar -r ${sampleName}.fa -t ${sampleName}index -as NNNTAAGC -ae LEFT_TAIL -d -m 16");
}

sub readsMapping{
	readFA("${sampleName}index.fasta");
	#system("makeblastdb -in ${sampleName}index.fasta -input_type fasta -dbtype nucl -parse_seqids -out ${sampleName}index");
	#system("blastn -query aga_matureUpdate032616.fa -db ${sampleName}index -task blastn -evalue 0.4 -outfmt '6 qseqid sseqid qlen pident length mismatch gapopen qstart qend sstart send evalue bitscore btop' -num_threads 8 -max_target_seqs 10000000 -out ${sampleName}indexmiBlast.out");
	processmiRNAmapping("${sampleName}indexmiBlast.out");
	#system("blastn -query ${sampleName}_transBlast.fa -db AgamP44Uv3 -task blastn -evalue 0.4 -outfmt '6 qseqid sseqid qlen pident length mismatch gapopen qstart qend sstart send evalue bitscore btop' -num_threads 8 -max_target_seqs 10 -out ${sampleName}index.Tout1");
	blastTranscriptFilter("${sampleName}index.Tout1","${sampleName}index.fasta");
}

#read in the fasta file of sample
sub readFA{
	my $file = shift;
    open(FA, $file) or die;
    my $id;
    while(my $ln=<FA>){
        chomp $ln;
        if($ln=~/^\>/){
            $id = $ln;
            $id=~s/\>//g;
        } else {
            $fastaSample{$id}=$ln;
        }
    }
    close(FA);
}

#filter the miRNA mappping results
sub processmiRNAmapping{
	my $miRNAblast = shift;
	open(BLAST,"$miRNAblast") || die;
	my $id;
	my %miRNAhit;
	while(my $ln=<BLAST>){
		chomp $ln;
		my @data = split(/\t/,$ln);
		if(countMismatch($data[13])>2){ next;}
		if(not exists $miRNAhit{$data[1]}){
			$miRNAhit{$data[1]} = $ln;
		} else {
			my $current_score = [split(/\t/,$miRNAhit{$data[1]})]->[12];
			if($current_score>$data[12]){
				$miRNAhit{$data[1]}=$ln;
			}
		}
	}
	close(BLAST);
	
	open(OUT,">${sampleName}_transBlast.fa") || die;
	open(FINAL,">${sampleName}.final") || die;
	foreach my $key(keys %fastaSample){
		if(exists $miRNAhit{$key}){
			my $read_len = length($fastaSample{$key});
			my @blast = split(/\t/,$miRNAhit{$key});
			my ($mistart,$miend);
			if($blast[4]/$read_len >=0.9){
				print FINAL $miRNAhit{$key},"\n";
			} else {
				if($blast[9]>$blast[10]){
					($mistart,$miend) = ($blast[10],$blast[9]);
				} else {
					($mistart,$miend) = ($blast[9],$blast[10]);
				}
				if(($mistart-1)>=16 and ($read_len-$miend)<($mistart-1)){
					print FINAL $miRNAhit{$key},"\n";
					print OUT '>',$key,"\n",substr($fastaSample{$key},0,($mistart-1)),"\n";
					$updateFa{$key} = substr($fastaSample{$key},0,($mistart-1));
				} elsif(($read_len-$miend)>=16 and ($read_len-$miend)>($mistart-1)){
					print FINAL $miRNAhit{$key},"\n";
					print OUT '>',$key,"\n",substr($fastaSample{$key},$miend,($read_len-$miend)),"\n";
					$updateFa{$key} = substr($fastaSample{$key},$miend,($read_len-$miend));
				} else {
					print OUT '>',$key,"\n",$fastaSample{$key},"\n";
					$updateFa{$key} = $fastaSample{$key};
				}
			}
		} else {
			print OUT '>',$key,"\n",$fastaSample{$key},"\n";
			$updateFa{$key} = $fastaSample{$key};
		}
	}
	close(OUT);
	close(FINAL);
}

#Analyze the blast results against Transcripts;
#Output read sequences which are not fully mapped to transcripts (<90% identity)
#try to map them on chromosome
sub blastTranscriptFilter{
	my $blast = shift;
	my $fa = shift;
	my %mapping;
	my %hitscore;
	my %read;
	#step1 rough filter
	open(BLAST, "$blast") || die;
	while(my $ln=<BLAST>){
		chomp $ln;
		my @data = split(/\t/,$ln);
		#for every blast result, pick all
		if(exists $mapping{$data[0]}){
			my @val = @{$mapping{$data[0]}};
			$mapping{$data[0]}[$#val+1] = $ln;
		} else {
			$mapping{$data[0]}[0] = $ln;
		}
	}
	#step2 count frequency of transcript with clearly and uniquly mapped reads
	#this will be used as priority if blast hits with same parameter
	my %freq;
	foreach my $key(sort keys %mapping){
		my @list = @{$mapping{$key}};
		my $score = [split(/\t/,$list[0])]->[12];
		my @index;
		my $j = 0;
		for(my $i=0;$i<=$#list;$i++){
			my @para = split(/\t/,$list[$i]);
			my $errorRate = countMismatch($para[13])/$para[4];
			if($errorRate<=0.1 and $para[12]>=$score){
				$index[$j++] = $i;
			}
		}
		if($j>0){
			foreach my $key(@index){
				my $id = [split(/\t/,$list[$key])]->[1];
				$freq{$id}++;
			}
		} else {
			next;
		}
	}
	#step3 assign each read with one loci
	my $mapcount = 0;
	my %cleanBlast;
	foreach my $key(sort keys %mapping){
		my @data = @{$mapping{$key}};
		$cleanBlast{$key}[0]=$data[0];
		if($#data>0){
			foreach my $alignment(@data[1..$#data]){
				my @fcomp = split(/\t/,$alignment);
				my @check = @{$cleanBlast{$key}};
				my $countOP = 0;    #check if their exists overlapping
				my $index = $#check;
				for(my $i=0;$i<=$#check;$i++){
					my @wcomp = split(/\t/,$check[$i]);
					if(overlap($fcomp[7],$fcomp[8],$wcomp[7],$wcomp[8])>=5){
						$countOP++;
						if($wcomp[12]<$fcomp[12]){
							$cleanBlast{$key}[$i]=$alignment;
							next;
						} elsif($wcomp[12]==$fcomp[12]){
							if($fcomp[1]=~/aga\-/ and $wcomp[1]!~/aga\-/){
								$cleanBlast{$key}[$i]=$alignment;
								last;
							}
							if(not exists $freq{$wcomp[1]}){
								$freq{$wcomp[1]} = 0.2;
							}
							if(not exists $freq{$fcomp[1]}){
								$freq{$fcomp[1]} = 0.2;
							}
							if($freq{$wcomp[1]}<$freq{$fcomp[1]} and $wcomp[1]!~/aga\-/){
								$cleanBlast{$key}[$i]=$alignment;
								last;
							}
						}
					}
				}
				if($countOP==0) {
					$cleanBlast{$key}[$index+1]=$alignment;
				}
			}
		}
	}
	open(BFINAL,">>${sampleName}.final") || die "can't open $!";
	foreach my $key(keys %cleanBlast){
		my @dat = @{$cleanBlast{$key}};
		foreach my $val(@dat){
			print BFINAL $val,"\n";
		}
	}
	close(BFINAL);
}

#count the insertion+deletion and mutation
sub countMismatch{
	my $string = shift;
	my @val = split(/\d+/,$string);
	my $mis = "";
	for(my $i=1;$i<=$#val;$i++){
		$mis.=$val[$i];
	}
	my $count = length($mis)/2;
	return($count);
}


#To identify chimeric reads contains miRNAs and other type of RNA
#overlapping <=4
#output the potential CLIP tags
sub chimericCall{
	my $blast = "${sampleName}.final";
	my $gap = shift;
	my %mapping;
	open(BLAST, "$blast") || die;
	my $id;
	while(my $ln=<BLAST>){
		chomp $ln;
		if($ln=~/^aga\-/){
			$id = [split(/\t/,$ln)]->[1];
		} else {
			$id = [split(/\t/,$ln)]->[0];
		}
		if (not exists $mapping{$id}){
			$mapping{$id}[0]=$ln;
		} else {
			my @inval = @{$mapping{$id}};
			$mapping{$id}[$#inval+1]=$ln;
		}
	}
	#call for the clash chemirc reads
	foreach my $key(keys %mapping){
		my @data = @{$mapping{$key}};
		if($#data==0){	next;}   #unique mapping, none of chemirc read exists
		my $mi="";
		my $mm="";
		my $bfTag = 'i'; #determine if the miRNA mapping happened in previous step or not
		for(my $i=0;$i<=$#data;$i++){
			if($data[$i]=~/aga\-/){
				my @tem = split(/\t/,$data[$i]);
				if($tem[1]=~/^aga\-/){
					$bfTag = 'a';
					($tem[0],$tem[1]) = ($tem[1],$tem[0]);
				} else {
					$bfTag = 'i';
				}
				$mi = join("\t",@tem);
				splice @data,$i,1;
				last;
			}
		}
		if($mi){
			my @current = split(/\t/,$mi);
			my $count = 0;
			foreach my $align(@data){

				my @test = split(/\t/,$align);
				#extract the mapped position for mm
				my $MMmatchSeq = substr($updateFa{$test[0]},$test[7]-1,abs($test[8]-$test[7])+1);
				my $MMmatchStart = index($fastaSample{$test[0]},$MMmatchSeq)+1;
				my $MMmatchEnd = $MMmatchStart+length($MMmatchSeq)-1;
				
				#extract the mapped position for mi
				my ($MImatchStart,$MImatchEnd);
				if($bfTag eq 'i'){
					if($current[9]>$current[10]){
						($MImatchStart,$MImatchEnd)= ($current[10],$current[9]);
					} else {
						($MImatchStart,$MImatchEnd)= ($current[9],$current[10]);
					}
				} else {
					($MImatchStart,$MImatchEnd)= ($current[7],$current[8]);
				}
				if(abs(overlap($MMmatchStart,$MMmatchEnd,$MImatchStart,$MImatchEnd))>$gap){
					next;
				} else {
					my $direction = 'F';
					if((($MMmatchStart+$MMmatchEnd)/2)<(($MImatchStart+$MImatchEnd)/2)){
						$direction = 'R';
					}
					$clash{$key}[$count++]=$mi."\t".$align."\t$direction";
				}
			}
		} else {
			next;
		}
	}
}

#cllapse the overlapping chimeric reads
sub resultOutput{
	my %MiMmPairs;
	foreach my $key(keys %clash){
		foreach my $value(@{$clash{$key}}){
			my @data = split(/\t/,$value);
			my $tagID = $data[1];
			my $mirna = $data[0];
			my $mmrna = $data[15];
			my $mmAlignStatus = $data[20]+$data[21];
			my $mmAlignStart = $data[23];
			my $mmAlignEnd = $data[24];
			my $mmStrand = '+';
			my $mmDirection = $data[28];
			my $miStrand = '+';
			my $miAlignStart = $data[9];
			my $miAlignEnd = $data[10];
			#determine the match strand
			if($mmAlignEnd < $mmAlignStart){
				$mmStrand = '-';
				($mmAlignEnd,$mmAlignStart)=($mmAlignStart,$mmAlignEnd);
			}
			if($miAlignEnd < $miAlignStart){
				$miStrand = '-';
			}
			#remove reverse reads-mRNA mapping
			if($mmStrand ne $miStrand){
				next;
			}
			if(!exists $MiMmPairs{$mirna}{$mmrna}){
				$MiMmPairs{$mirna}{$mmrna}[0]{'unique'} = 1;
				$MiMmPairs{$mirna}{$mmrna}[0]{'taglist'} = $tagID;
				$MiMmPairs{$mirna}{$mmrna}[0]{'AlignStatus'} = $mmAlignStatus;
				$MiMmPairs{$mirna}{$mmrna}[0]{'AlignStart'} = $mmAlignStart;
				$MiMmPairs{$mirna}{$mmrna}[0]{'AlignEnd'} = $mmAlignEnd;
				$MiMmPairs{$mirna}{$mmrna}[0]{'miStrand'} = $miStrand;
				$MiMmPairs{$mirna}{$mmrna}[0]{'mmStrand'} = $mmStrand;
				$MiMmPairs{$mirna}{$mmrna}[0]{'Direction'} = $mmDirection;
			} else {
				my @data = @{$MiMmPairs{$mirna}{$mmrna}};
				my $index = $#data;
				my $countCllapse = 0;
				for(my $j=0;$j<=$#data;$j++){
					my $currentStart = $MiMmPairs{$mirna}{$mmrna}[$j]{'AlignStart'};
					my $currentEnd = $MiMmPairs{$mirna}{$mmrna}[$j]{'AlignEnd'};
					my $overmatch = overlap($mmAlignStart,$mmAlignEnd,$currentStart,$currentEnd);
					if($overmatch>=0){
						$MiMmPairs{$mirna}{$mmrna}[$j]{'unique'}++;
						$MiMmPairs{$mirna}{$mmrna}[$j]{'taglist'}.='|'.$tagID;
						if($currentStart > $mmAlignStart){
							$MiMmPairs{$mirna}{$mmrna}[$j]{'AlignStart'} = $mmAlignStart;
						}
						if($currentEnd < $mmAlignEnd){
							$MiMmPairs{$mirna}{$mmrna}[$j]{'AlignEnd'} = $mmAlignEnd;
						}
						$countCllapse++;
						last;
					}
				}
				if($countCllapse==0){
					$MiMmPairs{$mirna}{$mmrna}[$index+1]{'unique'} = 1;
					$MiMmPairs{$mirna}{$mmrna}[$index+1]{'taglist'} = $tagID;
					$MiMmPairs{$mirna}{$mmrna}[$index+1]{'AlignStatus'} = $mmAlignStatus;
					$MiMmPairs{$mirna}{$mmrna}[$index+1]{'AlignStart'} = $mmAlignStart;
					$MiMmPairs{$mirna}{$mmrna}[$index+1]{'AlignEnd'} = $mmAlignEnd;
					$MiMmPairs{$mirna}{$mmrna}[$index+1]{'miStrand'} = $miStrand;
					$MiMmPairs{$mirna}{$mmrna}[$index+1]{'mmStrand'} = $mmStrand;
					$MiMmPairs{$mirna}{$mmrna}[$index+1]{'Direction'} = $mmDirection;
				}
			}
		}
	}
	open(OUTCLASH,">${sampleName}clash.test") || die;
	foreach my $key(sort keys %MiMmPairs){
		foreach my $sky(sort keys %{$MiMmPairs{$key}}){
			my @content = @{$MiMmPairs{$key}{$sky}};
			for(my $i=0;$i<=$#content;$i++){
				my %containValue = %{$content[$i]};
				my $biotype;
				if(not exists $geneType{$sky}){
					$biotype = "undetermined";
				}elsif($geneType{$sky} ne 'mRNA'){
					$biotype = $geneType{$sky};
				} else {
					$biotype = bioType($sky,$containValue{'AlignStart'}."\t".$containValue{'AlignEnd'});
				}
				if(not exists $geneDes{$sky}){
					$geneDes{$sky} = 'None annotation';
				}
				print OUTCLASH "$key"."\t".$sky."\t".$biotype,"\t",$containValue{'Direction'},"\t",$containValue{'unique'},"\t".$containValue{'AlignStart'}."\t".$containValue{'AlignEnd'}."\t".$containValue{'miStrand'}."\t".$containValue{'mmStrand'}."\t".$containValue{'taglist'}."\t".$geneDes{$sky}."\n";
			}
		}
	}
}

#To count the overlapping nucleotides betweed two alignments
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

#determine the biogenetic region by chromosome coordinates
sub bioType{
	my $id = shift;
	my $str = shift;
	my ($start, $end) = split(/\t/,$str);
	my ($a,$b,$c)= split(/\t/,$genePos{$id});
	my @utrcds;
	$utrcds[0] = 1;
	$utrcds[1] = $a;
	$utrcds[2] = $a+1;
	$utrcds[3] = $a+$b;
	$utrcds[4] = $a+$b+1;
	$utrcds[5] = $a+$b+$c;
	if($start>=$utrcds[0] and $start<=$utrcds[1]){
		return '5UTR';
	} elsif($start>=$utrcds[2] and $start<=$utrcds[3]){
		return 'CDS';
	} else{
		return '3UTR';
	}
}

#read in all the annotation
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
}
