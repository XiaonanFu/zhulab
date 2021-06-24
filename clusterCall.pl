#!/usr/bin/perl -w
use strict;

my %blast;
my %cluster;
my %cluster_peak;

my %minPos;
my %maxPos;

readBlast($ARGV[0]);
clusterAlign();
preprocess_cluster();
peakCalling();





sub readBlast{
	my $file = shift;
	open(BLAST,"$file") || die;
	my $testcout=0;
	while(my $ln=<BLAST>){
		chomp $ln;
		my @content = split(/\t/,$ln);
		my $id = $content[1];
		my $s = $content[9];
		my $e = $content[10];
		if($s>$e){
			($s,$e) = ($e,$s);
		}
		if(not exists $minPos{$id}){
			$minPos{$id}=$s;
		} elsif($minPos{$id}>$s){
			$minPos{$id}=$s;
		}
		if(not exists $maxPos{$id}){
			$maxPos{$id}=$e;
		} elsif($maxPos{$id}<$e){
			$maxPos{$id}=$e;
		}
		for(my $i=$s;$i<=$e;$i++){
			$blast{$id}{$i}++;
		}
	}
	close(BLAST);
}

sub clusterAlign{
	foreach my $key(sort keys %blast){
		my $cluster_index = 0;
		my $seg_index = 0;
		my $pre_cov;
		for(my $i=$minPos{$key};$i<=$maxPos{$key};$i++){
			if(not exists $blast{$key}{$i}){
				next;
			}
			if(not exists $blast{$key}{$i-1} and exists $blast{$key}{$i}){
				$cluster_index ++;
				$seg_index = 0;
				$pre_cov = $blast{$key}{$i};
			}
			if($blast{$key}{$i}==$pre_cov){
				$cluster{$key}{$cluster_index}{$seg_index}{$i} = $blast{$key}{$i};
			} else {
				$seg_index++;
				$cluster{$key}{$cluster_index}{$seg_index}{$i} = $blast{$key}{$i};
			}
			$pre_cov = $blast{$key}{$i};
		}
	}
}


sub preprocess_cluster{
	foreach my $key(keys %cluster){
		foreach my $sky(sort {$a<=>$b} keys %{$cluster{$key}}){
			foreach my $tky(sort {$a<=>$b} keys %{$cluster{$key}{$sky}}){
				my @order = sort {$a<=>$b} keys  %{$cluster{$key}{$sky}{$tky}};
				my $first = $order[0];
				my $last = $order[$#order];
				my $index = int(($first+$last)/2.0);
				my $height = $cluster{$key}{$sky}{$tky}{$first};
				if(not exists $cluster_peak{$key}{$sky}{'coord'}){
					$cluster_peak{$key}{$sky}{'coord'} = $index;
				} else {
					$cluster_peak{$key}{$sky}{'coord'} .= "\t".$index;
				}
				if(not exists $cluster_peak{$key}{$sky}{'height'}){
					$cluster_peak{$key}{$sky}{'height'} = $height;
				} else {
					$cluster_peak{$key}{$sky}{'height'} .= "\t".$height;
				}
				$cluster_peak{$key}{$sky}{'segs'}++; #define how manys segs are in the cluster
			}
		}
	}
}


sub peakCalling{
	foreach my $key(keys %cluster_peak){
		my $index = 0;
		foreach my $sky(sort {$a <=> $b} keys %{$cluster_peak{$key}}){
			#require more than 3 segs for cubic spline analysis
			if($cluster_peak{$key}{$sky}{'segs'}<3){
				next;
			}
			#print $key,"\t",$sky,"\t",$cluster_peak{$key}{$sky}{'segs'},"\n";
			my @xx = split(/\t/,$cluster_peak{$key}{$sky}{'coord'});
			my @yy = split(/\t/,$cluster_peak{$key}{$sky}{'height'});
			print $key,"\t";
			print join("\,",@xx),'|';
			print join("\,",@yy),"\n";
		}
	}
}

