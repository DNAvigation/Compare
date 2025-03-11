#!/usr/bin/perl

##############################################################################################
# 00.count.2021.2023.snv.pl
##############################################################################################
# 2025.02.
# goslak(Heo Jee Yeon; goslak@empal.com)
#############################################################################################
#ClinVar::make_ClinVar_vcf
#############################################################################################

#============================================================================================================
# 25.02.26
#============================================================================================================
# unique_2021_2023_clinvar
#============================================================================================================
sub unique_2021_2023_clinvar{
	print "=======================================================================================\n";
	print 1900+(localtime)[5]."-".(1+(localtime)[4])."-".(localtime)[3]."[".(localtime)[2].":".(localtime)[1].":".(localtime)[0]."]\n";

	$para{beforeY} = 20201226;
	$para{afterY} = 20231230;
	$para{in_fold} = "/home/dna/03_data/ClinVar/00_year";
	$para{prefix} = "clinvar_";
	
	&filter_2021_2023_snvs(\%para);
	#&diff_ClinVar_btw_version(\%para);

	print 1900+(localtime)[5]."-".(1+(localtime)[4])."-".(localtime)[3]."[".(localtime)[2].":".(localtime)[1].":".(localtime)[0]."]\n";
	print "=======================================================================================\n";
}
#============================================================================================================
&unique_2021_2023_clinvar();


#============================================================================================================
# 25.02.26
#============================================================================================================
# filter_2021_2023_snvs
#============================================================================================================
# IN_FILE
#------------------------------------------------------------------------------------------------------------
# clinvar_*.vcf
#------------------------------------------------------------------------------------------------------------
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
#1	1041249	rs113789806	C	T	.	.	RS=113789806;RSPOS=1041249;dbSNPBuildID=132;SSR=0;SAO=0;VP=0x0501600a0305140136100100;GENEINFO=AGRN:375790;WGT=1;VC=SNV;PM;SLO;REF;SYN;INT;R5;ASP;VLD;GNO;KGPhase1;KGPhase3;LSD;CLNALLE=1;CLNHGVS=NC_000001.11:g.1041249C>T;CLNSRC=ClinVar;CLNORIGIN=1;CLNSRCID=NM_198576.3:c.804C>T;CLNSIG=2;CLNDSDB=.;CLNDSDBID=.;CLNDBN=not_specified;CLNREVSTAT=single;CLNACC=RCV000116282.1
#============================================================================================================
sub filter_2021_2023_snvs{
	print 1900+(localtime)[5]."-".(1+(localtime)[4])."-".(localtime)[3]."[".(localtime)[2].":".(localtime)[1].":".(localtime)[0]."]\n";
	print "----------------------------------------------------------------------\n";
	print "* sub : filter_2021_2023_snvs\n";
	print "----------------------------------------------------------------------\n";
	
	my %CLNSIG = ("Benign"=>1,"Likely_benign"=>2,"Benign/Likely_benign"=>3,
				"Pathogenic"=>4,"Likely_pathogenic"=>5,"Pathogenic/Likely_pathogenic"=>6);
	
	my %CLNREVSTAT = ("practice_guideline"=>4, "reviewed_by_expert_panel"=>3, 
					"criteria_provided,_multiple_submitters,_no_conflicts"=>2, "criteria_provided,_single_submitter"=>1, 
					"no_interpretation_for_the_single_variant"=>0, "criteria_provided,_conflicting_interpretations"=>0,
					"no_assertion_criteria_provided"=>0, "no_assertion_provided"=>0);


	my (%dup)=();
	foreach $year (qw(20201226 20231230)){
		my $in = $_[0]->{in_fold}."/".$_[0]->{prefix}.$year.".vcf";
		open(IN, $in) or die "can't read $in: $!\n";
		
		my (%hash)=();
		while(<IN>){
			chomp;
			s/\r//g;
			
			if(/^\#/){
				next;
			}
			
			$hash{in}++;
			my @temp = split(/\t/, $_);
			#$hash{in_0_chr_.$temp[0]}++;
			#$hash{in_5_qual_.$temp[5]}++;
			#$hash{in_6_filter_.$temp[6]}++;
			
			if ($temp[0] !~ /^(?:[1-9]|1[0-9]|2[0-2]|X|Y|MT)$/){
				$hash{next_chr_.$temp[0]}++;
				next;
			}
			
			for $i (qw(3 4)){
				if($temp[$i] eq '.'){
					$hash{next_blank_.$i}++;
					next;
				}
			}
			
			if($temp[7]=~/CLNVCSO=SO:0001483/){
				$hash{cnt_snv}++;
				
				my $pos = join("_", @temp[0,1,3,4]);
				unless($dup{pos_.$pos}){
					$hash{dup1_.$year}++;
					$dup{pos_.$pos} = $year;
				}else{
					$hash{dup2_.$year}++;
					#print "!!! dup(pos):\t$pos\n";
					next;
				}
			}

			if ($temp[7]=~/CLNSIG=([^;]+)/){
				my $v = $1;
				if (exists $CLNSIG{$v}){
					$hash{ok_clnsig_exist}++;
					$hash{ok_clnsig_.$CLNSIG{$v}}++;
				}else{
					$hash{next_clnsig}++;
					next;
				}
			}
			
			if ($temp[7]=~/CLNREVSTAT=([^;]+)/){
				my $v = $1;
				if (exists $CLNREVSTAT{$v}){
					$hash{ok_CLNREVSTAT_exist}++;
					$hash{ok_CLNREVSTAT_.$CLNREVSTAT{$v}}++;
					if ($CLNREVSTAT{$v}>=2){
						$hash{ok_CLNREVSTAT_exist234}++;
					}
				}else{
					$hash{next_CLNREVSTAT}++;
					next;
				}
			}
		}
		close(IN);
	
		print "\tIn:\t$in\n"; 
		foreach $key (sort keys %hash){
			print "\t".$key."\t".$hash{$key}."\n";
		}
	}
	print "----------------------------------------------------------------------\n";
	print 1900+(localtime)[5]."-".(1+(localtime)[4])."-".(localtime)[3]."[".(localtime)[2].":".(localtime)[1].":".(localtime)[0]."]\n";
	print "\n\n";
}
#============================================================================================================

