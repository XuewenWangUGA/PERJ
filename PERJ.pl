#!/usr/bin/perl -w
use strict;

# function: to join two illumina pair-end reads from two files , with 5 end triming options
#input data are the name of two fastq files and options
#output data are joined paired Illumina read with same illumina ID
#programmed by Xuewen Wang in 2013 Aug. Any question, email xwwang@ymail.com


# Usage: perl PERJ.pl -l LeftFastqFile -r RightFastqFile -o OUTfilename [options]
# e.g. perl PERJ.pl -l TTGGATGG_1_50.fq -r TTGGATGG_2_50.fq -trimleft 1 -trimright 2 -format 1 -n 0 -o testOUT.txt
	&usageinfo();
	
	# get command options
	my %commandin = @ARGV;
	if ((scalar @ARGV)%2 != 0){print "arguments must be in pair";}
	my $refinl=$commandin{"-l"}; # file name containing left (-l) illumina reads in fastq format
	my $refinr=$commandin{"-r"}; # file name containing right (-r) illumina reads in fastq format	
	my $outputfile=$commandin{"-o"};# file name to save out results
	my $trimlen_l=$commandin{"-trimleft"}||0; # triming length for left seq 
	my $trimlen_r=$commandin{"-trimright"}||0; # triming length for right seq
	my $readformat=$commandin{"-format"}||1; # 1|2 , format of utput reads, 1 for fasta format beginning with > and 2 for fastq (illumina-like) format beginning with @ 
	my $fillgap=$commandin{"-n"}||0; #fill in Ns between left and right reads, default 0 (no fill)
	#my $refinr="TTGGATGG_2_50.fq";  # right seq file for testing, 12 reads

# file status check and preparation
	open (DNAhandl,$refinl)|| die (" left seq file opening failed, please check");
	open (OUT,">$outputfile")|| die (" result file writing failed, please check");
	open (OUTsat,">$outputfile.sts")|| die (" result file writing failed, please check");
	
#definition
	my $word="";
	my %seq_hash=();
	my $join_ct=0; #joined reads count
	my $r_ct=0; #reads count
	my $l_ct=0; #reads count
	my $seq_r="";
	my $seq_l="";
	my $qua_l="";
	my $qua_r="";

	(my $hashref_r,my $quaref_r)=&onlyfasta(\$refinr); #get the right seq
	my %fastaseq_r=%{$hashref_r};
	my %fastaqua_r=%{$quaref_r};
	my $signc='';
	if($readformat==1){$signc="\>";}
	elsif($readformat==2){$signc="\@";}
	my $lc=0;
	my @block=();
	
	#generate Ns
	my $NS="";
	for (1..$fillgap) { $NS=$NS."N"; }
	
	# using line controlling
	while ($word=<DNAhandl>){ 	
		chomp $word;
		$lc ++;
		my $sectn=$lc%4;
		if($sectn==1){$block[0]=$word;} #id
		if($sectn==2){$block[1]=$word;} #seq
	 if($sectn==0){$block[3]=$word; #qual				
		$block[0]=~ s/^@//;
		$r_ct++;
		(my $ID, my $description)=split('\s',$block[0]);
		
		#get seq
		$seq_l=&trim5end($block[1],$trimlen_l); #trim left seq	
		my $len_l=length $seq_l;
		$seq_r=&trim3end($fastaseq_r{$ID},$trimlen_r)if(exists $fastaseq_r{$ID}); #trim right reccom seq 
		my $len_r=length $seq_r;
		
		#get quality
		$qua_l=&trim5end($block[3],$trimlen_l); #trim left quality		
		$qua_r=&trim3end($fastaqua_r{$ID},$trimlen_r)if(exists $fastaqua_r{$ID}); #trim right reccom quality 
		
		#print $ID,"\t",$block[1],"\t","$block[3]\n";
		if($len_l>0 and $len_r>0){
			print OUT $signc,$ID, " join$len_l\|$len_r","bp","\n", $seq_l,$NS, $seq_r,"\n" ;	
			if ($readformat==2){
				print OUT "+\n";
				print OUT $qua_l,$qua_r,"\n";}
			$join_ct ++;				
			$seq_r=""; #reset right seq to blank
		}
	  }
	}
	
	&runtime(\*OUTsat);# for running log
	print OUTsat "options used: -l $refinl -r $refinr -trimleft $trimlen_l -trimright $trimlen_r  -format $readformat -o $outputfile.\n";	
	print OUTsat "Total reads in $refinl are $l_ct .\n";
	print OUTsat "Total reads in $refinr are $r_ct .\n";
	print OUTsat "Total joined reads in $outputfile are $join_ct .\n";
	
	
	close DNAhandl;			
	close OUT;	

############################
	sub onlyfasta{
		my ($fref)=@_;
		my $word_r='';
		my @block_r=();
		my $lc_r=0;
		my %seq_hash=();
		my %qua_hash=();
		open (DNAhand,${$fref})|| die (" ${$fref} seq file opening failed, please check");
		while ($word_r=<DNAhand>){
			chomp $word_r;
			$lc_r ++;
			my $sectn_r=$lc_r% 4;
		if($sectn_r==1){$block_r[0]=$word_r;} #id
		if($sectn_r==2){$block_r[1]=$word_r;} #seq
		 if($sectn_r==0){$block_r[3]=$word_r; #qual		
			
			$block_r[0]=~ s/^@//;
			$l_ct++;
			(my $ID, my $description)=split('\s',$block_r[0]);
			#seq
			$block_r[1]=reverse $block_r[1]; #reverse
			$block_r[1]=~tr/ACGTacgt/TGCAtgca/; #complement
			$seq_hash{$ID}=$block_r[1];
			#quality
			$block_r[3]=reverse $block_r[3]; #reverse			
			$qua_hash{$ID}=$block_r[3];
		 }
		}
		return \%seq_hash,\%qua_hash;		
		close DNAhand;
	}	
	
	sub trim5end{
	#usage: &trim5end(inputSeq,trimlen)
		my ($inputseq,$trimlen)=@_;
		my $seqlen=length $inputseq; 
		my $startpos=0;
		if($trimlen>0){			
			$startpos=$trimlen;
		}
		my $extractlength=$seqlen-$trimlen;
		my $seq_trimed=substr($inputseq,$startpos,$extractlength);
		return $seq_trimed;		
	}
	
	sub trim3end{
	#for the  3end. here is the right half seq of reverse complement strand	
		my ($inputseq,$trimlen)=@_;
		my $seqlen=length $inputseq;
		my $startpos=0;		
		my $extractlength=$seqlen-$trimlen;
		my $seq_trimed=substr($inputseq,$startpos,$extractlength);
		return $seq_trimed;		
	}

	sub runtime() {
		my $OUTfile=shift @_;
		my $local_time = gmtime();
		print {$OUTfile} "$0 was programmed by Xuewen Wang, for supporting email : xwvan\@yahoo.com\n";
		print {$OUTfile} "$0 was run and results were yielded at $local_time\n";		
	} # end sub
	
	sub usageinfo{
		my @usage=(); # showing content on how to use the programme
		$usage[0]="Functions: 
join fastq pair-end reads such as Illunima reads. Tested OK in Windows, MacOS, and Linux.\n";
		$usage[1]="for    help: perl $0 ; or more in manual.pdf \n";
		$usage[2]="for running: perl $0 -l LeftFastqFile -r RightFastqFile -trimleft lefTrimLength -trimright RightTrimLength -format 1forFasta_2forfastq -o OUTfilename\n";
		$usage[4]="Author: Xuewen Wang. For supporting,email xwvan\@yahoo.com\n";		
		$usage[5]="Programmed in year 2013 Aug.\n";
		$usage[3]="
options:	value

-l:		the fastq file name of the left (l)/forward reads. if file name has space, put \"\" to either side of the file name
-r:		the fastq file name of the right (r) /reverse reads
-trimleft	integer number, the length of nt to be removed in the 5 end of the left/forward reads. e.g. -trimleft 5. default is 0.
-trimright	integer number, the length of nt to be removed in the 5 end of the right/reverse reads. e.g. -trimright 5. Default is 0. 
-format		value is 1 or 2. value 1 will produce joined read in fasta format. value 2 will produce joined read in fastq format. Default is 1.
-o		value is the file name to store the results of joined reads
-n		positive integer number, to fill in N for -n value times between left and right reads, default 0 (no filling)

result files: 
 file 1:  storing the results of joined reads, file name is given by user
 file 2:  statistic file with information of each input file and output file. file name with a suffix .sts after the file 1 name\n\n"; 
 $usage[6]="\nSurppoted illumina-like fastq Identifier format: 
e.g.1
\@CCRI0219:133:D243CACXX:7:1101:20008:1931 1:N:0:
e.g.2
\@EAS139:136:FC706VJ:2:2104:15343:197393 1:Y:18:ATCACG
e.g.3
\@SRR001666.1_071112_SLXA-EAS1_s_7:5:1:817:345 length=36
\n";
	unless(@ARGV){print @usage; exit;} 
 } #end sub	
exit;