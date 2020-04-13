#! /usr/bin/env perl
use Bio::SeqIO; 
use Bio::Seq; 
# the script was modified to wor with the pipeline output
# the infile in this case is the pipeline output
my($cohort,$genome) = @ARGV;
chomp($cohort);
chomp($genome);

### static files ######
#$genome = "/shareDB/genome/human/human_g1k_v37.fasta";
$miRNA_data="/home/labs/hornsteinlab/chenc/annovar/3UTR_gain/lib/miRNA_seed_data.txt";
$genesStrandsInfoFile = "/home/labs/hornsteinlab/chenc/annovar/3UTR_gain/lib/genes_strand.txt";
### end of static files #######

# input output files
$annovar =$cohort."_annot.txt";
$outputF = $annovar.".1";
$variationFile = 'vars_3UTR.txt';

$S=-8; # index, how many bases around variation postision
$StrL=8;



## select 3UTR positions
## generate a file called vars_3UTR.txt
# criteria for SNP selection:
$UTR3=1;
$targetScan=0;
(%selected)=extract_SNPs_from_annovar($annovar,$variationFile);

## SNPs
## data of each SNP: chromosome (no chr), position, gene, change
## returns a hash called VAR
load_position($variationFile);
print "finished to load positions\n";
## We need to know the strand of each gene
## hash will be gene=> strand
load_gene_strand($genesStrandsInfoFile);

## load miRNA recognition site
load_miRNA_site($miRNA_data);
print "finished mir calculation\n";
## do the main calculation
main_calc();

# parse output of previous run= remove redundancy
# aggregate mirs per SNP
parse_output();

# write final output
write_final_output($annovar,$outputF);

sub write_final_output{
	my($input,$output)= @_;
	#read in annovar annotation. add a column to the end of file with 3UTR gain
	open(OUT,">$output")|| warn "can not write to out file of selected SNP annotation\n";
	
	
	open(IN,"$input")|| warn "no input from annovar\n";
	# parse header
	$line = <IN>;
	chomp($line);
	print OUT "$line\t3UTR_gain\n";
	$index = 0;
	my(@headerVar) = split(/\t/,$line);
	foreach $item (@headerVar){
		$item =~ s/\./_/g;
		$annotation{$item}=$index;
		$index++;
	}
	while($line = <IN>){
		chomp($line);
		$line1=$line; # we keep original line for output
		my(@data) = split(/\t/,$line1);
		$key="$data[$annotation{CHROM}]"."_"."$data[$annotation{POS}]"."_"."$data[$annotation{ALT}]";
		$mirS = '';
		
		foreach $mir (keys %{$mirs{$key}}){
			$mirS.= "$mir:$order[$mirs{$key}{$mir}{type}-1];";
		#$mir:$order[$mirs{$mutation}{$mir}{type}-1];
		}
	
		print OUT "$line$mirS\n";
		
	}
	close(IN);
	close(OUT);
	return();
}
sub parse_output{
	# parse output
	#order - 8mer,7mer-m8,7mer-1A
	@order = ('7mer-1A','7mer-m8_s','7mer-m8','8mer_s','8mer');
	$typeC{'8mer'}=5;
	$typeC{'8mer_s'}=4;
	$typeC{'7mer-m8'}=3;
	$typeC{'7mer-m8_s'}=2;
	$typeC{'7mer-1A'}=1;

	open(IN,"temp_3utr_gain.txt")|| warn "can not find temporary file\n";
	$line = <IN>;
	while($line = <IN>){
		chomp($line);
		#Gene,chr,position,alt_allele,var_type,mir,merSeq,WTSeq,MutSeq,merType,merPosition
		@data = split(/\,/,$line);
		$key = "$data[1]"."_"."$data[2]"."_"."$data[3]";
		if($mirs{$key}{$data[5]}{type} > 0 ){# if already defined, we keep the higher
		
			if($typeC{$data[9]} > $mirs{$key}{$data[5]}{type}){
			
				$mirs{$key}{$data[5]}{type}=$typeC{$data[9]};

			}
		
		}else{
			$mirs{$key}{$data[5]}{type}=$typeC{$data[9]};
			
		}
		
		
	}
	close(IN);
	open(OUT,">test_3UTR.txt");
	## now we aggregate data
	foreach $mutation (keys %mirs){
		print OUT "$mutation\t";
		foreach $mir (keys %{$mirs{$mutation}}){
			print OUT "$mir:$order[$mirs{$mutation}{$mir}{type}-1];";
		}
		print OUT "\n";
	}
	close(OUT);
	return();
}
######################################################################################################################
# sub with main calculation
sub main_calc{
	my $genomeSeq = Bio::SeqIO->new(-file   => "$genome",
									-format => 'fasta');			
	open(OUT,">temp_3utr_gain.txt");# temporary file with 3UTR gain
	print OUT "Gene,chr,position,alt_allele,var_type,mir,merSeq,WTSeq,MutSeq,merType,merPosition\n";

	while (my $seq = $genomeSeq->next_seq) {
	$chrNum= $seq->id,"\n";
		GO:foreach $variation (keys %{$VAR{$chrNum}}){
			#calculate 8 bases word coordinates
		
			$loc = $variation;
			$S=-8;
			@OriSeq = ();
			@OriSeq7 = ();
			@$OriSeq7_1=();
			
			@seqToCheck=();
			@seqToCheck7A=();
			
			GO1:for($i=1; $i<9;$i++){

				if($VAR{$chrNum}{$variation}{type} eq 'snp'){
					get_snp_seq($seq);
				}elsif($VAR{$chrNum}{$variation}{type} eq 'del'){
					get_del_seq($seq);
				}elsif($VAR{$chrNum}{$variation}{type} eq 'ins'){
					get_ins_seq($seq);
				}# end extract sequence
			
			
				# reverse complement if required
				if($genesStrand{$VAR{$chrNum}{$variation}{gene}} eq '-'){
					if($seqToCheck[$i] =~ /[AGCT]/){
						$seqToCheck[$i]=reverseSeq($seqToCheck[$i]);
					}
					if($OriSeq[$i] =~ /[AGCT]/){
						$OriSeq[$i] =reverseSeq($OriSeq[$i]);
					}
					
				}elsif($genesStrand{$VAR{$chrNum}{$variation}{gene}} eq '+'){
				}else{
					print "PROBLEM IN STRAND $genesStrand $VAR{$chrNum}{$variation}{gene}\n";
				}
			
				# extract 7A
				if($i==1){
					$seqToChek7A[$i] = substr($seqToCheck[$i],0,7);
					$OriSeq7[$i] = substr($OriSeq[$i],0,7);
					$OriSeq7_1[$i] = substr($OriSeq[$i],1,7);

				}elsif($i < 8){
					$seqToChek7A[$i] = substr($seqToCheck[$i],0,7);
					$OriSeq7[$i] = substr($OriSeq[$i],0,7);
					$OriSeq7_1[$i] = '';
		
				}else{
					$OriSeq7[$i] = substr($OriSeq[$i],0,7);
					$OriSeq7_1[$i] = '';
					$seqToChek7A[$i] = '';

				}
				$S++;
			}#end GO1
			
			# look what mirs bind to wildtype
			check_wt();
		
			# check mutated
			check_mut(\*OUT);
		}
			
	}#end chr load	
	close(OUT);
return();
}

sub reverseSeq{
	my($seq) = @_;
	
	$tempSeq = Bio::Seq->new(-seq => $seq);
	$tempSeq = $tempSeq->revcom();
	$newSeq = $tempSeq ->seq();
	
	return($newSeq);
}
sub load_gene_strand{
	my($inFile) = @_;
	open(IN,$inFile)|| warn "can not locate file with genes strand info\n";
	$line = <IN>;
	while ($line = <IN>){
		chomp($line);
		my($gene,$strand) = split(/\t/,$line);
		$genesStrand{$gene} =$strand;
	}
	close(IN);
	
	return();
	
}
sub load_miRNA_site{
	my($inFile) = @_;
	open(IN,$inFile)|| warn "can not locate miRNA seed details\n";
	$line = <IN>;
	open(OUT,">mir_processed.txt");
	while ($line = <IN>){
		chomp($line);
		($miRNA_ID,$miRNA,$Coodinates,$Strand,$seedSeq) = split(/\t/,$line);
		# seed seq has to be reverse complement
		$orig = $seedSeq;
		$seedSeq = reverseSeq($seedSeq);
		$m_seesSeq = substr($seedSeq,1);
		
		$miRNA{$miRNA_ID}{mer7A} = $seedSeq;
		$miRNA{$miRNA_ID}{mer8} = $seedSeq."A";
		$miRNA{$miRNA_ID}{mer7_1} = $m_seesSeq."A";
		print OUT "$miRNA_ID\tmer8\t$miRNA{$miRNA_ID}{mer8}\n",
		"$miRNA_ID\tmer7A\t$miRNA{$miRNA_ID}{mer7A}\n",
		"$miRNA_ID\tmer7_1\t$miRNA{$miRNA_ID}{mer7_1}\n";
			
	}
	close(IN);
	close(OUT);
	return();
}
sub load_position{
	my($inFile) = @_;
	open(IN,"$inFile")|| die "there is no variation file\n";
	while($line = <IN>){
		chomp($line);
		($chr,$start,$end,$ref,$alt,$gene,$type) = split(/\t/,$line);
		$chr =~ s/chr//;
		
		##($type eq 'snp') or 
	#	if(($type eq 'ins')){
			$VAR{$chr}{$start}{ref} = $ref;
			$VAR{$chr}{$start}{alt} = $alt;
			$VAR{$chr}{$start}{gene} = $gene;
			$VAR{$chr}{$start}{type} = $type;
	#	}
		
	}
	close(IN);
	return();
}
sub extract_SNPs_from_annovar{
	my($inFile,$output) = @_;
	$selectedLocs=();
	# selectedSNPs File
	open(OUT,">$output")|| die "can not write to out file of selected SNP annotation, quit process\n";
	
	
	open(IN,"$inFile")|| die "no input from annovar, quit process\n";
	# parse header
	$line = <IN>;
	
	$index = 0;
	my(@headerVar) = split(/\t/,$line);
	foreach $item (@headerVar){
		$item =~ s/\./_/g;
		$annotation{$item}=$index;
		$index++;
	}
	
	while($line = <IN>){
		chomp($line);
		my(@data) = split(/\t/,$line);
		$select = 0;
		
		if(($UTR3==1) and ($data[$annotation{Func_refGene}] eq 'UTR3')){
			if($targetScan==1){# to select UTR
				if($data[$annotation{TargetS16_sites}] ne "."){
					$select++;
				}
				if($data[$annotation{TargetS16_nonconsvSites_consvFam}] ne "."){
					$select++;
				}
				if($data[$annotation{TargetS16_noncosvSites_broadConsFam}] ne "."){
					$select++;
				}
				
				if(($select>0) and ($selectMotorN_target >0)){
					# parse TargetScan names
					%mirs=();
					if($data[$annotation{TargetS16_sites}] ne "."){
						my($name,$mir) = split(/=/,$data[$annotation{TargetS16_sites}]);
                        my(@mirs) = split(/,/,$mir);
                        foreach $temp (@mirs){
                                $mirs{$temp}++;
                        }

					}
					if($data[$annotation{TargetS16_nonconsvSites_consvFam}] ne "."){
						my($name,$mir) = split(/=/,$data[$annotation{TargetS16_nonconsvSites_consvFam}]);
                        my(@mirs) = split(/,/,$mir);
                        foreach $temp (@mirs){
                                $mirs{$temp}++;
                        }
					}
					if($data[$annotation{TargetS16_noncosvSites_broadConsFam}] ne "."){
						my($name,$mir) = split(/=/,$data[$annotation{TargetS16_noncosvSites_broadConsFam}]);
                        my(@mirs) = split(/,/,$mir);
                        foreach $temp (@mirs){
                                $mirs{$temp}++;
                        }
					}
				
					$select = 0;
					if($selectMotorN_target == 1){
						foreach $temp (keys %mirs){
							if($human_enriched{$temp} > 0){
								$select++;
							}
						}					
					}elsif($selectMotorN_target == 2){
						foreach $temp (keys %mirs){
							if($from_elik{$temp} > 0){
								$select++;
							}
						}					
					}elsif($selectMotorN_target == 3){
					
						foreach $temp (keys %mirs){
							if($extra_literature{$temp} > 0){
								$select++;
							}
						}
					}elsif($selectMotorN_target == 4){
						
						foreach $temp (keys %mirs){
						
							if($all_motor_neurons{$temp} > 0){
								$select++;
							}
						}	
					}
				}
				
			}else{
				$select++;
			}
		}
		if(($select >0) and ($data[$annotation{Gene_refGene}] !~/\,/)){# we select only autosomal chromosomes
			# determine SNPtype
			$SNPtype = '';
			if((length($data[$annotation{REF}]) ==1) and (length($data[$annotation{ALT}]) ==1)){
				$SNPtype='snp';
			}elsif(length($data[$annotation{ALT}]) >1 ){
				$SNPtype='ins';
			}elsif(length($data[$annotation{REF}]) > 1){
				$SNPtype='del';
			}
			
			# change indels into annovar format
			#if($SNPtype eq'ins'){
			#	$data[$annotation{ALT}]=substr($data[$annotation{ALT}],1);
			#	$data[$annotation{REF}]='-';
			#}elsif($SNPtype eq'del'){
			#	$data[$annotation{REF}]=substr($data[$annotation{REF}],1);
			#	$data[$annotation{ALT}]='-';
			
			#}
			print OUT "$data[$annotation{CHROM}]\t$data[$annotation{POS}]\t$data[$annotation{POS}]\t$data[$annotation{REF}]\t$data[$annotation{ALT}]\t$data[$annotation{Gene_refGene}]\t$SNPtype\n";
		}
	}
	close(IN);
	close(OUT);
	
	return(%selectedLocs);
	}
sub get_snp_seq{
	my($new_seq)= @_;
	$seq = bless($new_seq,Bio::Seq);
	$Sloc = $loc+$S+1;
	$Eloc = $Sloc+8-1;
	#extract sequence

	$seqToCheck[$i] = $seq->subseq($Sloc,$Eloc);
	$OriSeq[$i] = $seqToCheck[$i];

	#calculate the mutated sequence
	$mutPos = $StrL-$i+1;
	substr($seqToCheck[$i],$mutPos-1,1) = $VAR{$chrNum}{$variation}{alt};
	return();
}
sub get_del_seq{
	my($new_seq)= @_;
	$seq = bless($new_seq,Bio::Seq);
	$delL = length($VAR{$chrNum}{$variation}{ref})-1;
	$mutCount = 8-$delL+1;# number of possible sequences around the mutation
	
	#for the moment, we don't deal with dels larger than 1
	$OriSeq[$i] = '';
	$seqToCheck[$i] = '';
	if($delL <2 ){

	
		# we calculate from the 3' of the del
		$Sloc = $loc+$delL+$S+1;
		$Eloc = $Sloc+8-1;
		$OriSeq[$i] = $seq->subseq($Sloc,$Eloc);
		
		#calulate mutated sequence
		if($i<8){
			#preseq
			$nSloc = $Sloc-$delL+1;
			$nEloc = $loc;
			$preSeq = $seq->subseq($nSloc,$nEloc);

			# post seq
			$nSloc = $loc+$delL+1;
			$nEloc = $nSloc+$S+8;
			$postSeq = $seq->subseq($nSloc,$nEloc);
		
			$seqToCheck[$i] ="$preSeq$postSeq";
		}
		
		if($i > $mutCount){
			$OriSeq[$i] = '';
		}

	}
	
	return();
}					
sub get_ins_seq{
	my($new_seq)= @_;
	$seq = bless($new_seq,Bio::Seq);
	$delL = length($VAR{$chrNum}{$variation}{alt})-1;
	$var = substr($VAR{$chrNum}{$variation}{alt},1);
	$OriSeq[$i] = '';
	$seqToCheck[$i] = '';
	if($delL < 2){
		# extract wt
		$Sloc = $loc+$delL+$S;
		$Eloc = $Sloc+8-1;
		$OriSeq[$i] = $seq->subseq($Sloc,$Eloc);
		
		# extract preseq, postseq
		# explanations in excel
		$preSeq = '';
		$postSeq = '';
		
		$preS = $Sloc+1;
		$preE = $loc;
		$postS = $loc+1;
		$postE = $loc+$i-1;
		
		if($i == 1){
			$preSeq = $seq->subseq($preS,$preE);
			$postSeq = '';
		}elsif($i == 8){
			$preSeq = '';
			$postSeq = $seq->subseq($postS,$postE);
		}else{
			$preSeq = $seq->subseq($preS,$preE);
			$postSeq = $seq->subseq($postS,$postE);
		}
		$seqToCheck[$i] = "$preSeq"."$var"."$postSeq";
	}
	return();
}
sub check_wt{

	foreach $miRNA_ID (keys %miRNA){
		$found8{$miRNA_ID}=0;
		$found7A{$miRNA_ID}=0;
		$found7{$miRNA_ID}=0;
		for($i=1; $i<9;$i++){	
			if($OriSeq[$i]=~/$miRNA{$miRNA_ID}{mer8}/ ){
				$found8{$miRNA_ID}++ ;		
			}
			if(($OriSeq7[$i]=~/$miRNA{$miRNA_ID}{mer7_1}/ )){
				$found7A{$miRNA_ID}++;	
			}
			if(($OriSeq7[$i]=~/$miRNA{$miRNA_ID}{mer7A}/ )){
					$found7{$miRNA_ID}++;
			}
			if(($OriSeq7_1[$i]=~/$miRNA{$miRNA_ID}{mer7_1}/ )){
				$found7A{$miRNA_ID}++;
					
			}
			if(($OriSeq7_1[$i]=~/$miRNA{$miRNA_ID}{mer7A}/ )){
					$found7{$miRNA_ID}++;
			}
		}
	}# end check wildtype
	return();
}
sub check_mut{
	my($FH) = $_[0];
	foreach $miRNA_ID (keys %miRNA){
		$flag = 0;
		$score = '';
		for($i=1; $i<9;$i++){

			if(($miRNA{$miRNA_ID}{mer8} eq $seqToCheck[$i]) and($found8{$miRNA_ID} == 0)){
				
				if(($found7A{$miRNA_ID} > 0) or ($found7{$miRNA_ID} > 0)){
					$score = "_s";
				}else{
					$score = '';
				}
				print $FH "$VAR{$chrNum}{$variation}{gene},$chrNum,$variation,$VAR{$chrNum}{$variation}{alt},$VAR{$chrNum}{$variation}{type},$miRNA_ID",",$miRNA{$miRNA_ID}{mer8},$OriSeq[$i],$seqToCheck[$i],8mer$score,$i\n";

				$flag++;
			
			}elsif($miRNA{$miRNA_ID}{mer7A} eq $seqToChek7A[$i]){
				if(($found8{$miRNA_ID} >0) or($found7{$miRNA_ID} >0)){
				}elsif($flag == 0){
					if($found7A{$miRNA_ID} >0){
						$score = "_s";
					}else{
						$score = '';
					}
					print $FH "$VAR{$chrNum}{$variation}{gene},$chrNum,$variation,$VAR{$chrNum}{$variation}{alt},$VAR{$chrNum}{$variation}{type},$miRNA_ID",",$miRNA{$miRNA_ID}{mer7A},$OriSeq7,$seqToChek7A,7mer-m8$score,$i\n";
					$flag++;
				}
				#if(($found8{$miRNA_ID} >0) or($found7{$miRNA_ID} >0)){
			}elsif($miRNA{$miRNA_ID}{mer7_1} eq $seqToChek7A[$i]){
				if(($found8{$miRNA_ID} >0) or ($found7{$miRNA_ID} >0) or($found7A{$miRNA_ID} >0)){
				}elsif($flag==0){
					print $FH "$VAR{$chrNum}{$variation}{gene},$chrNum,$variation,$VAR{$chrNum}{$variation}{alt},$VAR{$chrNum}{$variation}{type},$miRNA_ID,$miRNA{$miRNA_ID}{mer7_1},$OriSeq7,$seqToChek7A,7mer-1A,$i\n";
				}
			}
		}# end screen all mutated options
	}# end screen all mirs

	
	return();
}