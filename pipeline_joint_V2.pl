#!
#!
# input => vcf file. must be tabix zipped and indexed. Otherwise= it is done by script
# the directory must contain the files: cases_id.txt, ctrl_id.txt
# Step 1: prepere vcf and counting
# step 2: annotation vcf
# step 3: extract table from vcf file
############################################
# module to load
# vcftools, bcftools, tabix=> allready in .bashrc file
$genome = "human_g1k_v37.fasta";
$annovarDir = "annovar";
$genomeV = "hg19";


my($VCF,$cohortName) = @ARGV;
chomp($VCF); 
chomp($cohortName);
$workVCF = $cohortName.".vcf.gz";
$annVCF=$workVCF.".hg19_multianno.txt";
$finalAnnotation =$cohortName."_annot.txt"; 
$VCF_ann = $cohortName.".vcf.gz.hg19_multianno.vcf";
#####
norm_vcf($VCF);# extract bi-allelic SNPs and normalize VCF file
count_allele_freq($workVCF);
run_annotation_pipeLine($workVCF);# run annovar
(%cases_counts)=load_counts("cases_counts.frq.count");
(%ctrl_counts)=load_counts("ctrl_counts.frq.count");
#######
vcf_toTable($VCF_ann);# extract annotation from VCF file. output table1.txt
### run prediction of 3UTRgain
$cmd = "perl $annovarDir/3UTR_gain/calculate_3UTR_impact.V2.pl $cohortName $genome";
system("$cmd");

####################################################################################
sub load_counts{
	# the subroutine loads counts of vcftools to memory. returns hash with counts
	my($infile) = @_;
	%counts= ();
	open(IN,"$infile")|| warn "can not read from file $infile\n";
	$line = <IN>;
	while($line=<IN>){
		chomp($line);
		my($CHROM,$POS,$N_ALLELES,$N_CHR,$a1,$a2)=split(/\t/,$line);
		$key = "$CHROM"."_"."$POS";
		($allele1,$count1) = split(/:/,$a1);
		($allele2,$count2) = split(/:/,$a2);
		$counts{$key}{$allele1}= $count1;
		$counts{$key}{$allele2}= $count2;
	}
	close(IN);
	return(%counts);
}
sub vcf_toTable{
	my($infile) = @_;
	local $" = "\t";
	
	$cmd="vcftools --vcf $infile ";
	$cmd .="--get-INFO Gene.refGene ";
	$cmd .="--get-INFO Func.refGene ",
	$cmd .="--get-INFO ExonicFunc.refGene ";
	$cmd .="--get-INFO AAChange.refGene ";
	$cmd .="--get-INFO avsnp144 ";
	$cmd .="--get-INFO 1000g2014oct_eur ";
	$cmd .="--get-INFO esp6500siv2_ea ";
	$cmd .="--get-INFO cg69 ";
	$cmd .="--get-INFO eXome_db ";
	$cmd .="--get-INFO phastConsElements46way ";
	$cmd .="--get-INFO genomicSuperDups ";
	$cmd .="--get-INFO SIFT_pred ";
	$cmd .="--get-INFO Polyphen2_HVAR_pred ";
	$cmd .="--get-INFO LRT_pred ";
	$cmd .="--get-INFO MutationTaster_pred ";
	$cmd .="--get-INFO MutationAssessor_pred ";
	$cmd .="--get-INFO FATHMM_pred ";
	$cmd .="--get-INFO MetaLR_pred ";
	$cmd .="--get-INFO GERP++_RS ";
	$cmd .="--get-INFO clinvar_20150629 ";
	$cmd .="--get-INFO ExAC_ALL ";
	$cmd .="--get-INFO targetScanS ";
	$cmd .="--get-INFO TargetS16_sites ";#
	$cmd .="--get-INFO TargetS16_conSites_poorConsFam ";#
	$cmd .="--get-INFO TargetS16_nonconsvSites_consvFam ";#
	$cmd .="--get-INFO TargetS16_noncosvSites_broadConsFam ";#
	$cmd .="--get-INFO preMiRNA ";
	$cmd .="--get-INFO matureMiRNA ";
	$cmd .="--get-INFO seedMiRNA ";
	$cmd .="--get-INFO wgEncodeRegDnaseClustered ";
	$cmd .="--get-INFO wgEncodeRegTfbsClustered ";
  $cmd .="--get-INFO GenCode3UTR ";
  $cmd .="--get-INFO PITA_ALL_f ";
  $cmd .="--get-INFO PITA_top_f";

	print "extacting data from vcf\n";
 print "$cmd\n";
	system("$cmd");
	
	# parsing out.INFO
	# final table is $finalAnnotation
	print "making final table. Data will be save in file $finalAnnotation\n";
	open(IN,"out.INFO")||warn "no out.INFO file\n";
	open(OUT,">$finalAnnotation")|| warn "can not write to final output\n";
	$line = <IN>;
	chomp($line);
	($chr,$loc,$ref,$alt,@header) = split(/\t/,$line);
	print OUT "$chr\t$loc\t$ref\t$alt\tcasesRefCount\tcasesAltcount\tctrlRefCount\tctrlAltCount\t@header\n";
	
	while($line = <IN>){
		chomp($line);
		($chr,$loc,$ref,$alt,@data) = split(/\t/,$line);
		$key = "$chr"."_"."$loc";
		print OUT "$chr\t$loc\t$ref\t$alt\t",
			"$cases_counts{$key}{$ref}\t$cases_counts{$key}{$alt}\t",
			"$ctrl_counts{$key}{$ref}\t$ctrl_counts{$key}{$alt}\t";
		foreach $item (@data){
			if($item =~/\\x3d|\\x3b/){
				$item =~ s/\\x3d/=/g;
				$item =~ s/\\x3b/;/g;
				print OUT "$item\t";
			}else{
				print OUT "$item\t";
			}
			
		}
		print OUT "\n";
	}
	close(IN);
	close(OUT);
	
	# clean
	system("rm out.INFO");
	return();
}
sub run_annotation_pipeLine{
        my($infile) = @_;

        $protocol = "refGene,avsnp144,1000g2014oct_eur,esp6500siv2_ea,cg69,phastConsElements46way,genomicSuperDups,dbnsfp30a,wgEncodeRegDnaseClustered,wgEncodeRegTfbsClustered,";
        $protocol .= "clinvar_20150629,exac03,targetScanS,preMiRNA,matureMiRNA,seedMiRNA,eXome_db,TargetS16_conSites_poorConsFam,TargetS16_nonconsvSites_consvFam,";
		$protocol .="TargetS16_noncosvSites_broadConsFam,TargetS16_sites,GenCode3UTR,PITA_Ts_run,PITA_top_f";

        $operation = "g,f,f,f,f,r,r,f,r,r,f,f,r,r,r,r,f,r,r,r,r,r,r,r";
        $args = "\'-splicing $splice_distance\',\'\',\'\',\'\',\'\',\'\',\'\',\'\',\'\',\'\',\'\',\'\',\'\',\'\',\'\',\'\'";
        $command = "$annovarDir/table_annovar.pl -remove --vcfinput $infile $annovarDir/humandb/ -buildver $genomeV -protocol $protocol -operation $operation";

         print "$command\n";
         system("$command");

        return();

}

sub count_allele_freq{
	my($inVCF) = @_;#vcf is expected to be bgzipped and index
	print "Analysis Step2: counts of cases vars and ctrl vars\n#############################################################\n";
	unless (-e "cases_id.txt"){
		print "Can not locat the file cases_id.txt!!!!\n";
	#	exit;
	}
	unless (-e "ctrl_id.txt"){
		print "Can not locat the file ctrl_id.txt!!!!\n";
	#	exit;
	}
	# generate cases vcf file
	system("vcftools --gzvcf $inVCF --out cases --keep cases_id.txt --recode");
	system("bgzip -f cases.recode.vcf");
	system("tabix -f cases.recode.vcf.gz");
	
	# count calls
	system("vcftools --gzvcf cases.recode.vcf.gz --out cases_counts --counts");
	
	# generate ctrl vcf file
	system("vcftools --gzvcf $inVCF --out ctrl --keep ctrl_id.txt --recode");
	system("bgzip -f ctrl.recode.vcf");
	system("tabix -f ctrl.recode.vcf.gz");
	
	# count calls
	system("vcftools --gzvcf ctrl.recode.vcf.gz --out ctrl_counts --counts");
	
	return();
}
sub norm_vcf{
	print "Analysis Step1: VCF normalization\n##################################\n";
	$outName = $cohortName.".vcf";
	if($VCF !~ /gz$/){
		print "bgzip -f $VCF\n";
		system("bgzip -f $VCF");
		$VCF .= ".gz";
		print "tabix -f $VCF\n";
		system("tabix -f $VCF");
	}
	
	# normalize VCF
	print "bcftools norm -f $genome -c s -o temp.vcf $VCF\n";
	system("bcftools norm -f $genome -c s -o temp.vcf $VCF");

	# extract biallelic SNPs
	print "vcftools --vcf temp.vcf --min-alleles 2 --max-alleles 2 --recode --out temp1\n";
	system("vcftools --vcf temp.vcf --min-alleles 2 --max-alleles 2 --recode --out temp1");
	system("mv temp1.recode.vcf $outName\n");
	
	# bgzip and index result
	print "bgzip -f $outName\n";
	system("bgzip -f $outName");
	$outName .= ".gz";
	print "tabix -f $outName\n";
	system("tabix -f $outName");
	
	# clean
	system("rm temp.vcf.gz");
	system("rm temp.vcf.gz.tbi");
	system("rm temp1.recode.vcf");
	
	return();
}
