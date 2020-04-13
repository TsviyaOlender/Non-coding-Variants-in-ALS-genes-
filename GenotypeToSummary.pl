#!
my($vcf) = @ARGV;
chomp($vcf);

open(IN,"$vcf")|| die "can not locate vcf file\n";
open(OUT,">genotypes_report.txt")|| die "can not write to out file\n";
do{
	$line = <IN>;
}until($line =~/^#CHROM/);

chomp($line);
my($CHROM,$POS,$ID,$REF,$ALT,$QUAL,$FILTER,$INFO,$FORMAT,@names) = split(/\t/,$line);
print OUT "SNP\tIndividuals\n";
while($line = <IN>){
	chomp($line);
	my($CHROM,$POS,$ID,$REF,$ALT,$QUAL,$FILTER,$INFO,$FORMAT,@geno) = split(/\t/,$line);
	$genotypes_ind = '';
	$i=0;
	foreach $ind (@geno){
	    (@data) = split(/\:/,$ind);
		if(($data[0] !~ /0\/0/) and ($data[0] !~ /\.\/\./)){
			$genotypes_ind .= "$names[$i]".";";

		}
		$i++;
	}
	print OUT "$CHROM\t$POS\t$REF\t$ALT\t$genotypes_ind\n";
}
close(IN);
close(OUT);
