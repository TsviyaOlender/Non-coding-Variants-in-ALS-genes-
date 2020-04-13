#!
$flag = 1;
$file1 = "linc_QC3_Jan132020_annot.txt.1";# original annotation file
$file2 = "linc_QC3_Jan132020_annot.txt.2";# filtration after step 1, call it something like cohort_annot.txt.2
$file3 = "linc_QC3_Jan132020_annot.txt.3";# filtration after step 2, final! call it something like cohort_annot.txt.3
if($flag == 1){
open(IN,"$file1")|| warn "can not find input file\n";
open(OUT,">$file2")|| warn "can not write to outputfile\n";
$line = <IN>;
print OUT "$line";

$countVars = 0;
$count_gains_1=0;
$count_gains_2=0;

while($line = <IN>){
  chomp($line);
  (@data) = split(/\t/,$line);
  # we analyze $data[42]- the output of 3UTR gain 
  #miR-548o-3p:7mer-1A;miR-1323:7mer-1A;miR-4802-3p:7mer-m8;
  (@gains) = split(/\;/,$data[42]);
  $new_prediction = '';
   
   if($data[42] =~/miR/){
     $count_gains_1++;
   }
   #filter
   $new_prediction = '';
  foreach $prediction (@gains){
    ($mir,$type) = split(/\:/,$prediction);
 #   if($type =~/^8mer/){
      $temp = $mir.":".$type.";";
      $new_prediction .=$temp;
  #  }
    
  }
  $data[42] = $new_prediction;
  if($data[42] =~/miR/){
     $count_gains_2++;
   }
   $newLine = '';
   foreach $annot (@data){
     $newLine .= "$annot\t";
   }
   print OUT "$newLine\n";
   $countVars++;
}
close(OUT);
close(IN);

print "varinat count = $countVars, With all 3UTR gains =$count_gains_1, only with 8-mers =$count_gains_2\n";
}

# step2
# count occurences
open(IN,"$file2")|| warn "can not find input file\n";
$line = <IN>;
while($line = <IN>){
  chomp($line);
  (@data) = split(/\t/,$line);
  $key = $data[0]."_".$data[1];
  $countMulti{$key}++;
}
close(IN);
#filter
open(IN,"$file2")|| warn "can not find input file\n";
open(OUT,">$file3")|| warn "can not write to outputfile\n";
$line = <IN>;
print OUT "$line";
while($line = <IN>){
  chomp($line);
  (@data) = split(/\t/,$line);
  $key = $data[0]."_".$data[1];
  if($countMulti{$key} > 1){
    #print "$key,$countMulti{$key}\n";
  }else{
    print OUT "$line\n";
  }
}
close(IN);
close(OUT);

system("rm $file2");

