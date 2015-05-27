#!/usr/bin/perl 
#use warnings;

sub sequenceGeneration
{
	my $randSeqLength = shift; 
			

	my @chars=('A','C','G','T');
	my $random_seq;
	
	foreach (1..$randSeqLength) 
	{
		
		$random_seq.=$chars[rand @chars];
		
	}
	
	
	return $random_seq;
}

sub motifGeneration
{
	my $randmotifLength = shift; 
			

	my @chars=('A','C','G','T');
	my $random_motif;
	
	foreach (1..$randmotifLength) 
	{
		
		$random_motif.=$chars[rand @chars];
		
	}
	
	
	return $random_motif;
}

my $noOfDatasets = 50;

for($w= 0; $w<$noOfDatasets;++$w){

$numofSeq = 20;  # NUMBER OF SEQUENCES
$Seqlength = 600; #SEQUENCE LENGTH
my @sequences;

for($a= 0; $a<$numofSeq;++$a)
{
	
	$sequences[$a]= & sequenceGeneration($Seqlength);
	
}


$motifLength = 10;   # LENGTH OF MOTIF
$numoferrors = 2;    # NUMBER OF PLANTED ERRORS

$randMotif = &motifGeneration($motifLength);
my @intervall = (0..($motifLength-1));

my @motive;

push(@motive, $randMotif);
my $originalMotif = $randMotif;

print $originalMotif."\n";

for($o= 0; $o<$numofSeq; ++$o){
	


    my @motifPos;
	my %alreadyTaken = ();
	my $used;
	my @matrix;




for($t=0; $t<$numoferrors; ++$t){
	
	$used = int(rand(scalar(@intervall)));
    

     redo if exists $alreadyTaken{$used};
     $alreadyTaken{$used}++;
     $motifPos[$t]= $used;
	
	
if(substr($randMotif, $motifPos[$t],1) eq 'A'){
		
		
		 @letters=('C','G','T');
		 $total=@letters;
		 $errbase = $letters[rand $total];
		 
		 substr($randMotif, $motifPos[$t],1)= $errbase ;
		
		
		 
	 }
	 elsif(substr($randMotif, $motifPos[$t],1) eq 'C'){
		 
		 
		 @letters=('A','G','T');
		 $total=@letters;
		 $errbase = $letters[rand $total];
		 
		 substr($randMotif, $motifPos[$t],1) = $errbase;
		
		
	 }
	 elsif(substr($randMotif, $motifPos[$t],1) eq 'G'){
		 
		  
		 @letters=('A','C','T');
		 $total=@letters;
		 $errbase = $letters[rand $total];
		 
		 substr($randMotif, $motifPos[$t],1) = $errbase;
		 
		
	 }
	 elsif(substr($randMotif, $motifPos[$t],1) eq 'T'){
		 
		
		
		 @letters=('A','C','G');
		 $total=@letters;
		 $errbase = $letters[rand $total];
		
		 substr($randMotif, $motifPos[$t],1) = $errbase;
		
		
		  
	 }
	
	
}	
    
    push(@motive, $randMotif);
    
    
    $randMotif = $originalMotif
   
} 
    
 
open(DATEI, ">L".$motifLength."D".$numoferrors.$w.".fasta"); 
open(DATEI2,">L".$motifLength."D".$numoferrors.$w."INFO.fasta");

$motive = @motive;
$sequen = @sequences;


 
for($c=0;$c<$numofSeq;++$c){
	
	
	$motifposi = int(rand($Seqlength-$motifLength));
	
	print DATEI ">seq".$c."\n";
	substr($sequences[$c], $motifposi,$motifLength) = $motive[$c];
	
	print DATEI $sequences[$c]."\n";
	
	print DATEI2 ">seq".$c."\n";
	$motif= $motive[$c];
	$range= $motifposi+$motifLength;
	
	print DATEI2 $motif."\t".$motifposi."-".$range."\n";
	
	
	}
    


close(DATEI);
close(DATEI2);

}



