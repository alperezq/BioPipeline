#!/usr/bin/perl
$argc = @ARGV;

if ($argc != 3)   
{
    print "usage: \n perl plotStatistics.pl comparisonfile improved_left_real improved_real_right\n";
    exit;
}
my $seq_file = $ARGV[0];
my $improvedLeft=$ARGV[1];
my $improvedRight=$ARGV[2];

open(seqIN, "<$seq_file")
    or die "Oops! can't open input file ".$seq_file."\n"; 
my @Seqs=<seqIN>;
close(seqIN);


my @values;

# read each sequence and scrumble its unit types
foreach $seq (@Seqs){
    if($seq=~/L: (\d+), R: (\d+), diff: (-*\d+)/){
	#print $seq;
	#print "Line: ".$1." ".$2." ".$3." \n";
	if(($1+$2)>0){
	    $value=$3/($1+$2);
	    #if(abs($value)< 0.0001){$value=0;}
	}
	#print "Line: ".$1." ".$2." ".$3." ".$value." \n";
	push(@values,$value);
    }    

}


@valuesSorted= sort { $a <=> $b } @values;
$arraysize=scalar(@valuesSorted);

$max=$valuesSorted[$arraysize-1];
$min=$valuesSorted[0];

$no_interval=10;

$interval=($max-$min)/$no_interval;
$step=$interval/$no_interval;

my $counter=0;
my $sum=$min+$interval;

$outputfile=$seq_file.".dat";
open(seqOut, ">$outputfile")
    or die "Oops! can't open file ".$outputfile."\n";

print seqOut "#--- $max, $min, $interval: ImprovedLeft: $improvedLeft, ImprovedRight: $improvedRight\n";

$maxCounter=0;

$vc=0;

print "Runs no: ".$arraysize." Min Score: ".$min." Max Score: ".$max."\n";
print "ImprovedLeft: $improvedLeft, ImprovedRight: $improvedRight\n";

$counter=0;
for($ic=0;$ic<$no_interval;$ic++){
    while(($valuesSorted[$vc]<$sum)&&($vc<$arraysize)){
	$counter++;
	$vc++;
    }
    print "Interval: ".$sum." : Conuter: ".$counter."\n";
    print seqOut $sum." ".$counter."\n";
    $sum=$sum+$interval;
    if($counter>$maxCounter){$maxCounter=$counter;}
    $counter=0;
}

#foreach $value (@valuesSorted){
#    if(($value<=($sum+$interval))&&($value>$sum))
#    {	$counter++;
#	if($counter>$maxCounter){$maxCounter=$counter;}
#    }
#    else if ($value>($sum+$interval){
#	print "Interval: ".$sum." : Conuter: ".$counter."\n";
#	print seqOut $sum." ".$counter."\n";
#	$sum=$sum+$interval;
#	$counter=1;
#    }    
#}

#if($counter!=1){
#    print "Interval: ".$sum." : Conuter: ".$counter."\n";
#    print seqOut $sum." ".$counter."\n";
#}

if(($improvedLeft+$improvedRight)!=0){
    $value=($improvedLeft-$improvedRight)/($improvedLeft+$improvedRight);
}
else{
    $value=0;
}
print seqOut "\n".$value." ".($maxCounter)."\n"; 


close(seqOut);

$maxCounter2=$maxCounter*(1.25);

$plotfile=$outputfile.".gp";
open(seqOut, ">$plotfile")
    or die "Oops! can't open input file ".$plotfile."\n"; ;

print seqOut "set terminal postscript landscape color \n";
print seqOut "set out "."'$plotfile.ps'"." \n";
print seqOut "set nokey \n";
print seqOut "set nogrid \n";
print seqOut "set xlabel \"X\" \n";
print seqOut "set ylabel \"Count\" \n";
print seqOut "set xrange [-1.2:1.2] \n";
print seqOut "set yrange [0:$maxCounter2] \n";
print seqOut "plot \"$outputfile\" index 0 w boxes \n"; 
close(seqOut);

$argstring="gnuplot $plotfile";
print "Runing: $argstring.\n";
$retval=system($argstring);


$argstring="gv $plotfile.ps";
print "Runing: $argstring.\n";
system($argstring);
if($retval!=0){
    print "Error running: $argstring\n";
    exit -1;
} 

exit;
