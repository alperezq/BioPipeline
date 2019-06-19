#!/usr/bin/perl
$argc = @ARGV;

if ($argc != 2)   
{
    print "usage: \n perl plotbreakeven.pl arlem_output_file_real arlem_output_file_random\n";
    exit;
}
my $seq_file1 = $ARGV[0];
my $seq_file2 = $ARGV[1];



$outputfile1=$seq_file1."BE.dat";

$argstring="breakeven.pl $seq_file1 $outputfile1";
$retcode=system($argstring);
if($retcode ne 0){
    print STDERR "\n failure: $argstring\n";
    exit -1;
}

$outputfile2=$seq_file2."BE.dat";
$argstring="breakeven.pl $seq_file2 $outputfile2";
$retcode=system($argstring);
if($retcode ne 0){
    print STDERR "\n failure: $argstring\n";
    exit -1;
}

# draw histogram of I



$plotfile=$outputfile1.".gp";

open(seqOut, ">$plotfile");
print seqOut "set terminal postscript landscape color \n";
print seqOut "set out "."'$plotfile.ps'"." \n";
print seqOut "set nokey \n";
print seqOut "set nogrid \n";
print seqOut "set xlabel \"X\" \n";
print seqOut "set ylabel \"Count\" \n";
print seqOut "set xrange [-1.2:1.2] \n";
print seqOut "set yrange [0:*] \n";
print seqOut "plot \"$outputfile1\" index 0 w boxes, \"$outputfile2\" index 0 w boxes \n"; 
close(seqOut);

$argstring="gnuplot $plotfile";
print "Runing: $argstring.\n";
$retcode=system($argstring);
if($retcode ne 0){
    print STDERR "\n failure: $argstring\n";
    exit -1;
}

$argstring="gv $plotfile.ps &";
print "Runing: $argstring.\n";
system($argstring);


exit;
