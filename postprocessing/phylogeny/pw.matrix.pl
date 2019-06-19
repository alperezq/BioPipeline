
#!/usr/bin/perl
$argc = @ARGV;

if ($argc != 2)   
{
    print "usage: \n perl  pw.matrix.pl ARLEMOutfile SequenceFile\n";
    exit;
}
my $infile = $ARGV[0];
my $seq_file = $ARGV[1];

my $imgSize;

if($argc>1){
#    $imgSize=$ARGV[3];
}
else{
#    $imgSize=256;
}

open(IN, "<$infile") or die "Unable to opern $infile\n";
my @lines=<IN>;
close(IN);

open(seqIN, "<$seq_file") or die "Unable to open $seq_file\n";
my @preSeqs=<seqIN>;
close(seqIN);

$outfile=$infile.".mat";
$treefile=$infile.".nwk";

open(distfile, ">$outfile") or die "Unable to open $outfile\n";


my @Seqs;
my $Seqs_names;
foreach $Seq (@preSeqs){
    if($Seq=~/>/){
	push (@Seqs_names,$Seq);
	next;
    }
    push (@Seqs,$Seq);
}
foreach $Seq (@Seqs){
    if($Seq=~/>/){	
	next;
    }
}

my @score_array;

my $seq_no=scalar @Seqs;

print distfile $seq_no."\n";

#exit; 
my $stringS="";
my $stringR="";





### part 1 read alignment
my $flag=0;

#my $counterS=0;
#my $counterR=0;
my $line_beg=0;
my $line_end=0;
my $comp_counter=0;

$line_n=0;
foreach $line2 (@lines){
    if($line2=~/operations/){
	$flag=1;
	$line_beg=$line_n;
	#next;
    }
    if($line2=~/Score of aligning/){
	$line2=~/Seq:(\d+).+Seq:(\d+)\s+=(\d+.*[\d+]*)/;
	$seq1_no=$1;
	$seq2_no=$2;
	$score=$3;
	$flag=1;
	$line_end=$line_n;	
	$stringS=$Seqs[$seq1_no];
	$stringR=$Seqs[$seq2_no];
	$comp_counter=$seq1_no."x".$seq2_no;
	#draw_alignment($line_beg,$line_end,$comp_counter);
#	print "Hallo ".$seq1_no." ".$seq2_no." ".$score."\n";
	push(@score_array,$score);
	$comp_counter++;
	#next;
    }
    
    $line_n++;
}
#print "HIIIIIIII Seq No. ".$seq_no."\n";
for($i=0;$i<$seq_no;$i++){
    if($Seqs_names[$i]){
	$Seqs_names[$i]=~s/\>//g;
	print distfile $Seqs_names[$i];
    }
    else{
	print distfile"S$i \n";
    }
    for($j=0;$j<$seq_no;$j++){
	
	if($i==$j){
	    print distfile "0\t";
	}
	else{
	    $i2=$i; $j2=$j;
	    if($i>$j){
		$tmp=$i2;
		$i2=$j2;
		$j2=$tmp;
	    }
	    $offset=($i2*(2*($seq_no-1)-$i2+1))/2;
	    $index=$offset+abs($j2-$i2)-1;
	    print distfile $score_array[$index]."\t";
	}
    }
    print distfile "\n";
}

close (distfile);

$argstring="./bionj_linux.x $outfile $treefile";
print $argstring."\n";
system($argstring);

$pages=$seq_no/60;
if($pages<1){$pages=1;}

if($seq_no<50){
    $fontsize=10;
}
elsif($seq_no<80){
    $fontsize=8;
}
else{
    $fontsize=6;
}

#$argstring="njplot.x -pc $pages -size $fontsize -psize 800x1040 $treefile";
$argstring="./njplot.x -pc $pages -size $fontsize  $treefile";
system($argstring);


$treefile=~/(.+\/)*(.+)/;
$filepath=$1;
$filename=$2;
if(!($filepath)){$filepath="./";}
print "File path: ".$filepath."\n";
print "File name: ".$filename."\n";

$filename=~/(.*?)\./;
$prefix=$1;


print "Prefix $prefix \n";


$argstring="mv $filepath$prefix.pdf $treefile.pdf";
print $argstring."\n";
system($argstring);


print "Ouput files: \n";
print "$outfile \n";
print "$treefile \n";
print "$treefile.pdf \n";
exit;
####





