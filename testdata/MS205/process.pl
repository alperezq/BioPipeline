#!/usr/bin/perl
$argc = @ARGV;

if ($argc != 1)   
{
    print "usage: \n perl scrummbleTest.pl sequence_file \n";
    exit;
}
my $seq_file = $ARGV[0];


#print "No of iterations: $numberofsimulations \n";

open(seqIN, "<$seq_file")
    or die "Oops! can't open input file ".$seq_file."\n"; 
my @Seqs=<seqIN>;
close(seqIN);


my @unique_types;
my @unique_type_length;
my $tmplen;
# read each sequence and scrumble its unit types

foreach $seq (@Seqs){
    chop $seq;
    my @types=split(/[\t+|\s+]/,$seq);
    $types[0]=~s/\s+//;
    $seqlen=scalar(@types);    
    if($seqlen==2){
	if(!($types[0]=~/a/i)){
	    next;
	}
	if(length($types[0])>1){
	    if($types[1]=~/aacggcact/i){
		print "> HG $types[1]\n";
		my @units=split(//,$types[0]);
		foreach $unit (@units){
		    print "$unit "
		}
		print "\n";
	    }
	}
    }
    
}

exit (0);
# HAP GROUP C gaacggcact
# HAP GROUP E gaacggccct
