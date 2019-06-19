#!/usr/bin/perl
$argc = @ARGV;

if ($argc != 2)   
{
    print "usage: \n perl processDoubleCopy.pl sequence_file costfile\n";
    exit;
}
my $seq_file = $ARGV[0];
my $costfile= $ARGV[1];

#print "No of iterations: $numberofsimulations \n";

open(seqIN, "<$seq_file")
    or die "Oops! can't open input file ".$seq_file."\n"; 
my @Seqs=<seqIN>;
close(seqIN);


## reading cost file
open(INcostfile, "<$costfile")
    or die "Oops! can't open cost file ".$costfile."\n"; 
my @lines=<INcostfile>;
close(INcostfile);

my $originaltypes_no=0;


#my @unique_types;
#my @unique_type_length;
my $tmplen;
my @uniqtypes;
my @modtypes;

foreach $line (@lines){
    if($line=~/Type no\s+(\d+)/){
	$originaltypes_no=$1;
    }
    if($line=~/Types\s+(.+)/){
	$line2=$1;
	#chop $line;
	@uniqtypes=split(/\s+/,$line2);
	#print "Hallo @uniqtypes size: ".scalar(@uniqtypes)."\n";
    }
    
 
}





# read each sequence and scrumble its unit types

my $startblok; # start of a block of tandem array
my $endblok;   # end block of a tandem array

my @finaltypes;

my @blockarrayS;
my @blockarrayE;

my $seqid=0;




foreach $seq (@Seqs){
    if($seq=~/>/){
	print $seq;
	next;
    }
    $seqid++;
    my $totalcollaps=0;
    chop $seq;    
    my @types=split(/[\t+|\s+]/,$seq);        
    $seqlen=scalar(@types);    
    my $startblockflag=0;

    
    for($i=0;$i<$seqlen;$i++){ # check unique types
	my $foundflag=0;
	for($kk=0;$kk<scalar(@uniqtypes);$kk++){
	    if($uniqtypes[$kk] eq $types[$i]){
		$foundflag=1;
		last;
	    }
	}
	if($foundflag==0){
	    #push(@uniqtypes,$types[$i]);
	    print stderr "Error: Unkown type $types[$i]\n";
	    exit(-1);
	}

    } #end check unique types


    for($i=3;$i<$seqlen;$i++){ # loop on the whole sequence
	if(($types[$i] eq $types[$i-2])
	   && ($types[$i-1] eq $types[$i-3])
	   && ($types[$i] ne $types[$i-1])){
	    
	    if($startblockflag==0){
		$startblockflag=1;
		$startblock=$i-3;
		$endblock=$i;		
	    }
	    else{		
		if($i-$endblock==2){
		    $endblock=$i;		 
		}
	    }
	}
	elsif($startblockflag==1){
	    push(@blockarrayS,$startblock);
	    push(@blockarrayE,$endblock);
	    $totalcollaps=$totalcollaps+($endblock-$startblock+1)/2;
#	    ($Mstartblock,$Mendblock)=($blockarrayS[scalar(@blockarrayS)-1],$blockarrayE[scalar(@blockarrayS)-1]);
#	    print "Before BushStartB $Mstartblock EndB $Mendblock\n";
	    $startblockflag=0;
	}
    }
#    print "Done $seqlen \n";
    push(@blockarrayS,$seqlen);
    push(@blockarrayE,$seqlen);

    my $blockcounter=0;
    my $ouputseqlen=0;
    for($j=0;$j<$seqlen;$j++){

	($Lstartblock,$Lendblock)=($blockarrayS[$blockcounter],$blockarrayE[$blockcounter]);
#	print "StartB $Lstartblock EndB $Lendblock\n";
	
	if($j<$Lstartblock){ # j is not within a block
	    print $types[$j]." ";
	    $ouputseqlen++;
	}
	else{ # we are at start of block	
	    
	    $typemod=$types[$j]."$types[$j+1]";	    
	    #$typemod="b";


            # check unique types
	    my $foundflag=0;
	    for($kk=0;$kk<scalar(@modtypes);$kk++){
		if($modtypes[$kk] eq $typemod){
		    $foundflag=1;
		    last;
		}
	    }
	    if($foundflag==0){
		#push(@uniqtypes,$typemod);
		push(@modtypes,$typemod);
	    }
		
	    #end check unique types



	    for($k=0;$k<($Lendblock-$Lstartblock);$k=$k+2){
		print $typemod." ";
		$ouputseqlen++;
	    }	    
	    $j=$Lendblock;
	    $blockcounter++;
	   
	}
	
    }
    ## quality check
    print "\n";

    print stderr "Seq: $seqid Outputlengh:  $ouputseqlen CollapseLen: $totalcollaps OriginlLen: $seqlen \n";
    if(($ouputseqlen+($totalcollaps))!=$seqlen){
	print "\n Error in  transforming sequence $seqid \n"; 
	exit;
	
    }
    $tmplen=scalar(@blockarrayS);
    for($i=0;$i<$tmplen;$i++){
	pop(@blockarrayS);
	pop(@blockarrayE);	
	
    }
    for($i=0;$i<$seqlen;$i++){
	pop(@types);
	
    }

  

}

print stderr "Done successful\n";
print stderr "Unique Types in final sequences: ";
for($kk=0;$kk<scalar(@uniqtypes);$kk++){
    print stderr "$uniqtypes[$kk] ";
    
}
print stderr "\n";


#########
# Handling cost file

open(outcostfile, ">$costfile.mod")
    or die "Oops! can't open input file $costfile.mod \n"; 

my @distmatrix;

$typeno= scalar(@uniqtypes)+ scalar(@modtypes);
$matrixline=0;
$linecount=0;
foreach $line (@lines){
    
    if($line=~/Type no./){	
	print outcostfile "# Type no. $typeno \n";       
    }
    elsif($line=~/Types/){
	print outcostfile "# Types ";
	for($kk=0;$kk<scalar(@uniqtypes);$kk++){
	    print outcostfile "$uniqtypes[$kk] ";    
	}
	for($kk=0;$kk<scalar(@modtypes);$kk++){
	    print outcostfile "$modtypes[$kk] ";
	}
	print outcostfile "\n";
    }
    elsif($line=~/Indel/){
	
	print outcostfile $line;
    }
    elsif($line=~/Dup/){
	print outcostfile $line;
    }
    elsif($line=~/matrix/){	
	print outcostfile $line; 
	$linecount++;
	$matrixline=0;
	for($kk=0;$kk<(scalar(@uniqtypes)-1);$kk++){	    
	    my @tokens=split(/\s+/,$lines[$linecount]);
	    @distrow;	    
	    push(@distmatrix,\@tokens);
	    # for debug
   	    #for($jj=0;$jj<scalar(@tokens);$jj++){
	    #	print outcostfile "$distmatrix[$matrixline][$jj] ";		
	    #}
	    # end for debug
	    $linecount++;
	    #print outcostfile "\n";
	    $matrixline++;
	    
	}	
	last;
    }
    $linecount++;
}
my $rowno=0;
my $colno=0;


foreach my $row(@distmatrix){
   foreach my $val(@$row){
      print outcostfile "$val ";
   }
   foreach my $modtoken(@modtypes){               
       $firsthalf=substr($modtoken,0,(length($modtoken)/2));
       $secondhalf=substr($modtoken,(length($modtoken)/2));
       #locate first & second half
       for($kk=0;$kk<scalar(@uniqtypes);$kk++){
	   if($uniqtypes[$kk] eq $firsthalf){
	       $indexFH=$kk;	      
	   }
	   if($uniqtypes[$kk] eq $secondhalf){
	       $indexSH=$kk;	      
	   }
       }
       # compute mutation between first and second half, 1 is added for dup
       if($indexFH<$indexSH){$mutationcost=1+$distmatrix[$indexFH][$indexSH-$indexFH-1];}
       else{$mutationcost=1+$distmatrix[$indexSH][$indexFH-$indexSH-1];}
       
       if($rowno==$indexFH){$dist1=0;}
       else
       {
	   if($rowno<$indexFH){
	       $dist1=$distmatrix[$rowno][$indexFH-$rowno-1];
	   }
	   else{
	       $dist1=$distmatrix[$indexFH][$rowno-$indexFH-1];
	   }
       }
       if($rowno==$indexSH){$dist2=0;}
       else
       {
	   if($rowno<$indexSH){
	       $dist2=$distmatrix[$rowno][$indexSH-$rowno-1];
	   }
	   else{
	       $dist2=$distmatrix[$indexSH][$rowno-$indexSH-1];
	   }
       }
 #      print "HLLO FH: $indexFH SH: $indexSH D1: $dist1 D2: $dist2 \n";
       if($dist1<$dist2){print outcostfile ($dist1+$mutationcost)." ";}
       else{print outcostfile ($dist2+$mutationcost)." ";}
       
   }
   $rowno++;
   print outcostfile"\n";
}

## printing distance to final unit
$flag=0;
foreach my $modtoken(@modtypes){               
    $flag=1;
    $firsthalf=substr($modtoken,0,(length($modtoken)/2));
    $secondhalf=substr($modtoken,(length($modtoken)/2));
    #locate first & second half
    for($kk=0;$kk<scalar(@uniqtypes);$kk++){
	if($uniqtypes[$kk] eq $firsthalf){
	    $indexFH=$kk;	      
	}
	if($uniqtypes[$kk] eq $secondhalf){
	    $indexSH=$kk;	      
	}
    }
    # compute mutation between first and second half, 1 is added for dup
    if($indexFH<$indexSH){$mutationcost=1+$distmatrix[$indexFH][$indexSH-$indexFH-1];}
    else{$mutationcost=1+$distmatrix[$indexSH][$indexFH-$indexSH-1];}
    $dist1=$distmatrix[$indexFH][scalar(@uniqtypes)-2-$indexFH];
    $dist2=$distmatrix[$indexSH][scalar(@uniqtypes)-2-$indexSH];
    if($dist1<$dist2){print outcostfile ($dist1+$mutationcost)." ";}
    else{print outcostfile ($dist2+$mutationcost)." ";}

}
if($flag==1){
    print outcostfile"\n";
}

## printing distances among modified units
for($i=0; $i<scalar(@modtypes);$i++){    
    $modtoken=$modtypes[$i];
    $firsthalf=substr($modtoken,0,(length($modtoken)/2));
    $secondhalf=substr($modtoken,(length($modtoken)/2));
    #locate first & second half
    for($kk=0;$kk<scalar(@uniqtypes);$kk++){
	if($uniqtypes[$kk] eq $firsthalf){
	    $indexFH=$kk;	      
	}
	if($uniqtypes[$kk] eq $secondhalf){
	    $indexSH=$kk;	      
	}
    }
    $flag=0; # to check if entered the loop
    for($j=$i+1; $j<scalar(@modtypes);$j++){  
	$flag=1;
  	$modtoken2=$modtypes[$j];
	$firsthalf2=substr($modtoken2,0,(length($modtoken2)/2));
	$secondhalf2=substr($modtoken2,(length($modtoken2)/2));
	#locate first & second half
	for($kk=0;$kk<scalar(@uniqtypes);$kk++){
	    if($uniqtypes[$kk] eq $firsthalf2){
		$indexFH2=$kk;	      
	    }
	    if($uniqtypes[$kk] eq $secondhalf2){
		$indexSH2=$kk;	      
	    }
	}

	# compute mutation between first and second half, 1 is added for dup
	if($indexFH<$indexFH2)
	{$dist1=$distmatrix[$indexFH][$indexFH2-$indexFH-1];}
	else{$dist1=$distmatrix[$indexFH2][$indexFH-$indexFH2-1];}

	if($indexSH<$indexSH2)
	{$dist2=$distmatrix[$indexSH][$indexSH2-$indexSH-1];}
	else{$dist2=$distmatrix[$indexSH2][$indexSH-$indexSH2-1];}
	print outcostfile ($dist1+$dist2+1)." ";
    }
    if($flag==1){print outcostfile"\n"};  
}
     



close(outcostfile);



exit (0);
# HAP GROUP C gaacggcact
# HAP GROUP E gaacggccct
