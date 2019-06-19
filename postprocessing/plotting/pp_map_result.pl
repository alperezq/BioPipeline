
#!/usr/bin/perl
$argc = @ARGV;

if ($argc != 3)   
{
    print "usage: \n perl  pp_map_result.pl MSATcompareOutfile SequenceFile outputfile_prefix\n";
    exit;
}
my $infile = $ARGV[0];
my $seq_file = $ARGV[1];
my $output_tex = $ARGV[2];

my $imgSize;

if($argc>1){
#    $imgSize=$ARGV[3];
}
else{
#    $imgSize=256;
}

open(IN, "<$infile");
my @lines=<IN>;
close(IN);

open(seqIN, "<$seq_file");
my @preSeqs=<seqIN>;
close(seqIN);
my @Seqs;
foreach $Seq (@preSeqs){
    if($Seq=~/>/){	
	next;
    }
    push (@Seqs,$Seq);
}
foreach $Seq (@Seqs){
    if($Seq=~/>/){	
	next;
    }
}

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
	$line2=~/Seq:(\d+).+Seq:(\d+)/;
	$seq1_no=$1;
	$seq2_no=$2;
	
	$flag=1;
	$line_end=$line_n;	
	$stringS=$Seqs[$seq1_no];
	$stringR=$Seqs[$seq2_no];

	$comp_counter=$seq1_no."x".$seq2_no;
	draw_alignment($line_beg,$line_end,$comp_counter);
	$comp_counter++;
	#next;
    }
    
    $line_n++;
}

exit;
####

sub draw_alignment{
    $begin=$_[0];
    $end=$_[1];
    $outputfilecounter=$_[2];

    my @seqS=split(" ",$stringS);
    my @seqR=split(" ",$stringR);	

    my $AliS;
    my $AliR;
    my $operation;
    my $line;
    for ($i=$begin;$i<=$end;$i++){
	
	$line=$lines[$i];
	if($line=~/operations/){
	    $flag=1;
	    next;
	}
	if(length($line)<3){
	    next;
	}
	if($flag==1){
	    if($line=~/Match.+\((\d+),\s(\d+)\)/){
		if($1==0) {next;}
		$AliS=$AliS." ".$seqS[$1-1];
		$AliR=$AliR." ".$seqR[$2-1];
	    }
	    if($line=~/dup.+S.+\[(\d+)..(\d+)\]/){
		$lb=$1;$rb=$2;
		if($line=~/Left/){$lb=$lb+1;}
		if($line=~/Right/){$rb=$rb-1;}
		for($k=$rb;$k>=$lb;$k--){
		    $AliS=$AliS." ".$seqS[$k-1];
		    $AliR=$AliR." -";
		}
		#$AliS=$AliS." ".$seqS[$1-1];
		#$AliR=$AliR." ".$seqR[$2-1];
	    }
	    if($line=~/dup.+R.+\[(\d+)..(\d+)\]/){
		$lb=$1;$rb=$2;
		if($line=~/Left/){$lb=$lb+1;}
		if($line=~/Right/){$rb=$rb-1;}
		for($k=$rb;$k>=$lb;$k--){
		    $AliR=$AliR." ".$seqR[$k-1];
		    $AliS=$AliS." -";
		}
		#$AliS=$AliS." ".$seqS[$1-1];
		#$AliR=$AliR." ".$seqR[$2-1];
	    }
	}
	$line_no++;
    }
    
    
    @tmp=split(/ /,$AliS);
    @tmp=reverse @tmp;
    
    $AliS2="\\\$ ";
    for($k=0;$k<scalar @tmp;$k++){
	$AliS2=$AliS2." ".$tmp[$k];	
    }

    @tmp=split(/ /,$AliR);
    @tmp=reverse @tmp;
    
    $AliR2="\\\$ ";
    for($k=0;$k<scalar @tmp;$k++){
	$AliR2=$AliR2." ".$tmp[$k];	
    }
#    $AliS2=reverse $AliS;
#    $AliR2=reverse $AliR;

#    $AliS2="\\\$ ".$AliS2;
#    $AliR2="\\\$ ".$AliR2;
    
    print $AliS2."\n";
    print $AliR2."\n";
    
    my @pos_alignmentS;
    my @pos_alignmentR;

    my $max_no_char=0;
## get pos in alignment    
    @words=split(/\s+/,$AliS2);
    $counter=0;
    $words_size=scalar @words;
    for($i=0;$i<$words_size;$i++){
	if($words[$i] ne "-"){	    
	    push(@pos_alignmentS,$i);
	    $counter++;
	    if((length($words[$i]) >$max_no_char)&&(!($words[$i]=~/\$/)))
	    {
		print" length: ".length($words[$i])."\n";
		$max_no_char=length($words[$i]);}
	}
	#print "Words: ".$words[$i]." length: ".length($words[$i])."\n";
    }
    
    #print "HIIIIIIIIIIIIIIII  ".$max_no_char."\n";
    @words=split(/\s+/,$AliR2);
    $counter=0;
    $words_size=scalar @words;
    for($i=0;$i<$words_size;$i++){
	if($words[$i] ne "-"){	    
	    push(@pos_alignmentR,$i);
	    $counter++;
	}
	if((length($words[$i]) >$max_no_char)&& (!($words[$i]=~/\$/)) )
	{
		$max_no_char=length($words[$i]);
	}

    }

# draw the arcs
# go through a second pass 
    my @curve_listS;
    my @curve_listR;
    $line_no2=0;
    $tmp_line;
    for ($i=$begin;$i<=$end;$i++){
	$line=$lines[$i];
	if($line=~/operations/){
	    next;
	}
	if(length($line)<3){
	    next;
	}

	if($line=~/dup.+S.+\[(\d+)..(\d+)\]/){
	    $j=$i+1;
	    $tmp_line=$lines[$j];
	    while($tmp_line=~/\#/){
		$tmp_line=~/(\d+)\s+\-\>\s+(\d+)/;
		my @tmp_curve=($pos_alignmentS[$1],$pos_alignmentS[$2]);
		push(@curve_listS,\@tmp_curve);
		$j++;
		$tmp_line=$lines[$j];		
	    }
	    $i=$j-1;
	}
	if($line=~/dup.+R.+\[(\d+)..(\d+)\]/){
	    $j=$i+1;
	    $tmp_line=$lines[$j];
	    while($tmp_line=~/\#/){
		$tmp_line=~/(\d+)\s+\-\>\s+(\d+)/;
		my @tmp_curve=($pos_alignmentR[$1],$pos_alignmentR[$2]);
		push(@curve_listR,\@tmp_curve);
		$j++;
		$tmp_line=$lines[$j];		
	    }
	    $i=$j-1;	 
	}


    }

    open(texfile, ">$output_tex.$outputfilecounter.tex");

    my $space_width=12*$max_no_char;
    
    my $char_width=11; # important for the curve
    my $fixed_left_shift=15;
    
    my $pic_size=0;

    
    $char_width=$space_width+2;
    $tmppp=$AliS2;    
    $tmppp=~s/\s+//g;     
    $pic_size=length($tmppp);
    $picture_width=$char_width*($pic_size)+$fixed_left_shift;
 
my  $i=0;
$curve_listS_size=scalar @curve_listS;
$curve_listR_size=scalar @curve_listR;
#    print "size $curve_listS_size \n";
#    print "pos $curve_listS[$i]->[0]\n";
#    print "pos $curve_listS[$i]->[1]\n";
    
#    print $picture_width."\n";
#### printing the tex file

    

    print texfile "\\documentclass[10pt]{article}\n";
    print texfile "\\usepackage{curves}\n";
    
    print texfile "\\addtolength{"."\\oddsidemargin}{-150pt}\n";
    print texfile "\\addtolength{\\evensidemargin}{-150pt}\n";
    print texfile "\\addtolength{\\topmargin}{-120pt}\n";
    print texfile "\\special{papersize=".$picture_width."pt,251pt}\n";
    print texfile "\\textwidth=".$picture_width."pt\n";
    print texfile "\\textheight=251pt\n";
    print texfile "\\pagestyle{empty}\n";
    
    print texfile "\\begin{document}\n";
    
    print texfile "\\setlength{\\unitlength}{1pt}\n";
    
    print texfile "\\begin{picture}($picture_width,245)\n";
    print texfile "\\linethickness{1pt}\n";
    
#    print texfile "\\put(7,140){\\mbox{{\\texttt{ \\begin{huge}X A T R a K K\\end{huge}}}}}\n";
#    print texfile "\\put(7,100){\\mbox{{\\texttt{ \\begin{huge}X A T R a K K\\end{huge}}}}}\n";

#    $AliS2=~s/\s+/\\hspace{10pt}/g;
#    $AliR2=~s/\s+/\\hspace{10pt}/g;
    
    # option1 fixed width characters
    #print texfile "\\put(7,140){\\mbox{{\\texttt{ \\begin{huge}$AliS2\\end{huge}}}}}\n";
    #print texfile "\\put(7,100){\\mbox{{\\texttt{ \\begin{huge}$AliR2\\end{huge}}}}}\n";
    # option 2 table
    $cellpos="";
    #print "hallo: ".$AliS2;
    $lenAli=length($AliS2);
    for($t=0;$t<length($AliS2)-3;$t++){
	$cellpos="|p{$space_width"."pt}".$cellpos;
    }
    $cellpos=$cellpos."p{$space_width"."pt}";
    $AliS2=~s/\s+$//;
    $AliR2=~s/\s+$//; # to remove last space
    #print $AliS2;
    $AliS2=~s/\s+/\&/g;
    $AliR2=~s/\s+/\&/g;
    print texfile "\\put(0,140){\n";
    print texfile "\\begin{huge}\n";
    print texfile "\\setlength{\\tabcolsep}{1pt}\n";
    print texfile "\\begin{tabular}{$cellpos}\n";
    print texfile "$AliS2\\\\ \n";
    print texfile "$AliR2";
    print texfile "\\end{tabular}\n";
    print texfile "\\end{huge}";
    print texfile "}\n";


#    print texfile "\\curve(198,160, 256,185, 308,160)\n";
    for($i=0;$i<$curve_listS_size;$i++){
	$pos1=($char_width)*$curve_listS[$i]->[0]+$char_width/2+$fixed_left_shift;
	$pos2=($char_width)*$curve_listS[$i]->[1]+$char_width/2+$fixed_left_shift;
	if($pos1<$pos2){
	    $pos3=$pos1+abs($pos2-$pos1)/2;
	}
	else{
	    $pos3=$pos2+abs($pos2-$pos1)/2;
	}
	
	$height=(abs($curve_listS[$i]->[1]-$curve_listS[$i]->[0]+1))*7+170;
	print texfile "\% $curve_listS[$i]->[0], $curve_listS[$i]->[1]\n";
	print texfile "\\curve($pos1,170,$pos3,$height, $pos2,170)\n";
    }

    for($i=0;$i<$curve_listR_size;$i++){
	$pos1=($char_width)*$curve_listR[$i]->[0]+$char_width/2+$fixed_left_shift;
	$pos2=($char_width)*$curve_listR[$i]->[1]+$char_width/2+$fixed_left_shift;
	if($pos1<$pos2){
	    $pos3=$pos1+abs($pos2-$pos1)/2;
	}
	else{
	    $pos3=$pos2+abs($pos2-$pos1)/2;
	}

	
	$height=120-((abs($curve_listR[$i]->[1]-$curve_listR[$i]->[0])+1)*7);
	print texfile "\% $curve_listR[$i]->[0], $curve_listR[$i]->[1]\n";
	print texfile "\\curve($pos1,120,$pos3,$height, $pos2,120)\n";
    }

    print texfile "\\end{picture}\n";
    print texfile "\\end{document}\n";
    
    close(texfile);
    
    #get file path
    $output_tex=~/(.*)\/(\w+)/;
    $filepath=$1;
    if(length($filepath)<1){
	$filepath="./";
    }
    
#    print "HALOOOOOOOOOOO: ".$filepath."\n";

    $argstring="latex --output-directory=$filepath $output_tex.$outputfilecounter.tex";

    system($argstring);

    $argstring="dvipng $output_tex.$outputfilecounter.dvi -o $output_tex.$outputfilecounter.png";
    system($argstring);

    $argstring="dvips $output_tex.$outputfilecounter.dvi -o $output_tex.$outputfilecounter.ps";
    system($argstring);

    $argstring="rm -f $output_tex.$outputfilecounter.log";
    system($argstring);
    $argstring="rm -f $output_tex.$outputfilecounter.aux";
    system($argstring);

    #print "HIIIIIIII2 ".$max_no_char."\n";
}
