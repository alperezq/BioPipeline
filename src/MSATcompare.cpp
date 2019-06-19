// MSATcompare.cpp : Defines the entry point for the console application.


#include <stdio.h>
#include <iostream>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <ctype.h>

#include "ThroughputAlignment.h"
//#include "DuplicationHistory.h"


void generate_random(int seq_number,int min_seq_length,int max_seq_length,int alphabet_size, int max_rep);
void print_header();
void print_header2();
void print_help();

using namespace std;

int main(int argc, char* argv[])
{

    //srand((unsigned)time( NULL ));
	srand((unsigned)time(0)); 
	//  timing

	//DuplicationHistory* dobj;
	char ch;
	char temp[500];
	char filename1[100];
	char filename[500];
	//char* filename;
	char filename2[1000];
	char randomfilename[1000];

	double score1,score2;
	double bensonleft,oursleft,bensonright,oursright;
	double alignscore1;
	double alignscore2;
	int only_left_dup_flag=0; // with right duplication
	int alphadep_flag=0; // 1= use alohabet dependent algorithm
	int rle_flag=0;
	int in_show_statistics_flag=0;
	int in_show_cost_file_flag=0;
/////
	char* codesfilename;
	char* costs_filename;
	
	
	
	/// flag to consider reverse or not
	int rev_seq_flag =0; //0=Do not consider reverse

	long i;
	int align_flag=0;
	int show_alignment_flag=0;
	int show_breakeven_flag=0;
	int flag=0;
	int allow_insertion_flag=0;
	if(argc<2)
	{
	    print_header();
	    print_help();
	    return 0;
	}
	

	/////////////// Parsing options /////////////////
	for(int j=1;j<argc;j++){
		if(strcmp(argv[j],"-f")==0){
			codesfilename=argv[j+1];	
			//cout << codesfilename << endl;
			j++;
			continue;
		}
		else if(strcmp(argv[j],"-onlyleft")==0){
			only_left_dup_flag=1;
			continue;
		}
		else if(strcmp(argv[j],"-onlyright")==0){
			only_left_dup_flag=1;
			rev_seq_flag =1;	
			continue;
		}
		else if(strcmp(argv[j],"-cfile")==0){
			costs_filename=argv[j+1];
			//cout << costs_filename << endl;
			j++;
			continue;
		}
		else if(strcmp(argv[j],"-alphadep")==0){
		    alphadep_flag=1;
		    continue;
		}
		else if(strcmp(argv[j],"-align")==0){
		    align_flag=1;
		}
		else if(strcmp(argv[j],"-showalign")==0){
		    show_alignment_flag=1;
		    align_flag=1;
		}
		else if(strcmp(argv[j],"-breakeven")==0){
		    show_breakeven_flag=1;
		    align_flag=1;
		}
		else if(strcmp(argv[j],"-insert")==0){
		    allow_insertion_flag=1;
		}
		else if(strcmp(argv[j],"-rle")==0){
		    rle_flag=1;
		}
		else if(strcmp(argv[j],"-showstat")==0){
		    in_show_statistics_flag=1;
		}
		else if(strcmp(argv[j],"-showcost")==0){
		    in_show_cost_file_flag=1;
		}
		else if(strcmp(argv[j],"-random")==0){
		    if(argc<j+6){
		      print_header();
		      cerr << "Error: Missing options for ranodm option\n";
		      exit(-1);
		    }
		    generate_random(atoi(argv[j+1]),atoi(argv[j+2]),atoi(argv[j+3]),atoi(argv[j+4]),atoi(argv[j+5]));
		    exit(0);
		}
		else if(strcmp(argv[j],"-evaluate")==0){
		    
		    ThroughputAlignment TheEvalObj;
		    if(argc<j+3){
		      print_header();
		      cerr << "Error: Missing options for  evaluation option\n";
		      exit(-1);
		    }
		    print_header2();
		    TheEvalObj.Evaluate(argv[j+1],argv[j+2]);
		    exit(0);
		}
		else{
		  print_header();
		  cerr << "Unkown option: " << argv[j] << " \n";
		  exit(-1);
		}
	}
	
	
	
	//start processing
	print_header2();



        // DuplicationHistory* Alignemntobj;

	ThroughputAlignment*  ThAligObj= new ThroughputAlignment;
	ThAligObj->set_align_flag(align_flag);
	ThAligObj->set_alphadep_flag(alphadep_flag);
	ThAligObj->set_show_alignment_flag(show_alignment_flag);
	ThAligObj->set_show_breakeven_flag(show_breakeven_flag);
	ThAligObj->set_allow_insertion_flag(allow_insertion_flag);
	ThAligObj->set_rle_flag(rle_flag);
	ThAligObj->set_show_flag(in_show_statistics_flag,in_show_cost_file_flag);
	//ThAligObj->RandomDataSet("Randomdataset.txt",37);
	//ThAligObj->RandomDataSet("Randomdataset.txt",37);


	
        // getchar();
	
	//strcpy(codesfilename,"Randomdataset.txt");
	//strcpy(codesfilename,"ReducedDataSet.txt");
	//strcpy(codesfilename,"MSY1codes.txt");
	string tempstr;
	if((only_left_dup_flag==1)&&(rev_seq_flag==0)){
	    tempstr.append(codesfilename);
	    tempstr.append(".align.onlyleft");
	}
	else if((only_left_dup_flag==1)&&(rev_seq_flag==1)){
	    tempstr.append(codesfilename);
	    tempstr.append(".align.onlyright");
	}
	else{
	    tempstr.append(codesfilename);
	    tempstr.append(".align");
	}
	
	//strcpy(filename,"AlignmentResult.txt");
	
	int retval=ThAligObj->AlldataSetAlign((char*)tempstr.c_str(),only_left_dup_flag,codesfilename,rev_seq_flag,costs_filename);

	delete ThAligObj;
	
	return retval;

}


void generate_random(int seq_number,int min_seq_length,int max_seq_length,int alphabet_size, int max_rep){

    string map_seq;
    int seq_range=max_seq_length-min_seq_length;
    float val=0;
    char* type=new char[alphabet_size];
    for(int j=0;j<seq_number;j++){	
	
	val=((float)rand())/(float) RAND_MAX;		
	
	int length=(int)(val*seq_range)+min_seq_length;
	//cout << "random "<< val <<" "<< length << " \n";
	for(int j=0;j<length+1;j++){
	    // choose an alphabet chracter
	    val=((float)rand())/RAND_MAX;		
	    int alphabet_char=(int)(val*alphabet_size);
	    // decide repetiton
	    val=((float)rand())/RAND_MAX;		
	    int rep_val=(int)(val*max_rep)+1;	    
	    int rep_val2=min(rep_val,length+1-j);
	    j=j+min(rep_val+1,rep_val2);
	    for(int k=0;k<rep_val2;k++){
		sprintf(type,"%d ",alphabet_char);
		//cout << "halloo "<< alphabet_char << " "<<type;
		map_seq.append((char*)type);
	    }
	    
	}
	//for(int i=0;i<map_seq.size();i++){
	cout << map_seq;
	    //}
	cout << endl;
	map_seq.clear();
    }

  
    //string random_score=
    ofstream out("random_score",ios::out);
    out << "# Type no. " << alphabet_size << endl;
    out << "# Types ";
    for(int i=0;i<alphabet_size;i++){
	out << i << " ";
    }
    out << endl;
    out << "# Indel align 20" << endl;
    out << "# Indel hist 20"<< endl;
    out << "# Dup 1"<< endl;
    out << "# matrix\n";
    for(int i=0;i<alphabet_size-1;i++){
	for(int j=i;j<alphabet_size-1;j++){
	    val=((float)rand())/RAND_MAX;		
	    int incost=int(val*10);
	    out << incost << " ";
	}
	out << endl;
    }    
    out.close();
}

void print_header(){
    cout << "*/////////////////////////////////////////////////////////////*\n";
    cout << "*                 This is ARLEM  version 1.0                  *\n";
    cout << "*          Copyright by Mohamed I. Abouelhoda  (C) 2007       *\n";
    cout << "*  Plz. cite Abouelhoda, Giegerich, Behzadi, and Steyaert,    *\n";
    cout << "*                    APBC 2008, and JBCB                      *\n";
    cout << "*/////////////////////////////////////////////////////////////*\n";
}

void print_header2(){
    cout << "*/////////////////////////////////////////////////////////////*\n";
    cout << "***************************************************************\n";
    cout << "********                                               ********\n";
    cout << "********                      ARLEM                    ********\n";
    cout << "********                                               ********\n";
    cout << "***************************************************************\n";
    cout << "*/////////////////////////////////////////////////////////////*\n";
    cout << "*                   This is ARLEM version 1.0                 *\n";
    cout << "*                                                             *\n";
    cout << "*        Copyright by Mohamed I. Abouelhoda  (C) 2007         *\n";
    cout << "*                                                             *\n";
    cout << "*       Unauthorized commercial usage and distribution        *\n";
    cout << "*              of this program is prohibited.                 *\n";
    cout << "*             Contact the author for a license.               *\n";
    cout << "*           Please report bugs and suggestions to             *\n";
    cout << "*                <mohamed.ibrahim@uni-ulm.de>                 *\n";  
    cout << "*   Plz. cite Abouelhoda, Giegerich, Behzadi, and Steyaert    *\n";
    cout << "*/////////////////////////////////////////////////////////////*\n";
    cout << "***************************************************************\n";
}
void print_help(){
    cout << "Usage: \n";
    cout << "> arlem -f map_file -cfile cost_file [options] \n";
    cout << "\n";
    cout << "SYNOPSIS: \n";
    cout << "       The program, without options, computes duplication history score\n";
    cout << "       originated from left and right unit for all input sequences. \n";
    cout << "OPTIONS:\n";
    cout << "-align            : compute pairwise alignment for all genomes \n";
    cout << "-onlyleft         : do not allow right duplications \n";
    cout << "-onlyright        : do not allow left duplications \n";
    cout << "-insert           : allow free insertions in the history and alignment \n";
    cout << "-alphadep         : use alphabet dependent algorithm \n";
    cout << "-rle              : use rle algorithm incorporating insetions,\n";
    cout << "                    useful if map conatains repetitions\n";
    cout << "-showalign        : display alignment, this option sets \"-align\" option \n";
    cout << "-breakeven        : display where in alignment half optimal score occurs \n";
    cout << "-showstat         : show statistics about input sequences\n";
    cout << "-showcost         : show cost file\n";
    cout << "-evaluate a b     : compare alignment results for 2 datasets, where \n";
    cout << "                    a is file containing alignments with left\n";
    cout << "                      and right duplications \n";
    cout << "                    b is file storing alignments without left/right duplications\n";
    cout << "-random a b c d e : generate random data, where \n";
    cout << "                    a is the number of map sequences \n";
    cout << "                    b is the min. sequence length\n";
    cout << "                    c is the max. sequence length\n";
    cout << "                    d is the alphabet_size\n";
    cout << "                    e is the max_number of unit repetitions\n";
}
