// ThroughputAlignment.cpp: implementation of the ThroughputAlignment class.
//
//////////////////////////////////////////////////////////////////////

//#include "stdafx.h"
#include "ThroughputAlignment.h"


#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string.h>


#include <stdlib.h>
#include <string>
#include <vector>
//#include <windows.h>
#include <time.h>

#ifndef max
#define max(a, b)  ((a)>=(b) ? (a) : (b))
#endif
#ifndef min
#define min(a, b)  ((a)<=(b) ? (a) : (b))
#endif

using namespace std;

#define SPECIALSYMBOLE "$"

#define BigValue2 99999

void fsplit (string line, const char* sep, vector<string> & words);
//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////


ThroughputAlignment::ThroughputAlignment()
{
	int i=0;
	//for(i=0;i<MaxSeqSize;i++){
	    //sequence[i]=0;
	    //outsequence[i]=0;
	//}
	//for(i=0;i<MaxAlphaSize;i++){
	    //strcpy(GivenAlphabet[i],"");
	//}
/*
	Mapping[0]='1';
	Mapping[1]='2';
	Mapping[2]='3';
	Mapping[3]='4';
	Mapping[4]='5';
	Mapping[5]='6';
	Mapping[6]='7';
	Mapping[7]='8';
	Mapping[8]='9';
	Mapping[9]='10';
	Mapping[10]='11';
	Mapping[11]='12';
*/
//	cout << "Constructed "<< endl;

	alphadep_flag=0; // flag to use alphabet dependent algorithm
	show_alignment_flag=0;
	show_breakeven_flag=0;
	allow_insertion_flag=0;
	rle_flag=0;

	show_statistics_flag=0;
	show_cost_file_flag=0;
}

ThroughputAlignment::~ThroughputAlignment()
{
    types_vector.clear();
    if(in_dist!=NULL){
	
	for(int i=0;i<type_no;i++){
	    delete [] in_dist[i];	    
	}
	delete [] in_dist;
    }
    in_dist=NULL;
    
}

///////////////////////////////////////////////////////////
// Computes all pairwise alignment for the maps in
// the input codesfilename.
// Input: filename is output file
//		: flag: to decide if to ignor left duplication or not
//		: in_rev_flag to use reverse of the sequnce or not
///////////////////////////////////////////////////////////

int ThroughputAlignment::AlldataSetAlign(char* filename,int left_dup_flag, char* codesfilename, int in_rev_flag,char* costs_file)
{

	int ourflag=in_rev_flag;
	char ch;
	char temp[5000];
	int score1,score2;
	int bensonleft,oursleft,bensonright,oursright;
	int alignscore1;
	int alignscore2;
	int i,j;
	
	TestingDupHist* Alignemntobj;
	TestingDupHist* HistObjArray;
	TestingDupHist* AligObj;
//	FILE* fptr;
//	FILE* fptr2;

	//cout << "Parse: "<< codesfilename <<" \n";


//	fptr2=fopen("AlignmentTiming.txt","w");
//	fptr=fopen(filename,"w");	
	long size=sizeof(TestingDupHist);
	//cout << "Memory: "<< size <<endl;
	
	read_cost_file(costs_file);

	
	get_file_statistics(codesfilename);		
		
	if(check_traingularity_of_cost()!=1){
	    //delete ThAligObj;
//	    fclose(fptr);
	    //fclose(fptr2);	    
	    return(-1);
	}
	

	//HistObjArray=(TestingDupHist*)malloc(sizeof(TestingDupHist)*number_of_sequences);
	HistObjArray=new TestingDupHist [number_of_sequences];		

	clock_t start, finish;
	double  duration;
	start = clock();


	

	// test with random data
	//if(flag==0)
	//this->RandomDataSet("Randomdataset.txt",30);
	
	cout << "#------------ Start Processing --------- \n" ;
	
	this->parse_expanded_datasetfile(codesfilename,NULL,HistObjArray,in_rev_flag);
	
	

	finish = clock();
	duration = (double)(finish - start) / CLOCKS_PER_SEC;
	cout << "Time for duplication history in sec.: " <<((double)(finish - start) / CLOCKS_PER_SEC)<< endl;
	if(align_flag==0){
	    
	    
	    for(long i=0;i<number_of_sequences;i++){
		//HistObjArray[i].FreeMEMORY();
		delete[] HistObjArray[i].sequence1;	    
	    }
	    	    
	    delete [] HistObjArray;
	    return (1);
	}

	cout << "#------------ Alignment Phase ----------- \n" ;
	

	finish = clock();
	duration = (double)(finish - start) / CLOCKS_PER_SEC;	
//	fprintf(fptr2,"Time for constructing duplication history: %2f sec == %2f min \n",((double)(finish - start) / CLOCKS_PER_SEC),
//		((double)(finish - start) / ((double) 60*CLOCKS_PER_SEC)));

	

	start = clock();

//	validsequencenumber=50;
	
	for(i=0;i<validsequencenumber-1;i++){
		for(j=i+1;j<validsequencenumber;j++){
		    if(show_alignment_flag)
		      cout << "Aligning Seq1: "<<i << "  Seq2: " <<j << endl;
		    AligObj= new TestingDupHist; 
		    AligObj->Intialize(in_dist);
		    if(show_alignment_flag){
			AligObj->set_show_alignment_flag(1,types_vector);
		    }
		    if(show_breakeven_flag){
			AligObj->set_breakeven_flag(1,types_vector);
		    }
		    if(left_dup_flag==0){			    
			if(rle_flag==0){
			    alignscore1=AligObj->ModAlignment(&HistObjArray[i],&HistObjArray[j]);			
			}
			else{
			    alignscore1=AligObj->AlignmentRLE(&HistObjArray[i],&HistObjArray[j],0);	
			}
		    }
		    else{
			if(rle_flag==0){
			    alignscore1=AligObj->AlignmentWithoutRightDup(&HistObjArray[i],&HistObjArray[j]);	
			}
			else{
			    alignscore1=AligObj->AlignmentRLE(&HistObjArray[i],&HistObjArray[j],1);
			}
		    }
		    
		    cout << "Score of aligning Seq:"<<i << ", Seq:" <<j << " ="<<alignscore1 << "\n";
		    cout << "#---------------------------------------- \n" ;		    
//		    fprintf(fptr,"%d %d Score %d\n", i,j,alignscore1);
		    		    
		    delete AligObj;
		}
	}



	finish = clock();
	duration = (double)(finish - start) / CLOCKS_PER_SEC;	
//	fprintf(fptr2,"Time for computing all alignment: %2f sec == %2f min \n",((double)(finish - start) / CLOCKS_PER_SEC),
//		((double)(finish - start) / ((double) 60*CLOCKS_PER_SEC)));
	cout << "Time for alignment in sec.: " <<((double)(finish - start) / CLOCKS_PER_SEC)<< endl;


	//cout << "cleaning " << number_of_sequences << endl;
	for(long i=0;i<number_of_sequences;i++){
	    //HistObjArray[i].FreeMEMORY();
	    delete[] HistObjArray[i].sequence1;	    
	}
	//free(HistObjArray);
	delete[] HistObjArray;
//	fclose(fptr);
//	fclose(fptr2);

	return(0);

}


void ThroughputAlignment::parsedatasetfile(char* filename,char* desfilename,TestingDupHist* HistArrayObj, int revflag)
{
}

char ThroughputAlignment::MakeMapping(char* token){

	if(strcmp(token,"0")==0){
		return('1');
	}
	else if(strcmp(token,"1")==0){
		return('2');
	}
	else if(strcmp(token,"1a")==0){
		return('3');
	}
	else if(strcmp(token,"2")==0){
		return('4');
	}
	else if(strcmp(token,"3")==0){
		return('5');
	}
	else if(strcmp(token,"3a")==0){
		return('6');
	}
	else if(strcmp(token,"4")==0){
		return('7');
	}
	else if(strcmp(token,"4a")==0){
		return('8');
	}
	else if(strcmp(token,"4*")==0){
		return('9');
	}
	else{
		cout << "CHARCHTER is not in the alphabet\n";
		return(0);
	}

}
void ThroughputAlignment::Evaluate(char* file1, char* file2)
{

	FILE* fptr1;
	FILE* fptr2;
//	FILE* outfptr;
	int i,j,score1,score2;
	int count=0;
	int aligncount=0;
	int twotimes=0;
	int threetimes=0;
	int fourtimes=0;
	int fivetimes=0;
	int sixtimes=0;
	int highertimes=0;
	int equal=0;
	int less=0;
	float ratio;

	//class myrecord{int seq1; int seq2; int score;};
	std::vector<string> all_records_file1;  
	std::vector<string> all_records_file2; 
	
	string line;
	//cout << "HALOOOOOOOOOOOOOOOO " << file1 << endl;
	ifstream myfile (file1);
	if (myfile.is_open())
	{
	    while (! myfile.eof() )
	    {
		getline (myfile,line);
		if(line.find("Score of aligning",0) != string::npos){
		    //cout << "Haloooooooo: " << line << endl;
		    all_records_file1.push_back(line);
		}
		
	    }
	    myfile.close();
	}
	else{
	    cerr << "Error: Unable to open file: "<< file1 << endl;;
	}

	ifstream myfile2 (file2);
	if (myfile2.is_open())
	{
	    while (! myfile2.eof() )
	    {
		getline (myfile2,line);
		if(line.find("Score of aligning",0) != string::npos){
		    //cout << "Haloooooooo: " << line << endl;
		    all_records_file2.push_back(line);
		}
		
	    }
	    myfile2.close();
	}
	else{
	    cerr << "Error: Unable to open file: "<< file2 << endl;;
	}


/*

	fptr1=fopen(file1,"r"); // file 1 is the complete model
	if(fptr1==NULL){
	    cerr << "Error: non existing file "<< file1<< endl;
	    return;
	}
	fptr2=fopen(file2,"r"); // file 2 is without
	if(fptr2==NULL){
	    cerr << "Error: non existing file "<< file2<< endl;
	    return;
	}
	
	outfptr=fopen(outfilename,"w");
	if(outfptr==NULL){
	    cerr << "Error: cannot open file "<< outfptr << endl;
	    return;
	}
	
	// reading first file
	char line[2000];
	while(feof(fptr1)==0){
	    
	}
	
	
	outfptr=fopen(outfilename,"w");
	if(outfptr==NULL){
	    cerr << "Error: cannot open file "<< outfptr << endl;
	    return;
	}
*/
//	while(feof(fptr1)==0){
	int kk=0;
	int size1=all_records_file1.size();
	int size2=all_records_file2.size();
	if(size1 !=size2){
	    cerr << "Error: it seems to compare unrelated files\n";
	    all_records_file1.clear();
	    all_records_file2.clear();
	    exit(1);
	}

	while(kk<size1){
	    // cout << all_records_file1[kk].c_str() << endl;
	    //kk++;
	    //continue;
	    //fscanf(fptr1,"%d %d Score %d\n", &i,&j,&score1);
	    //	fscanf(fptr2,"%d %d Score %d\n", &i,&j,&score2);
	    int bytes1=sscanf(all_records_file1[kk].c_str(),"Score of aligning Seq:%d, Seq:%d =%d", &i,&j,&score1);
	    int bytes2=sscanf(all_records_file2[kk].c_str(),"Score of aligning Seq:%d, Seq:%d =%d", &i,&j,&score2);

//	    if(bytes1<20) continue;
		if(score1<score2){
			count++;
			score1++;
			score2++;
			ratio=(float)score2/(float)score1;
			//fprintf(outfptr,"Sequence Pair: %d %d, Score1 %d, Score2 %d\n",i,j,score1-1,score2-1);
			cout << "Sequence Pair: "<< i << " " <<j << " with scores" << score1-1 << " and " << score2-1 << endl;
			if((1<ratio)&&(ratio<=2.5)){
				twotimes++;
			}
			else if((2.5<ratio)&&(ratio<=3.5)){
				threetimes++;
			}
			else if((3.5<ratio)&&(ratio<=4.5)){
				fourtimes++;
			}
			else if((4.5<ratio)&&(ratio<=5.5)){
				fivetimes++;
			}
			else if((5.5<ratio)&&(ratio<=6.5)){
				sixtimes++;
			}
			else{
				highertimes++;
			}
		}
		else if(score1==score2){
			equal++;
		}
		else{
			less++;

		}
		aligncount++;
		kk++;
		//fprintf(fptr1,"%d %d Score %d\n", &i,&j,&score);
	}

	cout <<"Number of total Alignment: "<< aligncount<<"\n";
	cout <<"No. of equal score alignment: "<< equal<<"\n";
	cout <<"No. of less score alignment. "<< less<<"\n";
	cout <<"Number of improved alignment: "<<count << "\n";
	cout <<"No. of two times Improvement: "<< twotimes<<"\n";
	cout <<"No. of three times Improvement: "<< threetimes<<"\n";
	cout <<"No. of four times Improvement: "<< fourtimes<<"\n";
	cout <<"No. of five times Improvement: "<< fivetimes<<"\n";
	cout <<"No. of six times Improvement: "<< sixtimes<<"\n";
	cout <<"No. of heigher times Improvement: "<< highertimes<<"\n";
	cout <<"No. of heigher times Improvement: "<< highertimes<<"\n";
/*	
	fprintf(outfptr," Number of total Alignment: %d\n",aligncount );
	fprintf(outfptr," Number of equal score alignment: %d\n",equal );
	fprintf(outfptr," Number of less score alignment: %d\n",less );
	fprintf(outfptr," Number of improved alignment: %d\n",count );
	fprintf(outfptr," Number of two times improvement: %d\n",twotimes );
	fprintf(outfptr," Number of three times improvement: %d\n",threetimes );
	fprintf(outfptr," Number of four times improvement: %d\n",fourtimes );
	fprintf(outfptr," Number of five times improvement: %d\n",fivetimes );
	fprintf(outfptr," Number of six times improvement: %d\n",sixtimes );
	fprintf(outfptr," Number of higher times improvement: %d\n",highertimes );

//	fclose(fptr1);
//	fclose(fptr2);

	fclose(outfptr);
*/
	all_records_file1.clear();
	all_records_file2.clear();

}


/////////////////////////////////////////////////
// For testing

void ThroughputAlignment::RandomDataSet(char* filename, int seq_number){
	
	FILE* fptr;
	float val;
	int length;
//	SYSTEM_INFO objinfo;
//	GetSystemInfo(&objinfo);
	float dupsize=15;
	float mindupsize=0;
	
	fptr=fopen(filename,"w");

	//	cout << "The random length:"<<length <<"\n" ;
	for(int j=0;j<seq_number;j++){

		val=((float)rand())/RAND_MAX;		
		length=(int)val*68+50;
		length=min(length,MaxSeqSize);
		int offset=1;
		for(int i=0;i<length;i=i+(int)val){			
			// get unit type
			int randmax=RAND_MAX;
			val=((float) rand()/(float)RAND_MAX);
			val=val*(float)75;
			// scale is 15	
			if(val<=(float)15){ //  unit 0
				//get length
				val=((float)rand())/RAND_MAX;
				val=val*dupsize+mindupsize;
				val=max(val,1);
				//val=min(length-i,val);
				if(((int)val+i)<length)
					fprintf(fptr,"(%s)%d","0",(int)val);
			}
			else if(val<=(float)22.5){ //unit 1
				//get length
				val=((float)rand())/RAND_MAX;
				val=val*dupsize+mindupsize;
				val=max(val,1);
				//val=min(length-i,val);
				if(((int)val+i)<length)
					fprintf(fptr,"(%s)%d","1",(int)val);
			}
			else if(val<=(float)30){ // unit 1a
				//get length
				//get length
				val=((float)rand())/RAND_MAX;
				val=val*dupsize+mindupsize;
				val=max(val,1);
				//val=min(length-i,val);
				if(((int)val+i)<length)
					fprintf(fptr,"(%s)%d","1a",(int)val);
			}
			else if(val<=(float)45){ // unit 2
				//get length
				val=((float)rand())/RAND_MAX;
				val=val*dupsize+mindupsize;
				val=max(val,1);
				//val=min(length-i,val);
				if(((int)val+i)<length)
					fprintf(fptr,"(%s)%d","2",(int)val);
			}
			else if(val<=52.5){ // unit 3
				//get length
				val=((float)rand())/RAND_MAX;
				val=val*dupsize+mindupsize;
				val=max(val,1);
				//val=min(length-i,val);
				if(((int)val+i)<length)
					fprintf(fptr,"(%s)%d","3",(int)val);
			}
			else if(val<=60){ // unit 3a
				//get length
				val=((float)rand())/RAND_MAX;
				val=val*dupsize+mindupsize;
				val=max(val,1);
				//val=min(length-i,val);
				if(((int)val+i)<length)
					fprintf(fptr,"(%s)%d","3a",(int)val);
			}
			else if(val<=65){ // unit 4
				//get length
				val=((float)rand())/RAND_MAX;
				val=val*dupsize+mindupsize;
				val=max(val,1);
				//val=min(length-i,val);
				if(((int)val+i)<length)
					fprintf(fptr,"(%s)%d","4",(int)val);
			}
			else if(val<=70){ // unit 4a
				//get length
				val=((float)rand())/RAND_MAX;
				val=val*dupsize+mindupsize;
				val=max(val,1);
				//val=min(length-i,val);
				if(((int)val+i)<length)
					fprintf(fptr,"(%s)%d","4a",(int)val);
			}
			else { // unit 4*
				//get length
				val=((float)rand())/RAND_MAX;
				val=val*dupsize+mindupsize;
				val=max(val,1);
				//val=min(length-i,val);
				if(((int)val+i)<length)
					fprintf(fptr,"(%s)%d","4*",(int)val);
			}			
				
		}
		fprintf(fptr,"\n");
	}
	fprintf(fptr,"\n");
	fclose(fptr);
}

////////////////////////////////////////////////////////////////////////////////////
void ThroughputAlignment::parse_expanded_datasetfile(char* filename,char* desfilename,TestingDupHist* HistArrayObj, int revflag)
{
	
	FILE *outfptr;
	int i=0;
	int count=0;
	int seq_id=0;
	int seq_len=0;
	int j=0;
	int unitnumber;
	int alphabetflag=0;
	int unitnumbermatrix[MaxAlphaSize];
	int alphbettypematrxi[MaxAlphaSize];
	int elementnumber=0;
	int minsequencelen=10000;
	int maxsequencelen=0;
	int averagesequencelen=0;
	
	

//	fptr=fopen(filename,"r");	
	ifstream fptr(filename,ios::in);
	if(fptr==0){
	    cerr <<"Error: cannot open input file " << filename << endl;
	    exit(-1);
	}

	//	outfptr=fopen("Statistics.txt","w");
	
//	if(fptr==NULL){
//		cout << "Non existing file"<< "\n";
//		return;
//	}
//	char token[1];
	char* codes_token;
	count=0; // count of the number of input sequences
	inputalphabet=0;

//	while(feof(fptr)==0){	
	string myline="";
	vector <string> words;
	int header_flag=0;
	while(getline(fptr,myline)){
	    if(myline.compare("\n")==0)continue;
	    if(myline.compare(" ")==0)continue;
	    if(myline.compare("")==0)continue;
	    // check for a header line
	    for(int tt=0;tt<myline.size();tt++){
		if(myline.c_str()[tt]=='>'){
		    //cout << "HEADDEEEEEEEEEEEEEEEEEEEER\n";
		    header_flag=1;
		    break;
		}
	    }
	    if(header_flag==1){
		header_flag=0;
		continue;
	    }
	    
	    words.clear();
	    words.push_back(SPECIALSYMBOLE);
	    
	    fsplit(myline, " ", words);	    	    
	      
	    seq_id++;	       
	    HistArrayObj[count].sequence1=new int [words.size()+2];	    

	    for(long j=0;j<words.size();j++){		
		// check if its existing in the alphabetlist and gets its id
		int found_flag=0;
		int k=0;
		for(k=0;k<types_vector.size();k++){		    
		    if(types_vector[k].compare(words[j])==0){
			found_flag=1;
			break;			
		    }			    
		}		
		//cout << "KKKKKKKKKKKKKKKKKK  " << k<< " size " << types_vector.size() << endl;
		if(found_flag==0){
		    cerr << "Error: Unkown character: " << words[j] << " nearly at line "<< count << endl;
		    exit(-1);
		}
		HistArrayObj[count].sequence1[j]=k; // i is the alphabet index, which is dealt with further
		//cout << HistArrayObj[count].sequence1[j];
	    }	    
	    
	    
	    seq_len=words.size();
	    HistArrayObj[count].reverseflag=revflag;
	    HistArrayObj[count].seqlength=seq_len;
	    HistArrayObj[count].sequence1[seq_len]=0;
	    HistArrayObj[count].fixed_dup_cost=duplication_cost;
	    HistArrayObj[count].indel_hist=indel_hist;
	    HistArrayObj[count].indel_align=indel_align;
	    HistArrayObj[count].inputalphabet=types_vector.size();
	    HistArrayObj[count].Alphasize=types_vector.size();
	    HistArrayObj[count].set_alphadep_flag(alphadep_flag);

	    HistArrayObj[count].allocate_array(seq_len,alphadep_flag);
	    

	    //fprintf(outfptr,"%d %d %s\n",count,seq_id,HistArrayObj[count].sequence1);

	    

	    HistArrayObj[count].Intialize(in_dist);

	    //cout <<  "HALOOOOOOOOOOOOOOOOOOOO \n";
	    HistArrayObj[count].getsequence(HistArrayObj[count].sequence1,seq_len);			

	    // process sequence
	    int tmp_val=0;
	    
	    if(alphadep_flag==0){
		if(rle_flag==1){
		    //cout << "Using RLE algorithm\n";
		    tmp_val=HistArrayObj[count].SymmBensonDongWithInsertionRLE();
		}
		else{
		    if(allow_insertion_flag){		    		    
			tmp_val=HistArrayObj[count].SymmBensonDongWithInsertion();			
		    }
		    else{		   
			tmp_val=HistArrayObj[count].SymmBensonDong();			    
		    }
		}
	    }
	    else{
		cout << "using Alphabet-dependent algorithm\n";		
		tmp_val=HistArrayObj[count].constructhistory(NULL);
	    }
	    count++;			
	    cout << "Processed Seq.: " << count << " Score: "<< tmp_val<<", including costs for special unit \"$\"\n";

	}

	validsequencenumber=count;

	fptr.close();
	//fclose(outfptr);
	words.clear();
}

void ThroughputAlignment::get_file_statistics(char* filename){

	
	int i=0;
	int count=0;
//	int seq_id=0;
	int j=0;	


	int minsequencelen=10000;
	int maxsequencelen=0;
	int averagesequencelen=0;
	int seq_len;

	
	if(show_statistics_flag){
	    cout << "*********************************************\n";
	    cout << "Getting file statisics of input map file "<<filename << endl;
	}

	
	char* codes_token;
	count=0;
	
	ifstream fptr(filename,ios::in);
	if(fptr==0){
	    cerr <<"Error: cannot open input file " << filename << endl;
	    exit(-1);
	}

	string myline="";
	
	vector <string> words;
	myline.clear();
	int header_flag=0;
	while(getline(fptr,myline)){
	    
	    if(myline.compare("\n")==0)continue;
	    if(myline.compare(" ")==0)continue;
	    if(myline.compare("")==0)continue;


	    if(myline.c_str()[0]=='>'){
		//cout << "HEADDEEEEEEEEEEEEEEEEEEEER\n";
		continue;
	    }
	    // check for a header line
	    for(int tt=0;tt<myline.size();tt++){
		if(myline.c_str()[tt]=='>'){
		    //cout << "HEADDEEEEEEEEEEEEEEEEEEEER\n";
		    header_flag=1;
		    break;
		}
	    }
	    if(header_flag==1){
		header_flag=0;
		continue;
	    }

	    //cout << "Seq: "<< myline << endl;
	    words.clear();
	    fsplit(myline," ", words);	    	    
	    //cout << "Output "<< words.size() <<endl;
	    for(long j=0;j<words.size();j++){		
		//cout << "haloo "<< words[j] << endl;
		// check if its existing in the alphabetlist and gets its id
		int found_flag=0;
		for(i=0;i<types_vector.size();i++){		    
		    if(types_vector[i].compare(words[j])==0){
			found_flag=1;
			break;			
		    }			    
		}		
		if(found_flag==0){
		    cerr << "Error: Unkown character: " << words[j] << " nearly at line "<<  count << endl;
		    exit(-1);
		}		
	    }	    
	    seq_len=words.size();
	    count++;
	    minsequencelen=min(minsequencelen,seq_len);
	    maxsequencelen=max(maxsequencelen,seq_len);
	    averagesequencelen=averagesequencelen+seq_len;		

	}						
	char sptr[1000];	
	min_seq_len=minsequencelen;
	max_seq_len=maxsequencelen;
	average_seq_len=int((float)averagesequencelen/(float)(count));
	number_of_sequences=count;		

	
	if(show_statistics_flag){
	    cout << "Number of input sequences:" << count << endl;
	    cout << "Min Sequence Length: " << minsequencelen << endl;	
	    cout << "Max Sequence Length: " << maxsequencelen << endl;
	    cout << "Average Sequence Length: " << (float)averagesequencelen/(float)(count) << endl;
	}
	fptr.close();
	words.clear();
		
}


void ThroughputAlignment::read_cost_file(char* costfilename){
	//FILE* fptr;

    type_no=0;
    int error_flag=0;
	ifstream in_file(costfilename,ios::in);
	string::size_type a = 0, e;
	
	vector <string> words;
	int matrix_flag=0;

	if(in_file==0){
		cerr << "Error: unable to open costs file \n";
		exit(-1);
	}
	string zeile;
	while(getline(in_file,zeile)){
	    if(zeile.compare("\n")==0)continue;
	    if(zeile.compare("")==0)break;
	    if(zeile.compare(" ")==0)continue;
	    if(zeile.length()==0)break;
	    words.clear();
	    fsplit(zeile, " ", words);
	    for(int i=0;i<words.size();i++){
		if(words[i].compare("no.")==0){
		    type_no=atoi(words[words.size()-1].c_str());
		    type_no++;
		    break;
		}
		if(words[i].compare("Types")==0){
		    for(int j=i+1;j<words.size();j++){
			types_vector.push_back(words[j]);
		    }
		    types_vector.push_back(SPECIALSYMBOLE);
		    break;
		}
		//if(words[i].compare("align")==0){
		//    indel_align=atoi(words[words.size()-1].c_str());
		//    break;
		//}
		//if(words[i].compare("hist")==0){
		//    indel_hist=atoi(words[words.size()-1].c_str());
		//    break;
		//}
		if(words[i].compare("Indel")==0){
		    indel_hist=atoi(words[words.size()-1].c_str());
		    indel_align=indel_hist;
		    break;
		}
		if(words[i].compare("Dup")==0){
		    duplication_cost=atoi(words[words.size()-1].c_str());
		    break;
		}
		if(words[i].compare("matrix")==0){
		    matrix_flag=1;
		    break;
		}
	    }
	    if(matrix_flag)break;
	}
	
	//// Create distance matrix
	in_dist=new int* [type_no];
	for(int i=0;i<type_no;i++){
	    in_dist[i]=new int [type_no];
	    for(int j=0;j<type_no;j++){
		in_dist[i][j]=0;
	    }
	}
	

		
	if(show_cost_file_flag){
	    //cout << "*********************************************\n";
	    cout << "# Cost file: " << costfilename << endl;
	    cout << "# Type no.: " << type_no << " (including one extra special symbol)"<<endl;
	    cout << "# Cost file: \n";
	    cout << "# Types: ";
	    for(int i=0;i<types_vector.size();i++){
		cout << types_vector[i] << " ";
	    }
	    cout << endl;
	    //	    cout << "# Indel cost in Alignment: " << indel_align << endl;
	    // cout << "# Indel cost in dupl. history:  " << indel_hist << endl;
	    cout << "# Indel:  " << indel_hist << endl;
	    cout << "Alphabet size: "<< types_vector.size()<< endl;
	    //cout << "List of alphabet characters: \n" ;
	    //for(int i=0;i<types_vector.size() ;i++){
	    //	cout<< types_vector[i].c_str() <<", ";
	    //}
	    //cout << endl;

	}
	/// reading scores
	int i=0; 
	
	
	while(getline(in_file,zeile)){
	    //cout << zeile << endl;
	    words.clear();
	    fsplit(zeile, " ", words);
	    if(words.size()!=(type_no-2-i)){
		error_flag=1;
		break;
	    }
	    in_dist[i][i]=0;
	    for(int j=0;j<words.size();j++){
		in_dist[i][j+i+1]=atoi(words[j].c_str());
		in_dist[j+i+1][i]=atoi(words[j].c_str());
		
		//in_dist[j][j]=0;
	    }
	    i++;
	    if(i==type_no-1){
		break;
	    }
	}
	// adding the cost for the special symbol
	for(int i=0;i<type_no;i++){	    
	    in_dist[i][type_no-1]=BigValue2;
	}
	for(int j=0;j<type_no;j++){
	    in_dist[type_no-1][j]=BigValue2;
	}
	in_dist[type_no-1][type_no-1]=0;

	if((type_no-1>i+1)||(error_flag)){
	    cerr << "Error: no of types is not consistent with distance matrix\n";
	    in_file.close();
	    exit(-1);
	}
	else{
	    if(show_cost_file_flag){
		cout << "# The cost matrix: " << endl;
		for(int i=0;i<type_no;i++){	    
		    for(int j=0;j<type_no;j++){
			cout << in_dist[i][j] << " ";
		    }
		    cout << endl;
		}
	    }
	    in_file.close();
	}
}

int ThroughputAlignment::check_traingularity_of_cost(){
    
    //return(1);
    for(int i=0;i<type_no;i++){	    
	for(int j=0;j<type_no;j++){
	    for(int k=0;k<type_no;k++){
		if(in_dist[i][k]+in_dist[k][j]<in_dist[i][j]){		    
		    cout << "The cost matrix does not satisfy the triangular equality\n";
		    cout << "dist(" << i <<","<<k<< ")+" << "dist(" << k <<","<<j<< ")<" << "dist(" << i <<","<<j<< ")=>"; 
		    cout << in_dist[i][k] << "+" << in_dist[k][j] <<"<" << in_dist[i][j] << "\n";
		    return(-1);
		}
	    }
	}
	
    }
    return(1);
}

void fsplit (string line, const char* sep, vector<string> & words) {
		
	string::size_type a = 0, e;
	while ((a = line.find_first_not_of (sep, a)) != string::npos) {
		e = line.find_first_of (sep, a);
		if (e != string::npos) {
			words.push_back (line.substr (a, e-a));
			a = e + 1;
		}
		else {
			words.push_back (line.substr (a));
			break;
		}
	}
}
