// TestingDupHist.cpp: implementation of the TestingDupHist class.
//
//////////////////////////////////////////////////////////////////////

//#include "stdafx.h"
#include "TestingDupHist.h"


#include <stdio.h>
#include <string.h>
#include <iostream>
#include <stdlib.h>
#include <vector>
//#include <windows.h>


#ifndef max
#define max(a, b)  ((a)>=(b) ? (a) : (b))
#endif
#ifndef min
#define min(a, b)  ((a)<=(b) ? (a) : (b))
#endif

using namespace std;

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

TestingDupHist::TestingDupHist()
{
	
    alphadep_flag=0;
/*
		// Intialize sequence array
	for (int i=0;i<seqlength;i++){
		sequence1[i]=0;
		sequence2[i]=0;
	}
*/
	reverseflag=0;
	fixed_dup_cost=0;
	indel_hist=0;
	indel_align=0;

	Amatrix=NULL;
	Bmatrix=NULL;
	histmatrix_out=NULL;
	histmatrix_in=NULL;
	show_alignment_flag=0;
	show_breakeven_flag=0;

	inverse_ptr=NULL;
	histmatrix_outRLE=NULL;
	RLE_sequence=NULL;
}

TestingDupHist::~TestingDupHist()
{ 
    
    
    if(histmatrix_out!=NULL){
	
	for(int x=0;x<seqlength+1;x++){
	    for(int y=0;y<seqlength+1;y++){		
		//cout << "cleaaaaaaaaaaaaaaaaaaaning\n";
		delete [] histmatrix_out[x][y];
		
	    }
	}
	for(int x=0;x<seqlength+1;x++){
	    delete[] histmatrix_out[x];	   
	}
	delete[] histmatrix_out;
    }
   
    if(histmatrix_in!=NULL){
	for(int x=0;x<seqlength+1;x++){
	    for(int y=0;y<seqlength+1;y++){	
		
		delete []  histmatrix_in[x][y];
	    }
	}
	for(int x=0;x<seqlength+1;x++){
	    delete[] histmatrix_in[x];	   
	}
	delete[] histmatrix_in;
    }

    // cleaning up
    if(inverse_ptr!=NULL){
	delete[] inverse_ptr;inverse_ptr=NULL;
    }
    if(histmatrix_outRLE!=NULL){
	for(int x=0;x<seqlengthRLE;x++){
	    for(int y=0;y<seqlengthRLE;y++){	
		delete [] histmatrix_outRLE[x][y];
	    }
	}		
	for(int x=0;x<seqlengthRLE;x++){
	    delete[] histmatrix_outRLE[x];	   
	}	
	delete[] histmatrix_outRLE; 
	histmatrix_outRLE=NULL;
    }    
    if(RLE_sequence!=NULL){
	delete[] RLE_sequence;RLE_sequence=NULL;
    }
    ///end cleaning up


/*

*/
    return;

}
void TestingDupHist::FreeMEMORY()
{ 
   
}

void TestingDupHist::allocate_array(int seq_len, int alphadep_flag){
   
    
    if(alphadep_flag==0){
	histmatrix_out=new int**[seq_len+1]; // for x
	for(int x=0;x<seq_len+1;x++){
	    histmatrix_out[x]=new int*[seq_len+1];	   
	}
	for(int x=0;x<seq_len+1;x++){
	    for(int y=0;y<seq_len+1;y++){	
		histmatrix_out[x][y]=new int[1];
	    }
	}
	
	for(int x=0;x<seq_len+1;x++){
	    for(int y=0;y<seq_len+1;y++){
		histmatrix_out[x][y][0]=BigValue;				
	    }	    
	}
	
    }
    else{	
	//cout << "halooooooooooooooooooo "<< inputalphabet << endl;
	histmatrix_out=new int**[seq_len+1]; // for x
	for(int x=0;x<seq_len+1;x++){
	    histmatrix_out[x]=new int*[seq_len+1];	   
	}
	for(int x=0;x<seq_len+1;x++){
	    for(int y=0;y<seq_len+1;y++){	
		histmatrix_out[x][y]=new int[inputalphabet+1];
	    }
	}
	
	for(int x=0;x<seq_len+1;x++){
	    for(int y=0;y<seq_len+1;y++){
		for(int z=0;z<inputalphabet+1;z++){
		    histmatrix_out[x][y][z]=BigValue;		
		}
	    }
	}
	histmatrix_in=new int**[seq_len+1]; // for x
	for(int x=0;x<seq_len+1;x++){
	    histmatrix_in[x]=new int*[seq_len+1];	   
	}
	for(int x=0;x<seq_len+1;x++){
	    for(int y=0;y<seq_len+1;y++){	
		histmatrix_in[x][y]=new int[inputalphabet+1];
	    }
	}
	
	for(int x=0;x<seq_len+1;x++){
	    for(int y=0;y<seq_len+1;y++){
		for(int z=0;z<inputalphabet+1;z++){
		    histmatrix_in[x][y][z]=BigValue;		
		}
	    }
	}
	
    }
    
}
//////////////////////////////////////////////////////////////////
void TestingDupHist::Intialize(int** in_dist)
{
	
   int i,j,k;	   
   dist=in_dist;
 
}
///////////////////////  My Code /////////////////////////////////

int TestingDupHist::SymmBensonDong(){

    int i,j,k,sl,sx;
	int tempval1=(int)0;
	int tempval2=(int)0;
	int tempval3=(int)BigValue;
	int tempval31=(int)0;
	int tempval32=(int)0;
	int tempval33=(int)0;
//	char originchar;

	for(int x=0;x<seqlength;x++){
	    for(int y=0;y<seqlength;y++){
		histmatrix_out[x][y][0]=BigValue;		
	    }
	}

	//cout << "Symmmmmmmmmmmmmmmmmmmm Dong: " << seqlength << endl;
	//FILE* fptr;
	//fptr=fopen("TempXFGresult.txt","w");

	for (i=1;i<=seqlength;i++){ // i is the interval length
		for(j=0;j+i-1<seqlength;j++){
			// the interval boundaries are [j..j+i-1]
			if(i==1){ //initialization
			    
			    for(k=j;k<j+i;k++){
				histmatrix_out[k][k][0]=(int)0;//dist[sequence1[k]][sx];
				
			    }
				
					
			}
			else if(i==2){
				sx=0;
				histmatrix_out[j][j+i-1][0]=dist[sequence1[j+i-1]][sequence1[j]];
				histmatrix_out[j+i-1][j][0]=dist[sequence1[j]][sequence1[j+i-1]];
//				fprintf(fptr,"[%d..%d] score: %2f\n",j,j+i-1,histmatrix_out[j][j+i-1][sx]);
			}
			else{
				
			
			
				sx=0;

				// our second correction for unsymmetric case
				// we have to consider to cases the interval is originated from
				// the leftmost or the right most.
				// case 1
					
				for(k=j+1;k<j+i-1;k++){
					// originated from the left side
					tempval1=histmatrix_out[j][k][0]+histmatrix_out[k][j+i-1][0];
					histmatrix_out[j][j+i-1][0]=min(tempval1,histmatrix_out[j][j+i-1][0]);
					// originated from the right side
					tempval1=histmatrix_out[k][j][0]+histmatrix_out[j+i-1][k][0];
					histmatrix_out[j+i-1][j][0]=min(tempval1,histmatrix_out[j+i-1][j][0]);

				}


				for(sl=j;sl<j+i-1;sl++){
					// originated from the left most
					tempval2=histmatrix_out[j][sl][0]+histmatrix_out[j+i-1][sl+1][0]+
						dist[sequence1[j+i-1]][sequence1[j]];
					histmatrix_out[j][j+i-1][0]=min(tempval2,histmatrix_out[j][j+i-1][0]);
					// originated from rightmost
					tempval3=histmatrix_out[j][sl][0]+histmatrix_out[j+i-1][sl+1][0]+
						dist[sequence1[j]][sequence1[j+i-1]];
					histmatrix_out[j+i-1][j][0]=min(tempval3,histmatrix_out[j+i-1][j][0]);

				}

				// end of second correction
				
				//fprintf(fptr,"[%d..%d] score left: %d, score right %d,\n",j,j+i-1,histmatrix_out[j][j+i-1][0],histmatrix_out[j+i-1][j][0]);
		
			}// else 
			//fprintf(fptr,"[%d..%d] score left: %d, score right %d,\n",j,j+i-1,histmatrix_out[j][j+i-1][0],histmatrix_out[j+i-1][j][0]);
		}
	}

	// Final Full Tree
//	fprintf(fptr,"An optimal final trees \n");
//	fprintf(fptr,"[%d..%d] score: %.2f \n",0,(seqlength-1),min(histmatrix_out[0][seqlength-1][0],histmatrix_out[seqlength-1][0][0]));
	
	
	bensonleft=histmatrix_out[0][seqlength-1][0];
	bensonright=histmatrix_out[seqlength-1][0][0];

	// added to consider length of duplication
	
	for(int x=0;x<seqlength-1;x++){
		for(int y=x+1;y<seqlength;y++){		
		    histmatrix_out[x][y][0]=histmatrix_out[x][y][0]+(abs((x-y))*fixed_dup_cost);
		    histmatrix_out[y][x][0]=histmatrix_out[y][x][0]+(abs((x-y))*fixed_dup_cost);		    
		    //    fprintf(fptr,"[%d..%d] score left: %d, score right %d,\n",x,y,histmatrix_out[x][y][0],histmatrix_out[y][x][0]);
		
		}
	}
//	fclose(fptr);
	//// finish addition	
	return(min(histmatrix_out[0][seqlength-1][0],histmatrix_out[seqlength-1][0][0]));
}


////////////////////////////////////////////////////////////////////
/////////////      and N^3 Implementation         //////////////////
int TestingDupHist::ModAlignment(TestingDupHist* Sptr,TestingDupHist* Rptr){
	
	int i,j,ls,lr,ts,tr;
	int seqSlength,seqRlength;
	int tempval1=BigValue;
	int tempval2=BigValue;

	int* seqS; 
	int* seqR;
	
	//Sptr->histmatrix_out;
	//histmatrix2=(int*)Rptr->histmatrix_out;
	seqS=Sptr->sequence1; 
	seqR=Rptr->sequence1;

	seqSlength=Sptr->seqlength;
	seqRlength=Rptr->seqlength;

	// initialize arrays
	if(1){
	    Amatrix=new int*[seqSlength+1];
	    Bmatrix=new int*[seqSlength+1];
	    for(int y=0;y<seqSlength+1;y++){
		Amatrix[y]=new int[seqRlength+1];
		Bmatrix[y]=new int[seqRlength+1];
	    }
	    for (i=0;i<seqSlength;i++){
		for(j=0;j<seqRlength;j++){
		    Amatrix[i][j]=(int)BigValue;
		    Bmatrix[i][j]=(int)BigValue;
		}
	    }
	}
	/// end intialize create array

	Amatrix[0][0]=min(dist[seqS[0]][seqR[0]],dist[seqR[0]][seqS[0]]);
	Bmatrix[0][0]=min(dist[seqS[0]][seqR[0]],dist[seqR[0]][seqS[0]]);
	
	for (i=0;i<seqSlength;i++){
		for(j=0;j<seqRlength;j++){
			
			
			if((i==0)&&(j==1)){
				int kk=0;
			}			
			// The match case
			if((i>0)&&(j>0))
			    Amatrix[i][j]=Amatrix[i-1][j-1]+min(dist[seqS[i]][seqR[j]],dist[seqR[j]][seqS[i]]);
			
			// leftside originated duplication
			for(ls=0;ls<i;ls++){				
				Amatrix[i][j]=min(Amatrix[ls][j]+Sptr->histmatrix_out[ls][i][0],Amatrix[i][j]);				
			}
			
			for(lr=0;lr<j;lr++){
				Amatrix[i][j]=min(Amatrix[i][lr]+Rptr->histmatrix_out[lr][j][0],Amatrix[i][j]);			
			}
			// Filling Bmatrix w.r.t. j
			
			if((i==0)&&(j==1)){
				int kk=0;
			}
			for(tr=0;tr<j;tr++){
				//if(i!=0){
					Bmatrix[i][j]=min(Amatrix[i][tr]+Rptr->histmatrix_out[j][tr+1][0],Bmatrix[i][j]);
				//}
				//else{
				//	Bmatrix[i][j]=min(Bmatrix[i][tr]+Rptr->histmatrix_out[j][0][0],Bmatrix[i][j]);
				//}
			}

			// rightside originated duplication
			//debug
			if((i==seqSlength-1)&&(j==seqRlength-1)){
				int kk=0;
			}
			if(i==0){
				Amatrix[i][j]=min(Rptr->histmatrix_out[j][0][0]+dist[seqS[i]][seqR[j]],Amatrix[i][j]);
				Amatrix[i][j]=min(Rptr->histmatrix_out[j][0][0]+dist[seqR[j]][seqS[i]],Amatrix[i][j]);
				
			}
			else if(j==0){
				Amatrix[i][j]=min(Sptr->histmatrix_out[i][0][0]+dist[seqS[i]][seqR[j]],Amatrix[i][j]);
				Amatrix[i][j]=min(Sptr->histmatrix_out[i][0][0]+dist[seqR[j]][seqS[i]],Amatrix[i][j]);

			}
			else{
				for(ts=0;ts<i;ts++){
				    Amatrix[i][j]=min(Bmatrix[ts][j]+Sptr->histmatrix_out[i][ts+1][0]+dist[seqS[i]][seqR[j]],Amatrix[i][j]);
				    Amatrix[i][j]=min(Bmatrix[ts][j]+Sptr->histmatrix_out[i][ts+1][0]+dist[seqR[j]][seqS[i]],Amatrix[i][j]); 	   	   
				}
				Amatrix[i][j]=min(Rptr->histmatrix_out[j][0][0]+Sptr->histmatrix_out[i][0][0]+dist[seqS[i]][seqR[j]],Amatrix[i][j]);
				Amatrix[i][j]=min(Rptr->histmatrix_out[j][0][0]+Sptr->histmatrix_out[i][0][0]+dist[seqR[j]][seqS[i]],Amatrix[i][j]);	   		
			}
			

		}
	}
	
	if(show_alignment_flag==1){
	    show_alignment(Sptr,Rptr,types_vector_ptr);
	}
	
	int optimal_Score=Amatrix[seqSlength-1][seqRlength-1];
	if(Amatrix!=NULL){
	    for(int y=0;y<seqSlength+1;y++){
		delete [] Amatrix[y];		
	    }
	    delete [] Amatrix;
	}
	
	if(Bmatrix!=NULL){
	    for(int y=0;y<seqSlength+1;y++){
		delete [] Bmatrix[y];	    
	    }
	    delete [] Bmatrix;
	}


	return(optimal_Score);
	

}
///////////////////////////////////////////////////////////////////77
void TestingDupHist::getsequence(int* in_sequence,int seqlength)
{
	
	int i=0;
	int count=0;

	for (i=0;i<seqlength;i++){
	    sequence1[i]=in_sequence[i];
	}
	//cout << "In getting sequence: " << seqlength << endl;
	/*
	for (i=0;i<seqlength;i++){
	    cout<< (int) sequence1[i]<<" ";
	}
	cout << endl;
	*/


	// to test the reverse
	if(reverseflag==1){
		cout << "NOTE: The reverse is ON"<<"\n";
		char temp;
		int j;
		j=seqlength-1;
		for (i=0;i<seqlength/2;i++){
			temp=sequence1[i];
			sequence1[i]=sequence1[j]; 
			sequence1[j]=temp;
				
			j--;
		}
	}
	// end the test

 	
}


/////////////      and N^4 Implementation         //////////////////
int TestingDupHist::AlignmentWithoutRightDup(TestingDupHist* Sptr,TestingDupHist* Rptr){
	
	int i,j,ls,lr,ts,tr;
	int seqSlength,seqRlength;
	int tempval1=BigValue;
	int tempval2=BigValue;

	int* seqS; 
	int* seqR;
	
	//Sptr->histmatrix_out;
	//histmatrix2=(double*)Rptr->histmatrix_out;
	seqS=Sptr->sequence1; 
	seqR=Rptr->sequence1;

	seqSlength=Sptr->seqlength;
	seqRlength=Rptr->seqlength;


		// initialize arrays
	if(1){
	    Amatrix=new int*[seqSlength+1];
	    Bmatrix=new int*[seqSlength+1];
	    for(int y=0;y<seqSlength+1;y++){
		Amatrix[y]=new int[seqRlength+1];
		Bmatrix[y]=new int[seqRlength+1];
	    }
	    for (i=0;i<seqSlength;i++){
		for(j=0;j<seqRlength;j++){
		    Amatrix[i][j]=(int)BigValue;
		    Bmatrix[i][j]=(int)BigValue;
		}
	    }
	}
	/// end intialize create array


/*
	// reverse the seqeunce
    // Reversing 
	char temp;
	j=seqSlength-1;
	for (i=0;i<seqSlength/2;i++){
		temp=seqS[i];
		seqS[i]=seqS[j]; 
		seqS[j]=temp;
		
		j--;
	}

	j=seqRlength-1;
	for (i=0;i<seqRlength/2;i++){
		temp=seqR[i];
		seqR[i]=seqR[j]; 
		seqR[j]=temp;		
		j--;
	}
	*/
	// End reversinng
    ////
	Amatrix[0][0]=min(dist[seqS[0]][seqR[0]],dist[seqR[0]][seqS[0]]);
	for (i=0;i<seqSlength;i++){
		for(j=0;j<seqRlength;j++){
			
			// The match case
			if((i>0)&&(j>0))
			Amatrix[i][j]=Amatrix[i-1][j-1]+min(dist[seqS[i]][seqR[j]],dist[seqR[j]][seqS[i]]);
					//debug
			if((i==0)&&(j==seqRlength-1)){
				int kk=0;
			}
			
			// leftside originated duplication
			for(ls=0;ls<i;ls++){
				Amatrix[i][j]=min(Amatrix[ls][j]+Sptr->histmatrix_out[ls][i][0],Amatrix[i][j]);
			}
			
			for(lr=0;lr<j;lr++){
				Amatrix[i][j]=min(Amatrix[i][lr]+Rptr->histmatrix_out[lr][j][0],Amatrix[i][j]);
			}
			// rightside originated duplication
			// the naive part
/*
			//debug
			if((i==seqSlength-1)&&(j==seqRlength-1)){
				int kk=0;
			}
			if(i==0){
				for(tr=0;tr<j;tr++){
					tempval1=dist[seqS[i]][seqR[j]]+Rptr->histmatrix_out[j][0][0];
					Amatrix[i][j]=min(tempval1,Amatrix[i][j]);
					tempval1=dist[seqR[j]][seqS[i]]+Rptr->histmatrix_out[j][0][0];
					Amatrix[i][j]=min(tempval1,Amatrix[i][j]);
				}
			}
			else if(j==0){
				for(ts=0;ts<i;ts++){
					tempval1=dist[seqS[i]][seqR[j]]+Sptr->histmatrix_out[i][0][0];
					Amatrix[i][j]=min(tempval1,Amatrix[i][j]);
					tempval1=dist[seqR[j]][seqS[i]]+Sptr->histmatrix_out[i][0][0];
					Amatrix[i][j]=min(tempval1,Amatrix[i][j]);
					
				}
			}
			else{
				for(ts=0;ts<=i;ts++){
					for(tr=0;tr<=j;tr++){
						if((tr==0)||(ts==0)){
							tempval1=dist[seqS[i]][seqR[j]]+Rptr->histmatrix_out[j][0][0]+Sptr->histmatrix_out[i][0][0];
							Amatrix[i][j]=min(tempval1,Amatrix[i][j]);
							tempval1=dist[seqR[j]][seqS[i]]+Rptr->histmatrix_out[j][0][0]+Sptr->histmatrix_out[i][0][0];
							Amatrix[i][j]=min(tempval1,Amatrix[i][j]);
						}

						else{
							tempval1=dist[seqS[i]][seqR[j]]+Rptr->histmatrix_out[j][tr][0]+Sptr->histmatrix_out[i][ts][0]+Amatrix[ts-1][tr-1];
							Amatrix[i][j]=min(tempval1,Amatrix[i][j]);
							tempval1=dist[seqR[j]][seqS[i]]+Rptr->histmatrix_out[j][tr][0]+Sptr->histmatrix_out[i][ts][0]+Amatrix[ts-1][tr-1];
							Amatrix[i][j]=min(tempval1,Amatrix[i][j]);
						}
					}
				}
			}
*/			

		}
	}
	if(show_alignment_flag==1){
	    show_alignment(Sptr,Rptr,types_vector_ptr);
	}
	int optimal_Score=Amatrix[seqSlength-1][seqRlength-1];
	if(Amatrix!=NULL){
	    for(int y=0;y<seqSlength+1;y++){
		delete [] Amatrix[y];		
	    }
	    delete [] Amatrix;
	}
	if(Bmatrix!=NULL){
	    for(int y=0;y<seqSlength+1;y++){
		delete [] Bmatrix[y];	    
	    }
	    delete [] Bmatrix;
	}

	return(optimal_Score);

}
/////////////////////////////////////////////////////////////////////
void TestingDupHist::intialize_type2distance(){

}

void TestingDupHist::intialize_type3distance(){

}

void TestingDupHist::intialize_typeBouzekridistance(){

	// 0 -->'0'
	dist[0][0]=0;dist[0][0]=0;           // a,a
	dist[0][1]=4;dist[1][0]=4;           // a,b     
	dist[0][2]=4;dist[2][0]=4;           // a,c
	dist[0][3]=4;dist[3][0]=4;           // a,d 
	dist[0][4]=4;dist[4][0]=4;           // a,e
	dist[0][5]=4;dist[5][0]=4;           // a,e
	dist[0][6]=4;dist[6][0]=4;           // a,e
	dist[0][7]=4;dist[7][0]=4;           // a,e
	dist[0][8]=4;dist[8][0]=4;           // a,e
	dist[0][9]=4;dist[9][0]=4;           // a,e

	// 1--> '1'
	dist[1][1]=0;dist[1][1]=0;          // b,b
	dist[1][2]=1;dist[2][1]=1;          // b,c     
	dist[1][3]=1;dist[3][1]=1;          // b,d
	dist[1][4]=1;dist[4][1]=1;          // b,e 
	dist[1][5]=2;dist[5][1]=2;          // b,b
	dist[1][6]=2;dist[6][1]=2;          // b,c     
	dist[1][7]=3;dist[7][1]=3;          // b,d
	dist[1][8]=3;dist[8][1]=3;          // b,e 
	dist[1][9]=4;dist[9][1]=4;          // b,e 

	
	// 2--> '1a'
	dist[2][2]=0;dist[2][2]=0;          // c,c
	dist[2][3]=2;dist[3][2]=2;          // c,d     
	dist[2][4]=2;dist[4][2]=2;          // c,e
	dist[2][5]=1;dist[5][2]=1;          // c,c
	dist[2][6]=3;dist[6][2]=3;          // c,d     
	dist[2][7]=2;dist[7][2]=2;          // c,e
	dist[2][8]=2;dist[8][2]=2;          // c,e
	dist[2][9]=4;dist[9][2]=4;          // c,e
	
	// 3--> '2'
	dist[3][3]=0;dist[3][3]=0;          // d,d
	dist[3][4]=2;dist[4][3]=2;          // d,e     
	dist[3][5]=3;dist[5][3]=3;          // d,d
	dist[3][6]=1;dist[6][3]=1;          // d,e     
	dist[3][7]=2;dist[7][3]=2;          // d,e     
	dist[3][8]=2;dist[8][3]=2;          // d,d
	dist[3][9]=4;dist[9][3]=4;          // d,e     

	// 4--> '3'
	dist[4][4]=0;dist[4][4]=0;          // d,d
	dist[4][5]=1;dist[5][4]=1;          // d,d
	dist[4][6]=1;dist[6][4]=1;          // d,d
	dist[4][7]=2;dist[7][4]=2;          // d,d
	dist[4][8]=2;dist[8][4]=2;          // d,d
	dist[4][9]=4;dist[9][4]=4;          // d,d

	// 5--> '3a'
	dist[5][5]=0;dist[5][5]=0;          // d,d
	dist[5][6]=2;dist[6][5]=2;          // d,d
	dist[5][7]=1;dist[7][5]=1;          // d,d
	dist[5][8]=1;dist[8][5]=1;          // d,d
	dist[5][9]=4;dist[9][5]=4;          // d,d

	// 6--> '4'
	dist[6][6]=0;dist[6][6]=0;          // d,d
	dist[6][7]=1;dist[7][6]=1;          // d,d
	dist[6][8]=1;dist[8][6]=1;          // d,d
	dist[6][9]=4;dist[9][6]=4;          // d,d

	// 7 --> '4a'
	dist[7][7]=0;dist[7][7]=0;          // d,d
	dist[7][8]=1;dist[8][7]=1;          // d,d
	dist[7][9]=2;dist[9][7]=2;          // d,d

	// 8 --> '4*'
	dist[8][8]=0;dist[8][8]=0;          // d,d
	dist[8][9]=2;dist[9][8]=2;          // d,d
	for(int i=0;i<10;i++){
		for(int j=0;j<10;j++){
			dist[i][j]=1*dist[i][j];
		}		
	}
}

void TestingDupHist::intialize_typeRtest(){

}

//////////////////////////////////////////////////////////////////////
/// This function constructs the duplication histroy
/// depending on the alphabet size
/////////////////////////////////////////////////////////////////
int TestingDupHist::constructhistory(char* filename){

    int i,j,k,sl,sx;
    int tempval1=0;
    int tempval2=0;
    int tempval3=BigValue;
    int tempval31=0;
    int tempval32=0;
    int tempval33=0;
    char originchar;


    for (i=1;i<seqlength;i++){ // i is the interval length
	for(j=0;j+i-1<seqlength;j++){
	    // the interval boundaries are [j..j+i-1]
	    if(i==1){ //initialization
		for(sx=0;sx<Alphasize;sx++){
		    for(k=j;k<j+i;k++){
			histmatrix_out[k][k][sx]=dist[sequence1[k]][sx];
			histmatrix_in[k][k][sx]=dist[sequence1[k]][sx];
			
		    }
		}
		
	    }
	    else{
		
		// Full Tree
		for(sl=j;sl<j+i;sl++){
		    //debug dd
		    if((j==4)&&(i==4)){
			int kk=0;
		    }
		    if(sl==j){
			tempval31=histmatrix_out[sl+1][j+i-1][sequence1[sl]];
		    }
		    else if(sl==j+i-1){			
			tempval31=histmatrix_out[j][sl-1][sequence1[sl]];
		    }
		    else{
			tempval31=histmatrix_out[j][sl-1][sequence1[sl]]+histmatrix_out[sl+1][j+i-1][sequence1[sl]];
		    }
		    histmatrix_in[j][j+i-1][sequence1[sl]]=
			min(histmatrix_in[j][j+i-1][sequence1[sl]],tempval31);
		    if(tempval3>histmatrix_in[j][j+i-1][sequence1[sl]]){
			tempval3=histmatrix_in[j][j+i-1][sequence1[sl]];
			originchar=sequence1[sl];
		    }
		}
		// External Rooted Tree
		
		for(sx=0;sx<Alphasize;sx++){
		    
		    if((j==1)&&(i==3)&&(sx==2)){
			int kk=0;
		    }
		    for(k=j;k<j+i-1;k++){
			tempval1=histmatrix_out[j][k][sx]+histmatrix_out[k+1][j+i-1][sx];
			histmatrix_out[j][j+i-1][sx]=min(tempval1,histmatrix_out[j][j+i-1][sx]);
		    }
		    for(sl=j;sl<j+i;sl++){
			tempval2=histmatrix_in[j][j+i-1][sequence1[sl]]+dist[sequence1[sl]][sx];
			histmatrix_out[j][j+i-1][sx]=min(tempval2,histmatrix_out[j][j+i-1][sx]);
			
		    }
		    //tempval1=min(tempval1,tempval2);		    		    
		}
		
		
//				fprintf(fptr,"[%d..%d] [%d] %2f\n",j,j+i-1,originchar,tempval3);
		
		
				tempval3=BigValue;
			}// else 
		
		}
	}

	// Final Full Tree
		j=0;
		i=seqlength;
		for(sl=j;sl<j+i;sl++){
//			fprintf(fptr,"The final trees \n");
			if(sl==j){
				tempval31=histmatrix_out[sl+1][j+i-1][sequence1[sl]];
			}
			else if(sl==j+i-1){
					tempval31=histmatrix_out[j][sl-1][sequence1[sl]];
			}
			else{
				tempval31=histmatrix_out[j][sl-1][sequence1[sl]]+histmatrix_out[sl+1][j+i-1][sequence1[sl]];
			}
				histmatrix_in[j][j+i-1][sequence1[sl]]=
					min(histmatrix_in[j][j+i-1][sequence1[sl]],tempval31);
				if(tempval3>histmatrix_in[j][j+i-1][sequence1[sl]]){
					tempval3=histmatrix_in[j][j+i-1][sequence1[sl]];
					originchar=sequence1[sl];
					tempval31=sl;
				}
//			fprintf(fptr,"[%d..%d] [%d] %.2f %d\n",j,j+i-1,sequence1[sl],histmatrix_in[j][j+i-1][sequence1[sl]],sl);
		}
//		fprintf(fptr,"An optimal final trees \n");
//		fprintf(fptr,"[%d..%d] [%d] %.2f %d\n",j,j+i-1,originchar,tempval3,tempval31);
	
//	fclose(fptr);
	oursleft=histmatrix_in[j][j+i-1][sequence1[0]];
	oursright=histmatrix_in[j][j+i-1][sequence1[seqlength-1]];

	// added to consider length of duplication
	
	for(int x=0;x<seqlength;x++){
		for(int y=0;y<seqlength;y++){
			for(sx=0;sx<Alphasize;sx++){
			    //    histmatrix_in[x][y][sx]=histmatrix_in[x][y][sx]+(abs((x-y)*fixed_dup_cost));
			    //histmatrix_out[x][y][sx]=histmatrix_out[x][y][sx]+(abs((x-y)*fixed_dup_cost));	
			    
			}
		}
	}
	//// finish addition
	//tempval3=tempval3+abs((seqlength-1)*fixed_dup_cost);
	return(tempval3);
}


void TestingDupHist::show_alignment(TestingDupHist* Sptr,TestingDupHist* Rptr,vector <string> types_vector_ptr ){

    	int i,j,ls,lr,ts,tr;
	int seqSlength,seqRlength;
	int tempval1=BigValue;
	int tempval2=BigValue;

	int* seqS; 
	int* seqR;

	cout << "\n#----------- The operations ------------\n";	
	
	//Sptr->histmatrix_out;
	//histmatrix2=(int*)Rptr->histmatrix_out;
	seqS=Sptr->sequence1; 
	seqR=Rptr->sequence1;

	seqSlength=Sptr->seqlength;
	seqRlength=Rptr->seqlength;
	string aligned_seq_S;
	string aligned_seq_R;
	anchor_pair anchor_pair_obj;
	
	vector <anchor_pair> anchor_list;
	//vector <string> types_vector_ptr;
	int anchor_flag=0;
	i=seqSlength-1;
	j=seqRlength-1;
	int found_flag=0;
	while((i!=-1)||(j!=-1)){	
	    
	    found_flag=0;
	    
	    // The match case
	    
	    if(((i==0)&&(j==0))&&(Amatrix[i][j]==min(dist[seqS[i]][seqR[j]],dist[seqR[j]][seqS[i]]))){
		cout << "Match             : " <<  "("<< i <<", "<<j<<")"<< " score: "<<Amatrix[i][j] << endl;
		i--;
		j--;	
		anchor_flag=1;
		found_flag=1;
		continue;
	    }
	    
	    if((i>0)&&(j>0)){
		if(Amatrix[i][j]==(Amatrix[i-1][j-1]+min(dist[seqS[i]][seqR[j]],dist[seqR[j]][seqS[i]]))){		
		    anchor_flag=1;
		    cout << "Match             : " <<  "("<< i <<", "<<j<<")"<< " score: "<<Amatrix[i][j] << endl;
		    i--;
		    j--;		
		    anchor_flag=1;
		    found_flag=1;
		    continue;
		}
	    }
	    
	    // The left duplication cases
	    for (ls=0;ls<i;ls++){
		if(ls<i){
		    if(Amatrix[i][j]==(Amatrix[ls][j]+Sptr->histmatrix_out[ls][i][0])){			
			anchor_flag=1;
			//cout << "Left dup. in S:     " << ls << "-->"<< i <<" , "<<j<<")"<< " score: "<<Amatrix[i][j] << endl;
			cout << "Left dup. in S    : " << "[" << ls << ".."<< i   << "]" << " score: "<<Amatrix[i][j] << endl;
			show_dup_scenario(Sptr,types_vector_ptr, ls, i,ls);
			i=ls;
			found_flag=1;
			break;
		    }
		}
	    }
	    
	    if(anchor_flag==1){
		anchor_flag=0;
		continue;
	    }
	    
	    for(lr=0;lr<j;lr++){
		if(lr<j){
		    if(Amatrix[i][j]==(Amatrix[i][lr]+Rptr->histmatrix_out[lr][j][0])){
			anchor_flag=1;
			//cout << "Left dup. in R:     " << lr <<" ("<< i <<" , "<<j<<")"<<  " score: "<<Amatrix[i][j] << endl;
			cout << "Left dup. in R    : "  << "[" << lr << ".."<<j << "] score: "<<Amatrix[i][j] << endl;
			show_dup_scenario(Rptr,types_vector_ptr, lr, j,lr);
			found_flag=1;
			j=lr;			    
			break;
		    }
		}
	    }
	    if(anchor_flag==1){
		anchor_flag=0;
		continue;
	    }

	    ls=0;
	    lr=0;
	    if(Amatrix[i][j]==(Sptr->histmatrix_out[i][0][0]+Rptr->histmatrix_out[j][0][0]
			       +min(dist[seqS[i]][seqR[j]],dist[seqR[j]][seqS[i]]))){	
		//cout << "Right dup. in S, R: " << ls <<" " <<lr << "("<< i <<" , "<<j<<")"<<  " score: "<<Amatrix[i][j] << endl;
		cout << "Match             : " << "("<< i <<", "<<j<<")" << " score: "<< Amatrix[i][j] << endl;
		if((ls+1)!=i){
		    cout << "Right dup. in S   : "  << "["<<ls+1 << ".." << i      << "]"<< " score: "<< Amatrix[i][j] << endl;
		    show_dup_scenario(Sptr,types_vector_ptr, ls+1, i,i);	      
		}
		if((lr+1)!=j){
		    cout << "Right dup. in R   : "  << "[" <<lr+1 << ".." << j     << "]"<< " score: "<< Amatrix[i][j] << endl;
		    show_dup_scenario(Rptr,types_vector_ptr, lr+1, j,j);
		}
		anchor_flag=1;
		found_flag=1;
		
		break;			   
	    }
	    
	    for (ls=0;ls<i;ls++){	
		for(lr=0;lr<j;lr++){	    					    	       		  	
		    if((ls<i)&&(lr<j)){
			if(Amatrix[i][j]==(Amatrix[ls][lr]+Sptr->histmatrix_out[i][ls+1][0]+Rptr->histmatrix_out[j][lr+1][0])
			   +min(dist[seqS[i]][seqR[j]],dist[seqR[j]][seqS[i]])){

			    //cout << "Right dup. in S, R: " << ls <<" " <<lr << "("<< i <<" , "<<j<<")"<<  " score: "<<Amatrix[i][j] << endl;
			    cout << "Match             : " << "("<< i <<", "<<j<<")" << " score: "<<Amatrix[i][j] << endl;
			    //cout << "Right dup. in S, R: " << ls+1 << "->" << i << " & " << lr+1 << "->" << j << " score: "<<Amatrix[i][j] << endl;
			    if((ls+1)!=i){
				cout << "Right dup. in S   : "  << "["<< ls+1 << ".." << i <<"]"     << " score: "<< Amatrix[i][j] << endl;
				show_dup_scenario(Sptr,types_vector_ptr, ls+1, i,i);
			    }
			    if((lr+1)!=j){
				cout << "Right dup. in R   : "  << "["<< lr+1 << ".." << j <<"]"   << " score: "<< Amatrix[i][j] << endl;
				show_dup_scenario(Rptr,types_vector_ptr, lr+1, j,j);
			    }
			    found_flag=1;
			    i=ls;
			    j=lr;			    
			    anchor_flag=1;
			    break;
			}
		    }
		} // end inner for loop
		if(anchor_flag==1){
		    anchor_flag=0;
		    break;
		}
	    } // end outter for
	    if(found_flag!=1){
		cout << "Not found_flag "<< i << " "<<j << endl;
		break;
	    }
	    
	}// end while loop
	
	cout << endl;
	cout <<"Note the unit at position 0 is the artificial unit \""<<types_vector_ptr[seqR[0]] << "\""<<endl;
	//cout<< "The alignment ancors: " << endl;
 	//for(int i=0;i<anchor_list.size();i++){
	    //  cout << "("<<anchor_list[i].i << ","<< anchor_list[i].j << ")" << endl;
	//}
	anchor_list.clear();
}

// The function reports a history for a given interval, and the root either left or right
void TestingDupHist::show_dup_scenario(TestingDupHist* Sptr,vector <string> types_vector_ptr, int lb, int rb, int root){
    // case 1 where the optimal score is the concatenation of two intervals [i..k][k..j]
    // case 2 where the optimal score is the concatenation of two intervals [i..k][k+1..j] and d(i,j)

    int i=lb;
    int j=rb;
    int k=0;
    int* seqS; 
    seqS=Sptr->sequence1; 
    //int seqSlength=Sptr->seqlength;
    if(i==j)return;
    //cout << "Duplication scenario:"<< i << " " << j << " " << root <<"indel " << Sptr->indel_hist << "\n";
    //return;
    
    if(root==lb){ // leftmost element
	if((j-i)<=1){
	    if((dist[seqS[i]][seqS[j]]+Sptr->fixed_dup_cost)==Sptr->histmatrix_out[i][j][0]){
		cout << "#                   "<< i << " -> " << j<< ", d(" <<i << "," << j<< ")=" <<dist[seqS[i]][seqS[j]]<< endl;
	    }
	    else if(Sptr->indel_hist==Sptr->histmatrix_out[i][j][0]){
		cout << "#                   indel"  << "  " << j <<endl;
	    }
	    else{
		cout << "## there is a problem, no clear operation: "<< Sptr->histmatrix_out[i][j][0] << " Indel " << indel_hist<<"\n";
	    }
	    return;
	}
	for (k=i+1;k<j;k++){ // case 2
	    if(Sptr->histmatrix_out[i][j][0]==(Sptr->histmatrix_out[i][k][0]+Sptr->histmatrix_out[k][j][0])){
		//cout << "Case1: (" << i << " , " << j<< ") from (" << i<< " , " << k << ") + (" << k<< " , " << j<<")"<< endl;
		show_dup_scenario(Sptr,types_vector_ptr,i,k,i);
		show_dup_scenario(Sptr,types_vector_ptr,k,j,k);
		return;
	    }	 
	}
	//cout << "Oops1: \n";
	for (k=i+1;k<j;k++){ // case 2
	    //cout << "!" <<Sptr->histmatrix_out[i][j][0] <<" "<< Sptr->histmatrix_out[i][k][0] << " " 
	    //	 << Sptr->histmatrix_out[j][k+1][0] << " "<<dist[seqS[i]][seqS[j]] << " "<<Sptr->fixed_dup_cost<< endl;
	    if(Sptr->histmatrix_out[i][j][0]==(Sptr->histmatrix_out[i][k][0]+Sptr->histmatrix_out[j][k+1][0]+dist[seqS[i]][seqS[j]]+Sptr->fixed_dup_cost)){
		//cout << "Case2: (" << i << " , " << j<< ") from (" << i<< 
		//" , " << k << ") + (" << k+1<< " , " << j<<")+d(" <<i << " , " << j<< ")" << endl;
		cout << "#                   " << i << " -> " << j<< ", d(" <<i << "," << j<< ")=" <<dist[seqS[i]][seqS[j]]<< endl;//Sptr->histmatrix_out[i][j][0]<<endl;
		show_dup_scenario(Sptr,types_vector_ptr,i,k,i);
		show_dup_scenario(Sptr,types_vector_ptr,k+1,j,j);
		return;
	    }	 
	}
	
	// check for indels
	for (k=i+1;k<j;k++){ // case 2
	    //cout << "!" <<Sptr->histmatrix_out[i][j][0] <<" "<< Sptr->histmatrix_out[i][k][0] << " " 
	    //	 << Sptr->histmatrix_out[j][k+1][0] << " "<<dist[seqS[i]][seqS[j]] << " "<<Sptr->fixed_dup_cost<< endl;
	    if(Sptr->histmatrix_out[i][j][0]==(Sptr->histmatrix_out[i][k][0]+Sptr->histmatrix_out[j][k+1][0]+Sptr->indel_hist)){
		//cout << "Case2: (" << i << " , " << j<< ") from (" << i<< 
		//" , " << k << ") + (" << k+1<< " , " << j<<")+d(" <<i << " , " << j<< ")" << endl;
		cout << "#                   indel"  << "  " << j <<endl;
		show_dup_scenario(Sptr,types_vector_ptr,i,k,i);
		show_dup_scenario(Sptr,types_vector_ptr,k+1,j,j);
		return;
	    }	 
	}
	
	//cout << "Oops2: \n";
	
	
    }
    else{ // rightmost element
	if((j-i)<=1){
	    //cout << "Case0: (" << i << " , " << j<< ") from d(" <<j << "," << i<< ")" <<endl;
	    if((dist[seqS[j]][seqS[i]]+Sptr->fixed_dup_cost)==Sptr->histmatrix_out[j][i][0]){
		cout << "#                   " << j << " -> " << i<< ", d(" <<j << "," << i<< ")=" <<dist[seqS[j]][seqS[i]]<< endl;//Sptr->histmatrix_out[j][i][0]<<endl;
	    }
	    else if(Sptr->indel_hist==Sptr->histmatrix_out[j][i][0]){
		cout << "#                   indel"  << "  " << i <<endl;
	    }
	    else{
		cout << "## there is a problem, no clear operation: "<< Sptr->histmatrix_out[i][j][0] << " Indel " << Sptr->indel_hist<< " dist "<<dist[seqS[j]][seqS[i]] <<"\n";
	    }
	}
	for (k=i+1;k<j;k++){ // case 2
	    if(Sptr->histmatrix_out[j][i][0]==(Sptr->histmatrix_out[j][k][0]+Sptr->histmatrix_out[k][i][0])){
		//cout << "Case1: (" << i << " , " << j<< ") from (" << j<< "," << k << ") + (" << k<< " , " << i<<")"<< endl;
		show_dup_scenario(Sptr,types_vector_ptr,i,k,k);
		show_dup_scenario(Sptr,types_vector_ptr,k,j,j);
		return;
	    }	 
	}
	for (k=i+1;k<j;k++){ // case 2
	    if(Sptr->histmatrix_out[j][i][0]==(Sptr->histmatrix_out[j][k+1][0]+Sptr->histmatrix_out[i][k][0]+dist[seqS[j]][seqS[i]]+Sptr->fixed_dup_cost)){
		//cout << "Case2: (" << i << " , " << j<< ") from (" << j<< 
		//" , " << k << ") + (" << k+1<< " , " << i<<")+d(" <<j << " , " << i<< ")" << endl;
		cout << "#                   " << j << " -> " << i<< ", d(" <<j << "," << i<< ")=" <<dist[seqS[j]][seqS[i]]<< endl;//Sptr->histmatrix_out[j][i][0]<<endl;
		show_dup_scenario(Sptr,types_vector_ptr,k+1,j,j);
		show_dup_scenario(Sptr,types_vector_ptr,i,k,i);
		return;
	    }	 
	}
	for (k=i+1;k<j;k++){ // case 2
	    if(Sptr->histmatrix_out[j][i][0]==(Sptr->histmatrix_out[j][k+1][0]+Sptr->histmatrix_out[i][k][0]+Sptr->indel_hist)){
		//cout << "Case2: (" << i << " , " << j<< ") from (" << j<< 
		//" , " << k << ") + (" << k+1<< " , " << i<<")+d(" <<j << " , " << i<< ")" << endl;
		cout << "#                   indel" << " " << i <<endl;
		show_dup_scenario(Sptr,types_vector_ptr,k+1,j,j);
		show_dup_scenario(Sptr,types_vector_ptr,i,k,i);
		return;
	    }	 
	}
    }
    
}
////////////////////////////////////////////////////////////////////////////////
/// Modified Benson-Dong based method allowing INSERTIONS
//////////////////////////////////////////////////////////////////////////7

int TestingDupHist::SymmBensonDongWithInsertion(){

    int i,j,k,sl,sx;
	int tempval1=(int)0;
	int tempval2=(int)0;
	int tempval3=(int)BigValue;
	int tempval31=(int)0;
	int tempval32=(int)0;
	int tempval33=(int)0;
//	char originchar;


	//cout << "MODIFIED Symmmmmmmmmmmmmmmmmmmm Dong: " <<indel_hist << endl;
	FILE* fptr;
//	fptr=fopen("TempXFGresult.txt","w");

	for (i=1;i<=seqlength;i++){ // i is the interval length
		for(j=0;j+i-1<seqlength;j++){
			// the interval boundaries are [j..j+i-1]
			if(i==1){ //initialization			    
			    for(k=j;k<j+i;k++){
				histmatrix_out[k][k][0]=(int)0;//dist[sequence1[k]][sx];				
			    }				
					
			}
			else if(i==2){
				sx=0;
				histmatrix_out[j][j+i-1][0]=min(dist[sequence1[j+i-1]][sequence1[j]]+fixed_dup_cost,indel_hist);
				histmatrix_out[j+i-1][j][0]=min(dist[sequence1[j]][sequence1[j+i-1]]+fixed_dup_cost,indel_hist);
//				fprintf(fptr,"[%d..%d] score: %2f\n",j,j+i-1,histmatrix_out[j][j+i-1][sx]);
			}
			else{
				
				// debug
				//if((j==0)&&(i==seqlength)){
				//	int kk=0;
				//}
				//sx=0;

				// our second correction for unsymmetric case
				// we have to consider to cases the interval is originated from
				// the leftmost or the right most.
				// case 1
			
				for(k=j+1;k<j+i-1;k++){
					// originated from the left side
					tempval1=histmatrix_out[j][k][0]+histmatrix_out[k][j+i-1][0];
					histmatrix_out[j][j+i-1][0]=min(tempval1,histmatrix_out[j][j+i-1][0]);
					// originated from the right side
					tempval1=histmatrix_out[k][j][0]+histmatrix_out[j+i-1][k][0];
					histmatrix_out[j+i-1][j][0]=min(tempval1,histmatrix_out[j+i-1][j][0]);

				}


				for(sl=j;sl<j+i-1;sl++){
					// originated from the left most
					tempval2=histmatrix_out[j][sl][0]+histmatrix_out[j+i-1][sl+1][0]+
						min(dist[sequence1[j+i-1]][sequence1[j]]+fixed_dup_cost,indel_hist);
					histmatrix_out[j][j+i-1][0]=min(tempval2,histmatrix_out[j][j+i-1][0]);
					// originated from rightmost
					tempval3=histmatrix_out[j][sl][0]+histmatrix_out[j+i-1][sl+1][0]+
						min(dist[sequence1[j]][sequence1[j+i-1]]+fixed_dup_cost,indel_hist);
					histmatrix_out[j+i-1][j][0]=min(tempval3,histmatrix_out[j+i-1][j][0]);

				}

				// end of second correction
				
//				fprintf(fptr,"[%d..%d] score left: %d, score right %d,\n",j,j+i-1,histmatrix_out[j][j+i-1][0],histmatrix_out[j+i-1][j][0]);
		
			}// else 
			//	fprintf(fptr,"[%d..%d] score left: %d, score right %d,\n",j,j+i-1,histmatrix_out[j][j+i-1][0],histmatrix_out[j+i-1][j][0]);
		}
	}

	// Final Full Tree
//	fprintf(fptr,"An optimal final trees \n");
//	fprintf(fptr,"[%d..%d] score: %.2f \n",0,(seqlength-1),min(histmatrix_out[0][seqlength-1][0],histmatrix_out[seqlength-1][0][0]));
	
//	fclose(fptr);
	bensonleft=histmatrix_out[0][seqlength-1][0];
	bensonright=histmatrix_out[seqlength-1][0][0];

	// added to consider length of duplication
	/*
	for(int x=0;x<seqlength;x++){
		for(int y=0;y<seqlength;y++){		
		    //histmatrix_out[x][y][0]=histmatrix_out[x][y][0]+(abs((x-y)*fixed_dup_cost));
		
		}
	}
	*/
	//// finish addition	
	return(min(histmatrix_out[0][seqlength-1][0],histmatrix_out[seqlength-1][0][0]));
}

////////////////////////////////////////////////////////////////////////////////
/// Modified Benson-Dong with RLE and INSERTION
/// The function works as follows:
// construct RLE of sequence 1, add inverse pointers
// construct history for RLE sequence,
// file back the original matrix 
//////////////////////////////////////////////////////////////////////////7

int TestingDupHist::SymmBensonDongWithInsertionRLE(){

    int i,j,k,sl,sx;
    int tempval1=(int)0;
    int tempval2=(int)0;
    int tempval3=(int)BigValue;
    int tempval31=(int)0;
    int tempval32=(int)0;
    int tempval33=(int)0;


    

    inverse_ptr=new int [seqlength]; // store for each element its index in the RLE_seq
    seqlengthRLE =0;
    

    histmatrix_outRLE=NULL;

    seqlengthRLE =compute_rle(seqlength);
    

    
//    FILE* fptr;
//    fptr=fopen("bensonresult.txt","w");

    

    for (i=1;i<=seqlengthRLE;i++){ // i is the interval length
	for(j=0;j+i-1<seqlengthRLE;j++){
	    // the interval boundaries are [j..j+i-1]
	    
	    if(i==1){ //initialization			    
		for(k=j;k<j+i;k++){		    
		    histmatrix_outRLE[k][k][0]=(int)0;//dist[sequence1[k]][sx];			    
		}				
		
	    }
	    else if(i==2){		
		histmatrix_outRLE[j][j+i-1][0]=min(dist[RLE_sequence[j+i-1]][RLE_sequence[j]]+fixed_dup_cost,indel_hist);
		histmatrix_outRLE[j+i-1][j][0]=min(dist[RLE_sequence[j]][RLE_sequence[j+i-1]]+fixed_dup_cost,indel_hist);
		//cout << "Dist  "<< i << " " << j << " " <<histmatrix_outRLE[j+i-1][j][0]  
		//<< " " << indel_hist 
		//    << endl;
//				fprintf(fptr,"[%d..%d] score: %2f\n",j,j+i-1,histmatrix_out[j][j+i-1][sx]);
	    }
	    else{	
		
		// our second correction for unsymmetric case
		// we have to consider to cases the interval is originated from
		// the leftmost or the right most.
		// case 1
		
		for(k=j+1;k<j+i-1;k++){
		    // originated from the left side
		    tempval1=histmatrix_outRLE[j][k][0]+histmatrix_outRLE[k][j+i-1][0];
		    histmatrix_outRLE[j][j+i-1][0]=min(tempval1,histmatrix_outRLE[j][j+i-1][0]);
		    // originated from the right side
		    tempval1=histmatrix_outRLE[k][j][0]+histmatrix_outRLE[j+i-1][k][0];
		    histmatrix_outRLE[j+i-1][j][0]=min(tempval1,histmatrix_outRLE[j+i-1][j][0]);
		    
		}
		
		
		for(sl=j;sl<j+i-1;sl++){
		    // originated from the left most
		    tempval2=histmatrix_outRLE[j][sl][0]+histmatrix_outRLE[j+i-1][sl+1][0]+
			min(dist[RLE_sequence[j+i-1]][RLE_sequence[j]]+fixed_dup_cost,indel_hist);
		    histmatrix_outRLE[j][j+i-1][0]=min(tempval2,histmatrix_outRLE[j][j+i-1][0]);
		    // originated from rightmost
		    tempval3=histmatrix_outRLE[j][sl][0]+histmatrix_outRLE[j+i-1][sl+1][0]+
			min(dist[RLE_sequence[j]][RLE_sequence[j+i-1]]+fixed_dup_cost,indel_hist);
		    histmatrix_outRLE[j+i-1][j][0]=min(tempval3,histmatrix_outRLE[j+i-1][j][0]);
		    
		}
		
		// end of second correction
		
		//	fprintf(fptr,"[%d..%d] score left: %d, score right %d,\n",j,j+i-1,histmatrix_outRLE[j][j+i-1][0],histmatrix_outRLE[j+i-1][j][0]);
		
	    }// else 
//	    fprintf(fptr,"[%d..%d] score left: %d, score right %d,\n",j,j+i-1,histmatrix_outRLE[j][j+i-1][0],histmatrix_outRLE[j+i-1][j][0]);
	}
    }


	
    

    // file back the original matrix 

    //   fprintf(fptr,"/**************** COST OF EXPANDED \n");    

    for (i=1;i<=seqlength;i++){ // i is the interval length
	for(j=0;j+i-1<seqlength;j++){
	    // the interval boundaries are [j..j+i-1]
	    if(i==1){ //initialization			    
		for(k=j;k<j+i;k++){		    
		    histmatrix_out[k][k][0]=(int)0;//dist[sequence1[k]][sx];			    
		}				
		
	    }    
	    else{
		histmatrix_out[j][j+i-1][0]=histmatrix_outRLE[inverse_ptr[j]][inverse_ptr[j+i-1]][0]+
		    (abs(i-1)*fixed_dup_cost-(abs(inverse_ptr[j]-inverse_ptr[j+i-1])*fixed_dup_cost));
		histmatrix_out[j+i-1][j][0]=histmatrix_outRLE[inverse_ptr[j+i-1]][inverse_ptr[j]][0]+
		    (abs(i-1)*fixed_dup_cost-(abs(inverse_ptr[j+i-1]-inverse_ptr[j])*fixed_dup_cost));
	    }
//	    fprintf(fptr,"[%d..%d] score left: %d, score right %d,\n",j,j+i-1,histmatrix_out[j][j+i-1][0],histmatrix_out[j+i-1][j][0]);
	}
    }

//    fclose(fptr);
    
    bensonleft=histmatrix_out[0][seqlength-1][0];
    bensonright=histmatrix_out[seqlength-1][0][0];
    
    
    // cleaning up
/*    
    delete[] inverse_ptr;inverse_ptr=NULL;
    
    for(int x=0;x<seqlengthRLE;x++){
	for(int y=0;y<seqlengthRLE;y++){	
	    delete histmatrix_outRLE[x][y];
	}
    }
   
    for(int x=0;x<seqlengthRLE;x++){
	delete[] histmatrix_outRLE[x];	   
    }

    delete[] histmatrix_outRLE; 
    histmatrix_outRLE=NULL;
    
    delete[] RLE_sequence;
  
  ///end cleaning up
  */    
    //// finish addition	
    return(min(histmatrix_out[0][seqlength-1][0],histmatrix_out[seqlength-1][0][0]));
}

int TestingDupHist::compute_rle(int seqlength){

    vector <int> tmp_rle_seq;
    tmp_rle_seq.clear();
    tmp_rle_seq.push_back(sequence1[0]);
    inverse_ptr[0]=0;	    
    int j=0;
    
    for(int i=1;i<seqlength;i++){
	if(sequence1[i]!=sequence1[i-1]){
	    tmp_rle_seq.push_back(sequence1[i]);	    
	    inverse_ptr[i]=tmp_rle_seq.size()-1;
	}
	else{
	    inverse_ptr[i]=inverse_ptr[i-1];
	}	
    }
    RLE_sequence=new int [tmp_rle_seq.size()];
    for(int i=0;i<tmp_rle_seq.size();i++){
	RLE_sequence[i]=tmp_rle_seq[i];		
    }


    

    //// construct histmatrix for rLE sequence
    histmatrix_outRLE=NULL;
    histmatrix_outRLE=new int ** [tmp_rle_seq.size()]; // for x
    if(histmatrix_outRLE==NULL){cerr << "Error: not enough memory"; exit(-1);}
    for(int x=0;x<tmp_rle_seq.size();x++){
	histmatrix_outRLE[x]=new int*[tmp_rle_seq.size()];
	if(histmatrix_outRLE[x]==NULL){cerr << "Error: not enough memory"; exit(-1);}
    }
    for(int x=0;x<tmp_rle_seq.size();x++){
	for(int y=0;y<tmp_rle_seq.size();y++){	
	    histmatrix_outRLE[x][y]=new int[1];
	    if(histmatrix_outRLE[x][y]==NULL){cerr << "Error: not enough memory"; exit(-1);}
	}
    }
  
    for(int x=0;x<tmp_rle_seq.size();x++){
	for(int y=0;y<tmp_rle_seq.size();y++){
	    histmatrix_outRLE[x][y][0]=BigValue;	    
	}	
    }
    //// end construct intialize histRLE
    
    return(tmp_rle_seq.size());
}
//////////////////////////////////////////////////////////////7
////////////////////////////////////////////////////////////////////
/////////////      and N^3 Implementation         //////////////////
int TestingDupHist::AlignmentRLE(TestingDupHist* Sptr,TestingDupHist* Rptr,int only_left_flag){
	
	int i,j,ls,lr,ts,tr;
	int seqSlength,seqRlength;
	int tempval1=BigValue;
	int tempval2=BigValue;

	int* seqS; 
	int* seqR;
	
	//Sptr->histmatrix_out;
	//histmatrix2=(int*)Rptr->histmatrix_out;
	seqS=Sptr->sequence1; 
	seqR=Rptr->sequence1;

	seqSlength=Sptr->seqlength;
	seqRlength=Rptr->seqlength;

	

	// initialize arrays
	if(1){
	    Amatrix=new int*[seqSlength+1];
	    Bmatrix=new int*[seqSlength+1];
	    for(int y=0;y<seqSlength+1;y++){
		Amatrix[y]=new int[seqRlength+1];
		Bmatrix[y]=new int[seqRlength+1];
	    }
	    for (i=0;i<seqSlength+1;i++){
		for(j=0;j<seqRlength+1;j++){
		    Amatrix[i][j]=(int)BigValue;
		    Bmatrix[i][j]=(int)BigValue;
		}
	    }
	}
	/// end intialize create array
	// construct anchor arrays for S and R
	// each entry in the array conatain a pointer to the position
	// where the anchor exists
	int * anchorS=new int [seqSlength];
	int * anchorR=new int [seqRlength];
	int tmpanchorpos=0;

	anchorS[0]=-1;	
	for(int i=1;i<seqSlength;i++){	    
	    if(seqS[i]==seqS[i-1]){
		anchorS[i]=anchorS[i-1];	
	    }
	    else{
		anchorS[i]=i-1;
	    }
	}
	anchorR[0]=-1;
	for(int i=1;i<seqRlength;i++){	    
	    if(seqR[i]==seqR[i-1]){
		anchorR[i]=anchorR[i-1];	
	    }
	    else{
		anchorR[i]=i-1;
	    }
	}
	

	Amatrix[0][0]=min(dist[seqS[0]][seqR[0]],dist[seqR[0]][seqS[0]]);
	Bmatrix[0][0]=min(dist[seqS[0]][seqR[0]],dist[seqR[0]][seqS[0]]);
	
	for (i=0;i<seqSlength;i++){
		for(j=0;j<seqRlength;j++){
						
			if((i==0)&&(j==1)){
				int kk=0;
			}			
			// The match case
			if((i>0)&&(j>0))
			    Amatrix[i][j]=Amatrix[i-1][j-1]+min(dist[seqS[i]][seqR[j]],dist[seqR[j]][seqS[i]]);

			//cout << "****** Matirx"<< i <<j << " "<< Amatrix[i][j] << endl;
			// leftside originated duplication
			// left duplication in S case i inside a run
			// left duplication in S case i anchor (last element of a run) a run
			if((i>0)&&(seqS[i]!=seqS[i-1])){ // to check if i is withing a run
			    if(i>0){			    
				ls=i-1;
				Amatrix[i][j]=min(Amatrix[ls][j]+Sptr->histmatrix_out[ls][i][0],Amatrix[i][j]);		       			    
			    }
			    for(ls=anchorS[i];ls>=0;ls=anchorS[ls]){			    
				Amatrix[i][j]=min(Amatrix[ls][j]+Sptr->histmatrix_out[ls][i][0],Amatrix[i][j]);			    	   		
			    }
			}
			else if (i>0){
			    ls=i-1;
			    Amatrix[i][j]=min(Amatrix[ls][j]+Sptr->histmatrix_out[ls][i][0],Amatrix[i][j]);	
			}

			// left duplication in R case j inside a run
			// left duplication in R case j anchor (last element of a run) a run
			// check if within a run
			if((j>0)&&(seqR[j]!=seqR[j-1])){
			    if(j>0){			    
				lr=j-1;
				Amatrix[i][j]=min(Amatrix[i][lr]+Rptr->histmatrix_out[lr][j][0],Amatrix[i][j]);			    
			    }			
			    for(lr=anchorR[j];lr>=0;lr=anchorR[lr]){
				Amatrix[i][j]=min(Amatrix[i][lr]+Rptr->histmatrix_out[lr][j][0],Amatrix[i][j]);
			    }
			}
			else if(j>0){
			    lr=j-1;
			    Amatrix[i][j]=min(Amatrix[i][lr]+Rptr->histmatrix_out[lr][j][0],Amatrix[i][j]);	
			}

			if(only_left_flag==1)continue; // ignore right duplications

			// Filling Bmatrix w.r.t. j		
			if((j>0)&&(seqR[j]!=seqR[j-1]))
			{
			    for(tr=j;tr>=1;tr=anchorR[tr]){			    
				// the following line removed to show it is enough
				// to iterate over leftmost unit of each block only
				//Bmatrix[i][j]=min(Amatrix[i][tr-1]+Rptr->histmatrix_out[j][tr][0],Bmatrix[i][j]);
				if(anchorR[tr]>=0){
				    if(anchorR[anchorR[tr]]>=0){
					Bmatrix[i][j]=min(Amatrix[i][ anchorR[anchorR[tr]] ]+
							  Rptr->histmatrix_out[j][anchorR[anchorR[tr]]+1][0],Bmatrix[i][j]);
				    }
				}										
			    }
			}
			if(j>1) 
			{			    			    
			    Bmatrix[i][j]=min(Amatrix[i][j-2]+Rptr->histmatrix_out[j][j-1][0],Bmatrix[i][j]);
			}
			if(j>0){
			    Bmatrix[i][j]=min(Amatrix[i][j-1]+Rptr->histmatrix_out[j][j][0],Bmatrix[i][j]);
			}
			
			// rightside originated duplication		
			if(i==0){
				Amatrix[i][j]=min(Rptr->histmatrix_out[j][0][0]+dist[seqS[i]][seqR[j]],Amatrix[i][j]);
				Amatrix[i][j]=min(Rptr->histmatrix_out[j][0][0]+dist[seqR[j]][seqS[i]],Amatrix[i][j]);
				
			}
			else if(j==0){
				Amatrix[i][j]=min(Sptr->histmatrix_out[i][0][0]+dist[seqS[i]][seqR[j]],Amatrix[i][j]);
				Amatrix[i][j]=min(Sptr->histmatrix_out[i][0][0]+dist[seqR[j]][seqS[i]],Amatrix[i][j]);

			}			
			else{
			    int within_flag=0;
			    if((i>1)&&((seqS[i]!=seqS[i-1])))//||(seqS[i-2]!=seqS[i-1])))
			    {
				for(ts=i;ts>=0;ts=anchorS[ts]){
				    if(ts>0){
					// removed to show that it is enought to iterate
					// over the leftmost unit only
					//Amatrix[i][j]=min(Bmatrix[ts-1][j]+Sptr->histmatrix_out[i][ts][0]+dist[seqS[i]][seqR[j]],Amatrix[i][j]);
					//Amatrix[i][j]=min(Bmatrix[ts-1][j]+Sptr->histmatrix_out[i][ts][0]+dist[seqR[j]][seqS[i]],Amatrix[i][j]); 	
				    }
				    if(anchorS[ts]>=0){
					if(anchorS[anchorS[ts]]>=0){
					    Amatrix[i][j]=min(Bmatrix[anchorS[anchorS[ts]]][j]+Sptr->histmatrix_out[i][anchorS[anchorS[ts]]+1][0]+dist[seqS[i]][seqR[j]],Amatrix[i][j]);
					    Amatrix[i][j]=min(Bmatrix[anchorS[anchorS[ts]]][j]+Sptr->histmatrix_out[i][anchorS[anchorS[ts]]+1][0]+dist[seqR[j]][seqS[i]],Amatrix[i][j]);
					    
					}
				    }
				}
				within_flag=1;
			    }
			    if(i>1)
			    {
				Amatrix[i][j]=min(Bmatrix[i-2][j]+Sptr->histmatrix_out[i][i-1][0]+dist[seqS[i]][seqR[j]],Amatrix[i][j]);
				Amatrix[i][j]=min(Bmatrix[i-2][j]+Sptr->histmatrix_out[i][i-1][0]+dist[seqR[j]][seqS[i]],Amatrix[i][j]); 
				
				Amatrix[i][j]=min(Bmatrix[i-1][j]+Sptr->histmatrix_out[i][i][0]+dist[seqS[i]][seqR[j]],Amatrix[i][j]);
				Amatrix[i][j]=min(Bmatrix[i-1][j]+Sptr->histmatrix_out[i][i][0]+dist[seqR[j]][seqS[i]],Amatrix[i][j]);
				
			    }
			    if(i>0){				
				Amatrix[i][j]=min(Bmatrix[i-1][j]+Sptr->histmatrix_out[i][i][0]+dist[seqS[i]][seqR[j]],Amatrix[i][j]);
				Amatrix[i][j]=min(Bmatrix[i-1][j]+Sptr->histmatrix_out[i][i][0]+dist[seqR[j]][seqS[i]],Amatrix[i][j]);
			    }
			    
			    if(within_flag==1){
				within_flag=0;
				continue;
			    }
			    // to consider if within a run and boundary cases

			    
			    Amatrix[i][j]=min(Amatrix[anchorS[i]][anchorR[j]]+Rptr->histmatrix_out[j][anchorR[j]+1][0]
					      +Sptr->histmatrix_out[i][anchorS[i]+1][0]+dist[seqS[i]][seqR[j]],Amatrix[i][j]);
			    Amatrix[i][j]=min(Amatrix[anchorS[i]][anchorR[j]]+Rptr->histmatrix_out[j][anchorR[j]+1][0]
					      +Sptr->histmatrix_out[i][anchorS[i]+1][0]+dist[seqR[j]][seqS[i]],Amatrix[i][j]); 

			    Amatrix[i][j]=min(Amatrix[anchorS[i]+1][anchorR[j]+1]+Rptr->histmatrix_out[j][anchorR[j]+2][0]
			    		      +Sptr->histmatrix_out[i][anchorS[i]+2][0]+dist[seqS[i]][seqR[j]],Amatrix[i][j]);
			    Amatrix[i][j]=min(Amatrix[anchorS[i]+1][anchorR[j]+1]+Rptr->histmatrix_out[j][anchorR[j]+2][0]
			    		      +Sptr->histmatrix_out[i][anchorS[i]+2][0]+dist[seqR[j]][seqS[i]],Amatrix[i][j]); 
			    if((anchorS[i]!=0)&&(anchorR[j]!=0)){
				Amatrix[i][j]=min(Amatrix[anchorS[i]-1][anchorR[j]-1]+Rptr->histmatrix_out[j][anchorR[j]][0]
						  +Sptr->histmatrix_out[i][anchorS[i]][0]+dist[seqS[i]][seqR[j]],Amatrix[i][j]);
				Amatrix[i][j]=min(Amatrix[anchorS[i]-1][anchorR[j]-1]+Rptr->histmatrix_out[j][anchorR[j]][0]
						  +Sptr->histmatrix_out[i][anchorS[i]][0]+dist[seqR[j]][seqS[i]],Amatrix[i][j]); 
			    }
			
			    Amatrix[i][j]=min(Rptr->histmatrix_out[j][0][0]+Sptr->histmatrix_out[i][0][0]+dist[seqS[i]][seqR[j]],Amatrix[i][j]);
			    Amatrix[i][j]=min(Rptr->histmatrix_out[j][0][0]+Sptr->histmatrix_out[i][0][0]+dist[seqR[j]][seqS[i]],Amatrix[i][j]);
		
			    if(i>1){
				
				Amatrix[i][j]=min(Amatrix[i-2][j-1]+Sptr->histmatrix_out[i][i-1][0]+dist[seqS[i]][seqR[j]],Amatrix[i][j]);
				Amatrix[i][j]=min(Amatrix[i-2][j-1]+Sptr->histmatrix_out[i][i-1][0]+dist[seqR[j]][seqS[i]],Amatrix[i][j]); 			       
			    }
			    if(j>1){
				Amatrix[i][j]=min(Amatrix[i-1][j-2]+Rptr->histmatrix_out[j][j-1][0]+dist[seqS[i]][seqR[j]],Amatrix[i][j]);
				Amatrix[i][j]=min(Amatrix[i-1][j-2]+Rptr->histmatrix_out[j][j-1][0]+dist[seqR[j]][seqS[i]],Amatrix[i][j]); 
				
			    }		       			
			    if((i>1)&&(j>1)){
				Amatrix[i][j]=min(Amatrix[i-2][j-2]+Rptr->histmatrix_out[j][j-1][0]+Sptr->histmatrix_out[i][i-1][0]
						  +dist[seqS[i]][seqR[j]],Amatrix[i][j]);
				Amatrix[i][j]=min(Amatrix[i-2][j-2]+Rptr->histmatrix_out[j][j-1][0]+Sptr->histmatrix_out[i][i-1][0]
						  +dist[seqR[j]][seqS[i]],Amatrix[i][j]); 
			    }
			  
			    //end considering if within a run and boundary cases
			}
			

		}
	}
	
	if(show_alignment_flag==1){
	    show_alignment(Sptr,Rptr,types_vector_ptr);
	}
	
	if(show_breakeven_flag==1){
	    print_breakeven(Sptr,Rptr,types_vector_ptr);
	}
	
	int optimal_Score=Amatrix[seqSlength-1][seqRlength-1];
	if(Amatrix!=NULL){
	    for(int y=0;y<seqSlength+1;y++){
		delete [] Amatrix[y];		
	    }
	    delete [] Amatrix;
	}
	
	if(Bmatrix!=NULL){
	    for(int y=0;y<seqSlength+1;y++){
		delete [] Bmatrix[y];	    
	    }
	    delete [] Bmatrix;
	}

	delete [] anchorS;
	delete [] anchorR;
	
	return(optimal_Score);
	

}


////////////////////////////////////////////////////

//Print breakeven point

void TestingDupHist::print_breakeven(TestingDupHist* Sptr,TestingDupHist* Rptr,vector <string> types_vector_ptr ){

    	int i,j,ls,lr,ts,tr;
	int seqSlength,seqRlength;
	int tempval1=BigValue;
	int tempval2=BigValue;

	int* seqS; 
	int* seqR;

	seqS=Sptr->sequence1; 
	seqR=Rptr->sequence1;

	seqSlength=Sptr->seqlength;
	seqRlength=Rptr->seqlength;
	string aligned_seq_S;
	string aligned_seq_R;
	anchor_pair anchor_pair_obj;
	
	vector <anchor_pair> anchor_list;
	int anchor_flag=0;
	i=seqSlength-1;
	j=seqRlength-1;
	int found_flag=0;
	
	int optimalscore=Amatrix[i][j];
	float breakevenscore=(float)Amatrix[i][j]/(float)2;
	
	int reached_flag=0;
	float reach_i=seqSlength-1; 
	float reach_j=seqRlength-1;
	
	while((i!=-1)||(j!=-1)){	
	    
	    found_flag=0;
	    
	    // The match case
	    
	    if(((i==0)&&(j==0))&&(Amatrix[i][j]==min(dist[seqS[i]][seqR[j]],dist[seqR[j]][seqS[i]]))){
	      //cout << "Match             : " <<  "("<< i <<", "<<j<<")"<< " score: "<<Amatrix[i][j] << endl;
	      	if((Amatrix[i][j]<=breakevenscore)&&(reached_flag==0)){
		  reached_flag=1;
		  reach_i=i;
		  reach_j=j;
		}
		
		i--;
		j--;	
		anchor_flag=1;
		found_flag=1;
		continue;
	    }
	    
	    if((i>0)&&(j>0)){
		if(Amatrix[i][j]==(Amatrix[i-1][j-1]+min(dist[seqS[i]][seqR[j]],dist[seqR[j]][seqS[i]]))){		
		    anchor_flag=1;
		    
		    //cout << "HIIIIIIIIIII Match Entering loop\n";
		    if((Amatrix[i-1][j-1]<=breakevenscore)&&(reached_flag==0)){
		      reached_flag=1;
		      reach_i=i;
		      reach_j=j;
		    }
		    //int countij=0;
	
			/*
			  while((Amatrix[i][j]==Amatrix[i+1][j+1])&&(i<seqSlength)&&(j<seqRlength)){
			  countij++;
			  i++; j++;
			  }
			*/
		    int countij=0;
		    int k1=(int)reach_i;
		    int k2=(int)reach_j;  // stopped here mohamed
		    while((k1>1)&&(k2>1)&&(Amatrix[k1][k2]==Amatrix[k1-1][k2-1])&&(Amatrix[k1][k2]==breakevenscore)){
		      countij++;
		      k1--; k2--;
		    }
		    if(reach_i>0.5)
		      reach_i=reach_i-max(0,(countij/(float)2));
		    if(reach_j>0.5)
		      reach_j=reach_j-max(0,(countij/(float)2));
		    //cout << "HIIIIIIIIIII Match Done Entering loop\n";
		    i--;
		    j--;				    
		    anchor_flag=1;
		    found_flag=1;
		    continue;
		}
	    }
	    
	    // The left duplication cases
	    for (ls=0;ls<i;ls++){
		if(ls<i){
		    if(Amatrix[i][j]==(Amatrix[ls][j]+Sptr->histmatrix_out[ls][i][0])){			
			anchor_flag=1;
			//cout << "Left dup. in S    : " << "[" << ls << ".."<< i   << "]" << " score: "<<Amatrix[i][j] << endl;
			//show_dup_scenario(Sptr,types_vector_ptr, ls, i,ls);
			if((1)&&(Amatrix[ls][j]<breakevenscore)&&(reached_flag==0)){
			  reached_flag=1;
			  // consider fixed_dup_cost
			  // use score diff, yields better estimate
			  
			  //reach_i=min(seqSlength,(ls+1+(Sptr->histmatrix_out[ls][i][0]/(2*Sptr->fixed_dup_cost))));         
			  reach_i=min(seqSlength,(i-((Amatrix[i][j]-breakevenscore)/(Sptr->fixed_dup_cost))));      
			  reach_j=j+1;
			  
			}
			/*
			if((0)&&(Amatrix[ls][j]<breakevenscore)&&(reached_flag==0)){
			  reached_flag=1;
			  // consider fixed_dup_cost
			  reach_i=ls+(breakevenscore-Amatrix[ls][j])/Sptr->fixed_dup_cost;			  
			  reach_j=j;
			}
			*/
		
			i=ls;
			found_flag=1;
			break;
		    }
		}
	    }
	    
	    if(anchor_flag==1){
		anchor_flag=0;
		continue;
	    }
	    
	    for(lr=0;lr<j;lr++){
		if(lr<j){
		    if(Amatrix[i][j]==(Amatrix[i][lr]+Rptr->histmatrix_out[lr][j][0])){
			anchor_flag=1;
			//cout << "Left dup. in R    : "  << "[" << lr << ".."<<j << "] score: "<<Amatrix[i][j] << endl;
			//show_dup_scenario(Rptr,types_vector_ptr, lr, j,lr);			
			  
			if((1)&&(Amatrix[i][lr]<breakevenscore)&&(reached_flag==0)){
			  reached_flag=1;
			  // consider fixed_dup_cost
			  //reach_j=min(seqRlength,(lr+1+(Rptr->histmatrix_out[j][lr+1][0]/(2*Sptr->fixed_dup_cost))));      
			  // use score diff, yields better estimate
			  reach_j=min(seqRlength,(j-((Amatrix[i][j]-breakevenscore)/(Sptr->fixed_dup_cost))));      
			  reach_i=i+1;
			  
			}
			  /*

			if((0)&&(Amatrix[i][lr]<breakevenscore)&&(reached_flag==0)){
			  reached_flag=1;
			  // consider fixed_dup_cost
			  reach_j=lr+(breakevenscore-Amatrix[i][lr])/Sptr->fixed_dup_cost;			  
			  reach_i=i;
			}
			*/
			found_flag=1;
			j=lr;			    
			break;
		    }
		}
	    }
	    if(anchor_flag==1){
		anchor_flag=0;
		continue;
	    }

	    ls=0;
	    lr=0;
	    if(Amatrix[i][j]==(Sptr->histmatrix_out[i][0][0]+Rptr->histmatrix_out[j][0][0]
			       +min(dist[seqS[i]][seqR[j]],dist[seqR[j]][seqS[i]]))){	
	      //cout << "Match             : " << "("<< i <<", "<<j<<")" << " score: "<< Amatrix[i][j] << endl;
		if((ls+1)!=i){
		  //cout << "Right dup. in S   : "  << "["<<ls+1 << ".." << i      << "]"<< " score: "<< Amatrix[i][j] << endl;
		  //show_dup_scenario(Sptr,types_vector_ptr, ls+1, i,i);	      
		}
		if((lr+1)!=j){
		  //cout << "Right dup. in R   : "  << "[" <<lr+1 << ".." << j     << "]"<< " score: "<< Amatrix[i][j] << endl;
		  //show_dup_scenario(Rptr,types_vector_ptr, lr+1, j,j);
		}
		
	
		//cout << "HIIIIIIIIIII RD Entering loop\n";
		if((1)&&(Amatrix[ls][lr]<breakevenscore)&&(reached_flag==0)){
		  reached_flag=1;
		  // consider fixed_dup_cost
		  // consider fixed_dup_cost
		  int scorediff=(Amatrix[i][j]-breakevenscore);
		  int tmptt=((i-ls)+(j-lr));
		  int i_share, j_share;
		  if(tmptt==0) {i_share=0.5; j_share=0.5;} 
		  else {
		    i_share=(i-ls)/tmptt;
		    j_share=(j-lr)/tmptt;
		  }
		  reach_i=min(seqSlength,(i-((scorediff*i_share)/Sptr->fixed_dup_cost)));
		  reach_j=min(seqRlength,(j-((scorediff*j_share)/Sptr->fixed_dup_cost)));
		  

		  /*
		  reach_j=min(seqRlength,(lr+1+
					  (Rptr->histmatrix_out[j][lr+1][0]
					   /(2*Sptr->fixed_dup_cost))));      
		  reach_i=min(seqSlength,(ls+1+
					  (Sptr->histmatrix_out[i][ls+1][0]
					  /(2*Sptr->fixed_dup_cost))));
		  */
		}

		anchor_flag=1;
		found_flag=1;
		
		
		break;			   
	    }
	    
	    for (ls=0;ls<i;ls++){	
		for(lr=0;lr<j;lr++){	    					    	       		  	
		    if((ls<i)&&(lr<j)){
			if(Amatrix[i][j]==(Amatrix[ls][lr]+Sptr->histmatrix_out[i][ls+1][0]+Rptr->histmatrix_out[j][lr+1][0])
			   +min(dist[seqS[i]][seqR[j]],dist[seqR[j]][seqS[i]])){

			  //    cout << "Match             : " << "("<< i <<", "<<j<<")" << " score: "<<Amatrix[i][j] << endl;
			    if((ls+1)!=i){
			      //cout << "Right dup. in S   : "  << "["<< ls+1 << ".." << i <<"]"     << " score: "<< Amatrix[i][j] << endl;
			      //show_dup_scenario(Sptr,types_vector_ptr, ls+1, i,i);
			    }
			    if((lr+1)!=j){
			      //cout << "Right dup. in R   : "  << "["<< lr+1 << ".." << j <<"]"   << " score: "<< Amatrix[i][j] << endl;
			      //show_dup_scenario(Rptr,types_vector_ptr, lr+1, j,j);
			    }	
		 
		
			    if((1)&&(Amatrix[ls][lr]<breakevenscore)&&(reached_flag==0)){
			      reached_flag=1;
			      // consider fixed_dup_cost
			      int scorediff=(Amatrix[i][j]-breakevenscore);
			      int tmptt=((i-ls)+(j-lr));
			      int i_share, j_share;
			      if(tmptt==0) {i_share=0.5; j_share=0.5;} 
			      else {
				i_share=(i-ls)/tmptt;
				j_share=(j-lr)/tmptt;
			      }
			      reach_i=min(seqSlength,(i+1-((scorediff*i_share)/Sptr->fixed_dup_cost)));
			      reach_j=min(seqRlength,(j+1-((scorediff*j_share)/Sptr->fixed_dup_cost)));
			      /*
			      reach_j=min(seqRlength,(lr+1+
						      (Rptr->histmatrix_out[j][lr+1][0]
						       /(2*Sptr->fixed_dup_cost))));      
			      reach_i=min(seqSlength,(ls+1+
						      (Sptr->histmatrix_out[i][ls+1][0]
						       /(2*Sptr->fixed_dup_cost))));
			      */
			    }


			    found_flag=1;
			    i=ls;
			    j=lr;			    
			    anchor_flag=1;
			    break;
			}
		    }
		} // end inner for loop
		if(anchor_flag==1){
		    anchor_flag=0;
		    break;
		}
	    } // end outter for
	    if(found_flag!=1){
	      //cout << "Not found_flag "<< i << " "<<j << endl;
		break;
	    }
	    
	}// end while loop

	// some corrections to break
	//i=(int)reach_i;
	//j=(int)reach_j;

	// we can subtract one to account for dummy character
	cout << "Breakeven i= " << (float)reach_i/(float)(seqSlength) << " Breakeven j= " 
	     << (float)reach_j/(float)(seqRlength) << " i=" << reach_i << " j= "<< reach_j << endl;
	if(anchor_list.size()>0)
	anchor_list.clear();
}
