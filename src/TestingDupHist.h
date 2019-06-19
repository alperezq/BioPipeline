// TestingDupHist.h: interface for the TestingDupHist class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_TESTINGDUPHIST)
#define AFX_TESTINGDUPHIST


//#include "DuplicationHistory.h"
#include <string>
#include <vector>

#include <stdio.h>
#include <string.h>
#include <iostream>
#include <stdlib.h>
#include <vector>

#define MaxSeqSize 130
//#define Alphasize 20
#define MaxCodeSize 5
#define MaxAlphaSize 30
#define BigValue 500000

#ifndef min
#define min(a, b)  ((a)<=(b) ? (a) : (b))
#endif

using namespace std;

class TestingDupHist  
{
public:
	TestingDupHist();
	virtual ~TestingDupHist();
	public:
	// Attributes
	//char sequence1[MaxSeqSize];
	int* sequence1;
	//char sequence2[MaxSeqSize];
	int* sequence_rle;

	//int dist[Alphasize][Alphasize];
	int** dist;
///	int histmatrix_out[MaxSeqSize][MaxSeqSize][Alphasize];

	int*** histmatrix_out;

///	int histmatrix_in[MaxSeqSize][MaxSeqSize][Alphasize];
	int*** histmatrix_in;

//	int Amatrix[MaxSeqSize][MaxSeqSize];
	int** Amatrix;

//	int Bmatrix[MaxSeqSize][MaxSeqSize];
	int** Bmatrix;

	int bensonleft,bensonright;
	int oursleft,oursright;
	int seqlength;
	int inputalphabet;
	int Alphasize;
	int reverseflag;

	int fixed_dup_cost;
	int indel_hist;
	int indel_align;
	vector <string> types_vector_ptr;

	int show_alignment_flag;
	// Functions 
	void readsequence(char*);
	void getsequence(int* in_sequence,int);
	
	void allocate_array(int seq_len, int alphadep_flag);
	void FreeMEMORY();
	////////////////////// History Construction Function ///////////
	int SymmBensonDong();
	int SymmBensonDongWithInsertion();
	int constructhistory(char* filename);
	//////////////////////History for RLE sequences
	int SymmBensonDongWithInsertionRLE();
	int compute_rle(int);

	int*** histmatrix_outRLE;
	int* RLE_sequence;
	int* inverse_ptr;
	int seqlengthRLE;

	int AlignmentRLE(TestingDupHist* Sptr,TestingDupHist* Rptr, int only_left_flag);
	///////////////////////////////////////////////////////////////
	int ModAlignment(TestingDupHist* Sptr,TestingDupHist* Rptr);
	void Intialize(int**);
	
	int AlignmentWithoutRightDup(TestingDupHist* Sptr,TestingDupHist* Rptr);
	void intialize_type2distance();
	void intialize_type3distance();
	void intialize_typeBouzekridistance();
	void intialize_typeRtest();

	/////////// display functions //////////////////
	void show_alignment(TestingDupHist* Sptr,TestingDupHist* Rptr,vector <string>);
	void set_show_alignment_flag(int in_show_alignment_flag, vector <string> in_types_vector_ptr)
	    {show_alignment_flag=in_show_alignment_flag;
	    types_vector_ptr=in_types_vector_ptr;};

	void show_dup_scenario(TestingDupHist* Sptr,vector <string> types_vector_ptr, int, int,int);


	char show_breakeven_flag;
	void set_breakeven_flag(int in_breakeven_flag, vector <string> in_types_vector_ptr)
	    {show_breakeven_flag=in_breakeven_flag;
	    types_vector_ptr=in_types_vector_ptr;};
	void print_breakeven(TestingDupHist* Sptr,TestingDupHist* Rptr,vector <string>);
	
	/////////////
	int alphadep_flag;
	void set_alphadep_flag(int in_alphadep_flag){alphadep_flag=in_alphadep_flag;};

	

};
typedef struct { int i; int j;}anchor_pair;

#endif // !defined(AFX_TESTINGDUPHIST_H__D690904B_5E9F_40BB_B374_F112C6B2E04B__INCLUDED_)
