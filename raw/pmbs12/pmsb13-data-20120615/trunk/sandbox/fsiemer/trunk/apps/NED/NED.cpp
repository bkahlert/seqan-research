#include <iostream>
#include <fstream>
#include <cstdio>
#include <vector>

#include <seqan/graph_types.h>
#include <seqan/graph_algorithms.h>
#include <seqan/stream.h>
#include <seqan/basic.h>
#include <seqan/index.h>
#include <seqan/file.h>
#include <seqan/align.h>

using namespace seqan;

//***** READ_FASTA *****//
template<typename TString>
TString read_fasta(char* const &filename) {
	std::ifstream fasta(filename, std::ios_base::in | std::ios_base::binary);
    if (!fasta.good())
        return 1;
    RecordReader<std::ifstream, SinglePass<> > reader(fasta);
    // Define variables for storing the sequences and sequence ids.
    CharString id;
    TString seq;
    // Read FASTA file.
    while (!atEnd(reader)) {
        if (readRecord(id, seq, reader, Fasta()) != 0) {
            std::cerr << "ERROR reading FASTA.\n";
            return 1;
        }
    }

	return seq;
}

//***** GET_COST *****//
template<typename TScore, typename TAlp>
TScore get_cost(TAlp const &a, TAlp const &b, TScore const &match_cost, TScore const &missm_cost) {
	return (a == b) ? match_cost : missm_cost;
}

//***** ADD *****//
template<typename TScore, typename TScore2>
TScore ADD(TScore const &a, TScore2 const &b) {
	return (a == MaxValue<TScore>::VALUE || b == MaxValue<TScore>::VALUE) ? MaxValue<TScore>::VALUE : (a + b);
}

//***** MIN/ MAX *****//
template<typename TScore>
TScore MAX(TScore const &a, TScore const &b) {
	return (a > b) ? a : b;
}
template<typename TScore>
TScore MIN(TScore const &a, TScore const &b) {
	return (a < b) ? a : b;
}
template<typename TScore>
TScore MIN(TScore const &a, TScore const &b, TScore const &c) {
	return MIN(a, MIN(b, c));
}

//***** MIN_PATH *****//
//this function is a regular MIN(x,y,z) but it also sets the ALI[] matrix by reference
template<typename TScore>
TScore MIN_PATH(TScore const &a, TScore const &b, TScore const &c, char &ALI) {
	ALI = 0;

	if(MIN(a, c, b) == a) {
		ALI = 1;
		return a;
	}
	else if(MIN(a, c, b) == b) {
		ALI = 2;
		return b;
	}

	return c;
}

//***** ned *****//
template<typename TScore, typename TString>
float ned(TString const &x, TString const &y, TScore const &match, TScore const &missm, TScore const &gap) {
	return ned(x, y, match, missm, gap, false);
}
template<typename TScore, typename TString>
float ned(TString const &x, TString const &y, TScore const &match, TScore const &missm, TScore const &gap, bool const &traceback) {
	if(length(y) <= 0 || length(x) <= 0) {
		std::cerr << "ERROR\tFunction 'ned()' requires 2 sequences.\n";
		return 0;
	} 
	if(match < 0 || missm < 0 || gap < 0) {
		std::cerr << "ERROR\tAll scores/ costs in function 'ned()' need to be positive numbers.\n";
		return 0;
	}
	if(length(y) > length(x)) {
		return ned(y, x, match, missm, gap, traceback);
	}

	typedef unsigned short T_DP;
	T_DP INF = MaxValue<T_DP>::VALUE;

	TScore maxcost = TScore(INF)/(2 * length(x));
	if(match > maxcost || missm > maxcost || gap > maxcost) {
		std::cerr << "ERROR\tAll costs/ scores in function 'ned()' shouldn't exceed " << maxcost << ".\n";
		return 0;
	}

	int xlen = length(x);
	int ylen = length(y);
	int n = xlen+1;
	int m = (xlen+1) * (ylen+1);

	T_DP* D = 0;          //the DP matrix
	int* PATH = 0;        //the traceback path: last element -> previous element -> ... -> start
	char* ALI = 0;        //the alinment for the traceback:   1 : gap in y,   2: gap in x,   0: match/ missmatch
	bool* GAP_x = 0;      //used to print the aligment graph
	bool* GAP_y = 0;      //used to print the aligment graph

	long long size = long(xlen+1) * long(ylen+1) * long(xlen+ylen+3);
	D = new T_DP[size];

	if(traceback) {
		PATH = new int[size];
		ALI  = new char[size];

		PATH[ (0) + (0) *n+ (0) *m] = 0;
		ALI[  (0) + (0) *n+ (0) *m] = 0;
	}

	D[ (0) + (0) *n+ (0) *m] = 0;
	D[ (0) + (0) *n+ (1) *m] = INF;
	
	for(int j=1;j<=ylen;++j) {
		//initialize y
		D[ (0) + (j) *n+ (j-1) *m] = INF;
		D[ (0) + (j) *n+ (j) *m] = ADD( D[ (0) + (j-1) *n+ (j-1) *m], gap );
		D[ (0) + (j) *n+ (j+1) *m] = INF;

		if(traceback) {
			PATH[ (0) + (j) *n+ (j) *m] = ( (0) + (j-1) *n+ (j-1) *m);
			ALI[  (0) + (j) *n+ (j) *m] = 2;
		}
	} // endfor: j

	for(int i=1;i<=xlen;++i) {
		//initialize x
		D[ (i) + (0) *n+ (i-1) *m] = INF;
		D[ (i) + (0) *n+ (i) *m] = ADD( D[ (i-1) + (0) *n+ (i-1) *m], gap );
		D[ (i) + (0) *n+ (i+1) *m] = INF;

		if(traceback) {
			PATH[ (i) + (0) *n+ (i) *m] = ( (i-1) + (0) *n+ (i-1) *m);
			ALI[  (i) + (0) *n+ (i) *m] = 1;
		}

		for(int j=1;j<=ylen;++j) {
			D[ (i) + (j) *n+ (MAX(i,j)-1) *m] = INF;
			for(int k=MAX(i,j);k<=(i+j);++k) {
				if(traceback) {
					// get the best score for the next DP-entry and set the traceback path
					D[ (i) + (j) *n+ (k) *m] = MIN_PATH(
						ADD( D[ (i-1) + (j) *n+ (k-1) *m], gap ), 
						ADD( D[ (i) + (j-1) *n+ (k-1) *m], gap ), 
						ADD( D[ (i-1) + (j-1) *n+ (k-1) *m], get_cost(x[i-1], y[j-1], match, missm) ),
						ALI[  (i) + (j) *n+ (k) *m]
					);
					if(ALI[  (i) + (j) *n+ (k) *m] == 1) { PATH[ (i) + (j) *n+ (k) *m] = ( (i-1) + (j) *n+ (k-1) *m); }
					else if(ALI[  (i) + (j) *n+ (k) *m] == 2) { PATH[ (i) + (j) *n+ (k) *m] = ( (i) + (j-1) *n+ (k-1) *m); }
					else { PATH[ (i) + (j) *n+ (k) *m] = ( (i-1) + (j-1) *n+ (k-1) *m); }
				} // endif: traceback
				else {
					// get the best score for the next DP-entry
					D[ (i) + (j) *n+ (k) *m] = MIN(
						ADD( D[ (i-1) + (j) *n+ (k-1) *m], gap ), 
						ADD( D[ (i) + (j-1) *n+ (k-1) *m], gap ), 
						ADD( D[ (i-1) + (j-1) *n+ (k-1) *m], get_cost(x[i-1], y[j-1], match, missm) )
					);
				} // endif: no traceback
			} // endfor: k
			D[ (i) + (j) *n+ (i+j+1) *m] = INF;
		} // endfor: j
	} // endfor: i

	float out = MaxValue<float>::VALUE;  // the output value
	int best_k = -1;                     // length of the aligment, used for traceback

	for(int k=xlen;k<=(xlen+ylen);++k) {
		if((D[ (xlen) + (ylen) *n+ (k) *m] < INF) && out > float(D[ (xlen) + (ylen) *n+ (k) *m]) / float(k) ) {
			out = float(D[ (xlen) + (ylen) *n+ (k) *m]) / float(k);
			best_k = k;
		}
	} // endfor: k

	delete[] D;

	if(traceback) { // print traceback:
		typedef Align<TString, ArrayGaps> TAlign;
		typedef Row<TAlign>::Type TRow;

		TString str1 = x;
		TString str2 = y;
		TAlign align;
		resize(rows(align), 2);
		assignSource(row(align,0),str1);
		assignSource(row(align,1),str2);

		TRow &row1 = row(align,0);
		TRow &row2 = row(align,1);

		//those two arrays help to build the aligment graph in the correct order:
		GAP_x = new bool[best_k];
		GAP_y = new bool[best_k];

		int path_i = (xlen) + (ylen) *n+ (best_k) *m;
		for(int k=(best_k-1);k>=0;--k) { //reverse order: walk the traceback
			if(ALI[path_i] == 1) {
				GAP_x[k] = false;
				GAP_y[k] = true;
			}
			else if(ALI[path_i] == 2) {
				GAP_x[k] = true;
				GAP_y[k] = false;
			}
			else {
				GAP_x[k] = false;
				GAP_y[k] = false;
			}
			path_i = PATH[path_i];
		} // endfor: reverse order
		for(int k=0;k<=(best_k-1);++k) { //correct order: insert gaps
			if(GAP_x[k]) insertGap(row1,k);
			if(GAP_y[k]) insertGap(row2,k);
		} // endfor: correct order

		std::cout << "\nNED alignment:\n" << align;

		if(true) { // show global alignment (default)
			TString str3 = x;
			TString str4 = y;
			TAlign alignG;
			resize(rows(alignG), 2);
			assignSource(row(alignG,0),str3);
			assignSource(row(alignG,1),str4);
			Score<int> scoringScheme( -1 * int(match), -1 * int(missm), -1 * int(gap), -1 * int(gap) );
			int score = globalAlignment(alignG, scoringScheme);

			std::cout << "Global alignment:\n" << alignG;
		}
	} // endif: traceback

	delete[] GAP_x;
	delete[] GAP_y;
	delete[] PATH;
	delete[] ALI;

	return out;
}

int main(int argc, char** const argv) { 
	//Ini:
	typedef String<Dna> TString; 
	double start = sysTime();
	//std::cout << "BUILD#" << 62 << "\n\n";

	//Check parameters:
	bool traceback = true;
	TString seq1, seq2;
	if (argc == 3){
		seq1 = read_fasta< TString >(argv[1]);
		seq2 = read_fasta< TString >(argv[2]);
		traceback = false;
	}
	else if(argc == 4){
		seq1 = read_fasta< TString >(argv[1]);
		seq2 = read_fasta< TString >(argv[2]);
	}
	else {
		seq1 = "ACCCGGT";
		seq2 = "AAACTTG";
		std::cout << "Invalid argument count -> fastaFILE1 fastaFILE2 <optional>SHOWalignment\n--- Initializing defaults ---\nSequences: '" << seq1 << "' and '" << seq2 << "', show alignment\n";
	}

	//Calculation:
	float r = ned(seq1, seq2, 0, 3, 2, traceback);
	std::cout << "NED: " << r << "\n";

	//Time:
	std::cout << "\n\tTime: " << sysTime()-start << "s\n\n";

	return 0;
}