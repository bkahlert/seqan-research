// ==========================================================================
//                               motifHeader.h
// ==========================================================================
// Copyright (c) 2006-2012, Knut Reinert, FU Berlin
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ==========================================================================
// Author: Vishalini Vimalakanthan <vishalini.v@fu-berlin.de>
// ==========================================================================

#ifndef SANDBOX_SHALINI6_APPS_MOTIFSEARCH_MOTIFHEADER_H_
#define SANDBOX_SHALINI6_APPS_MOTIFSEARCH_MOTIFHEADER_H_


#include <iostream>
#include <seqan/find_motif.h>
#include <seqan/find.h>
#include <fstream>
#include <seqan/stream.h>
#include <seqan/align.h>
#include <seqan/refinement.h>
#include <seqan/basic.h>
#include <seqan/graph_types.h>
#include "..\..\..\core\apps\stellar\stellar.h"
#include "..\..\..\core\apps\stellar\stellar_types.h"
#include "..\..\..\core\apps\stellar\stellar_output.h"
#include "..\..\..\core\apps\stellar\stellar_extension.h"
#include "..\..\..\core\include\seqan\graph_msa.h"
#include <seqan/index.h>
#include <vector>


using namespace seqan;

struct MotifNewAlgorithm_;
typedef Tag< MotifNewAlgorithm_> const  MotifNewAlgorithm;



//////////////////////////////////////////////////////////////////////////////
// MotifFinder - MotifFindingAlgorithm Spec
//
// eps:= epsilon
// minL:= minLength  
// x:= xDrop
//
// 
//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename TRng>
class MotifFinder<TValue,  MotifNewAlgorithm, TRng>
{

public:

	typedef unsigned int TSize;
    typedef String<TValue> TString;
    typedef String<TString> TStrings;
	
    Holder<TRng> _rng;
	double epsilon;
	const int minLength;
	TSize xDrop;
	
	MotifFinder():
	_rng(defaultRng<TRng>()),
    epsilon(0),
	minLength(0),
	xDrop(0)

	{
      
    }
	
	MotifFinder(double const & eps_, 
				const int & minL_,
				TSize const & x_):
        _rng(defaultRng<TRng>()),
		epsilon(eps_),
		minLength(minL_),
		xDrop(x_)
		
        
	{
    
	}

	MotifFinder(MotifFinder const & other_):
        _rng(other_.rng),
		epsilon(other_.epsilon),
		minLength(other_.minLength),
		xDrop(other_.xDrop)
		
	{
    
	}
	~MotifFinder()
	{
    
	}

	MotifFinder const &
	operator = (MotifFinder const & other_)
	{
    
		if(this!=&other_)
		{
			this->epsilon = other_.epsilon;
			this->minLength = other_.minLength;
			this->xDrop = other_.xDrop;
			
		}

		return *this;
	}



}; // class MotifFinder<TValue, MotifNewAlgorithm>



//////////////////////////////////////////////////////////////////////////////
// Main function
//////////////////////////////////////////////////////////////////////////////


//findMotif is the main function that is employed to search for
//hidden motifs
//for the initiation it needs a finder, a dataset containing DNA sequences
//and the sequences model Oops()
template<typename TSeqType, typename TStringSet, typename TModel, typename TRng>
void
findMotif(MotifFinder<TSeqType, MotifNewAlgorithm, TRng> & finder, 
		  TStringSet & dataset,
		  TModel seq_model)
{
    
	typedef String<Dna> TSequence;
	typedef CharString TId;
	typedef StringSet<TSequence> TStringSet;
	typedef StellarMatch<TSequence, TId> TMatch;
	StringSet<QueryMatches<TMatch > > matches; //container to store STELLAR matches
	resize(matches, length(dataset));
		

	// ----------------------------------------------------------------------------
	// STEP 1:
	// finding local segment matches using STELLAR
	// ----------------------------------------------------------------------------

	 //computation of local matches for ALL pairs of input sequences
	 for(unsigned r = 0; r<length(dataset)-1; ++r)
	 {

		TSequence &query = dataset[r];
		StringSet<TSequence, Dependent<> > remainingSeqs;

		for(unsigned s = r+1; s<length(dataset);++s)
		{

			appendValue(remainingSeqs, dataset[s]);
		}
			
	
			callStellar(query,remainingSeqs, matches, finder.epsilon, finder.minLength, finder.xDrop);
		

	  }

		// ----------------------------------------------------------------------------
		// STEP 2:
		// matchrefinement of found matches
		// ----------------------------------------------------------------------------
			

	
			typedef typename Graph<Alignment<StringSet<TSequence, Dependent<> > > > TAlignGraph;
			TAlignGraph outgraph(dataset);
			clear(outgraph);

			refinement(matches, dataset, outgraph); 
			

		// ----------------------------------------------------------------------------
		// STEP 3:
		// progressive alignment
        // incorporates:
		// STEP 4:
        // motif extraction
		// STEP 5:
        // determination of a consensus sequence
		// ----------------------------------------------------------------------------

			const int mL = finder.minLength;
			progressiveAlign(outgraph, dataset, mL);

		
			
			
			
		
	 

}
//////////////////////////////////////////////////////////////////////////////




//////////////////////////////////////////////////////////////////////////////
//Subfunctions
//////////////////////////////////////////////////////////////////////////////

template<typename TSequence, typename TStringSet, typename TMatch, typename TSize>
void 
callStellar(TSequence &query, 
			TStringSet &remainingSeqs, 
			StringSet<QueryMatches<TMatch > > &matches,
	        double const & eps_, 
		    const int & minL_,
			TSize const & x_)
{
	
	typedef String<Dna> TSequence;
	typedef StringSet<TSequence, Dependent<> > TStringSet;

	
	//computation of smin (adopted from stellar.cpp ll.320-324)
	int errMinLen = (int) floor(eps_ * minL_);
	int n = (int) ceil((errMinLen + 1) / eps_);
	int errN = (int) floor(eps_ * n);
	unsigned smin = (unsigned) _min(ceil((double)(minL_-errMinLen)/(errMinLen+1)),
	                                    ceil((double)(n-errN)/(errN+1)));
	
	unsigned qGram = smin;  

	//Pattern
	typedef Index<StringSet<TSequence, Dependent<> >, IndexQGram<SimpleShape, OpenAddressing> >  TQGramIndex;
	TQGramIndex qgramIndex(remainingSeqs);
	resize(indexShape(qgramIndex), qGram);
	Pattern<TQGramIndex, Swift<SwiftLocal> > sPattern(qgramIndex);
	indexRequire(qgramIndex, QGramSADir());
	
	//Finder
	typedef Finder<TSequence, Swift<SwiftLocal> > TFinder;
	TFinder sFinder(query);

	
	//calling STELLAR
	stellar(sFinder, sPattern, eps_, minL_, x_, matches, AllLocal());

}


/////////////////////////////////////////////////////////////////////////////////

//STEP 2: refinement of detected local matches
//overlapping matches are refined
template<typename TMatch, typename TStringSet, typename TSequence>
void 
refinement(StringSet<QueryMatches<TMatch > > &matches, 
				   TStringSet &dataset, 
				   Graph<Alignment<StringSet<TSequence, Dependent<> > > > &outgraph)
{

	typedef String<Dna> TSequence;
	typedef StringSet<TSequence, Dependent<> > TStringSet;
	typedef typename Size<typename TMatch::TAlign>::Type TSize;
	typedef Align<TSequence,ArrayGaps> TAlign;      
	typedef typename Iterator<String<TMatch> >::Type TIterator;	 
	

	String<TAlign> alignments;
	
	//conversion to alignments by iterating over the container of matches
	for (TSize j = 0; j < length(matches); j++) 
	{
	
		QueryMatches<TMatch> &queryMatches = value(matches, j);

		TIterator it = begin(queryMatches.matches);
	    TIterator itEnd = end(queryMatches.matches);

	 
		while(it!=itEnd)
		{
		   TAlign align;
	     
			appendValue(align.data_rows, (*it).row1);
			appendValue(align.data_rows, (*it).row2);
			appendValue(alignments, align);

			++it;
	     }

	  }
	
	
	typedef Graph<Alignment<StringSet<TSequence, Dependent<> > > > TAlignGraph;
	TAlignGraph passgraph(dataset);

	//call matchRefinement function, a SeqAn::T-Coffee component
	matchRefinement(alignments,dataset,passgraph);
	
	clear(alignments);
	
	outgraph = passgraph; // passed on to progressiveAlign() function

}

////////////////////////////////////////////////////////////////////////////////


template<typename TStringSet, typename TSequence>
void 
progressiveAlign(Graph<Alignment<StringSet<TSequence, Dependent<> > > > &outgraph, 
					 TStringSet &dataset,
					 const int & minL_
					 )
{
	
	typedef Graph<Alignment<StringSet<TSequence, Dependent<> > > > TAlignGraph;

	//call triplet extension 
	//any two or more STELLAR matches that share a high degree of similarity
	//are connected by the end of this step
	tripletLibraryExtension(outgraph);

	//computation of pairwise distance matrix using k-mer count
	typedef String<double> TDistanceMatrix;
	TDistanceMatrix distanceMatrix;
	getDistanceMatrix(outgraph, distanceMatrix,KmerDistance());
	
	//building a guide tree using UPGMA
	typedef Graph<Tree<double> > TGuideTree;
    TGuideTree  guideTree;
    upgmaTree(distanceMatrix, guideTree);

	TAlignGraph gOut(dataset);


	//call progressive alignment function of SeqAn::T-Coffee
    progressiveAlignment(outgraph, guideTree, gOut);
	
	//call subfunction that extracts motifs
	extractMotifs(gOut,dataset,minL_);
	
	clear(outgraph);
	clear(guideTree);
}


/////////////////////////////////////////////////////////////////////////

template<typename TStringSet,typename TSequence>
void 
extractMotifs(Graph<Alignment<StringSet<TSequence, Dependent<> > > > &gOut, 
			  TStringSet &dataset, 
			  const int & minL
			  )
{
	
	typedef String<Dna> TSequence;
	typedef Graph<Alignment<StringSet<TSequence, Dependent<> > > > TAlignGraph;
	typedef typename VertexDescriptor<TAlignGraph>::Type TVertexDescriptor;
	typedef typename Iterator<TAlignGraph, OutEdgeIterator>::Type TOutEdgeIterator;
	typedef typename Iterator<TAlignGraph, VertexIterator>::Type TVertexIterator;
	typedef typename Iterator<String<TVertexDescriptor>, Rooted>::Type TStringIterator;
	typedef typename Iterator<String<int>, Rooted>::Type TStrIntIterator;
	
	
	TVertexDescriptor nilVertex = getNil<TVertexDescriptor>();
	
	//information about extracted motifs are stored in:
	StringSet<String<Dna> > motifs;
	String<int> motifid;
	String<int> motifposS;
	String<int> motifposE;


	//MERGING CONNECTED COMPONENTS

	//iterating over the refined and aligned alignment graph
	for(TVertexIterator itVertex(gOut);!atEnd(itVertex);++itVertex) 
	{   
		//if no merging is need, extract motif instances right away
		//criterion of 0.7*minL added after first couple of test runs
		if((outDegree(gOut,*itVertex)==(length(dataset)-1)) && (fragmentLength(gOut, *itVertex)>= 0.7*minL)){ 

		   appendValue(motifs, label(gOut,*itVertex));
		   appendValue(motifid, sequenceId(gOut, *itVertex));    
		   appendValue(motifposS, fragmentBegin(gOut, *itVertex));  //starting position in sequence
		   appendValue(motifposE, fragmentLength(gOut,*itVertex));  //end position in sequence
		}
		
		//extract only, if the motif is conserved across all input sequences
		if(outDegree(gOut,*itVertex)==(length(dataset)-1)) 
		{
			
		   //search for an immediate neighbor
		   TVertexDescriptor neighbor = findVertex(gOut, sequenceId(gOut, *itVertex),fragmentBegin(gOut,*itVertex)+fragmentLength(gOut,*itVertex));
		   if (neighbor == nilVertex)  continue; 
		      
		     //to be joined, current vertex and its neighbor should have the same number of outgoing edges
		     if(outDegree(gOut,*itVertex) == outDegree(gOut, neighbor)) 
		     {
				 
			     TOutEdgeIterator e1(gOut, *itVertex);
			     TOutEdgeIterator e2(gOut, neighbor);

				 String<TVertexDescriptor> k;
				 String<TVertexDescriptor> n;

				 appendValue(k, *itVertex); 
				 appendValue(n, neighbor);

				 int countNeighbor = 0;

				 //iteration over outgoing edges of current vertex and its neighbor
				 //targetvertices are stored
			     while(!atEnd(e1), !atEnd(e2))
			     {

				    TVertexDescriptor t1 = targetVertex(e1);
				    TVertexDescriptor t2 = targetVertex(e2);

					appendValue(k, t1);
					appendValue(n, t2);

					//check whether each pair of connected components of current vertex and neighbor
					//are to be found in the same sequence and if they are immediate neighbors
					if(sequenceId(gOut, t1) == sequenceId(gOut, t2) && (fragmentBegin(gOut,t1)+fragmentLength(gOut,t1))== fragmentBegin(gOut,t2))
						++countNeighbor;

					goNext(e1);
				    goNext(e2);

				 }


				 //if the degree of current vertex corresponds to countNeighbor (a means, used to count number of "real" neighbors
				 //store information about its location -> seqId, posStart, posEnd
				 if(outDegree(gOut, *itVertex) == countNeighbor)
				 {
					 String<int> seqIdAndPosA; //seqid & posStart
					 String<int> posE; //posEnd

					 TStringIterator itk = begin(k);
					 TStringIterator itn = begin(n);

					 for(itk,itn;!atEnd(itk),!atEnd(itn);++itk,++itn)
					 {
							 appendValue(seqIdAndPosA, sequenceId(gOut, *itk));    //seqId
							 appendValue(seqIdAndPosA, fragmentBegin(gOut, *itk)); //starting position of the vertex at the far left
							 appendValue(posE, fragmentLength(gOut,*itk)+fragmentLength(gOut,*itn)); //knotenlänge
							 
						
						     removeVertex(gOut,*itk); // vertex at the far left is not required anymore, hence removed
						 
					 }
	
					
					
					 
					//search for further neighbors
					while(findVertex(gOut, sequenceId(gOut, neighbor),fragmentBegin(gOut,neighbor)+fragmentLength(gOut,neighbor))!=nilVertex)
						  
					{
						
						TVertexDescriptor nextneigh = findVertex(gOut, sequenceId(gOut, neighbor),fragmentBegin(gOut,neighbor)+fragmentLength(gOut,neighbor)); 
						
						
						//they too, should share the same number of outgoing edges
						if(outDegree(gOut,neighbor) == outDegree(gOut, nextneigh))
						{
							
							TOutEdgeIterator e3(gOut, neighbor);
							TOutEdgeIterator e4(gOut, nextneigh);

							String<TVertexDescriptor> k1;
							String<TVertexDescriptor> n1;

							appendValue(k1, neighbor);
							appendValue(n1, nextneigh);

							int countNeighbor2 = 0;

							
							while(!atEnd(e3), !atEnd(e4))
							{

								TVertexDescriptor t3 = targetVertex(e3);
								TVertexDescriptor t4 = targetVertex(e4);
								appendValue(k1, t3);
								appendValue(n1, t4);

								if(sequenceId(gOut, t3) == sequenceId(gOut, t4) && (fragmentBegin(gOut,t3)+fragmentLength(gOut,t3))== fragmentBegin(gOut,t4))
									++countNeighbor2;

								goNext(e3);
								goNext(e4);
							}
							
							 
							
							TStringIterator itn1 = begin(n1);
							TStrIntIterator itPosE = begin(posE);

							//if they are real neighbors, position information is stored
							//just the end position (fragment length) of new neighbor is needed
							if(outDegree(gOut, neighbor) == countNeighbor2)
							{
								    
									for(itn1, itPosE;!atEnd(itn1),!atEnd(itPosE);++itn1,++itPosE)
									{
										//add fragment length of new neighbor to the old end position
										assignValue(posE,position(itPosE), *itPosE +fragmentLength(gOut, *itn1));  
										
									}
									
								
								
								neighbor = nextneigh; //set neighbor to new neighbor to search for further neighbors
								
								TStringIterator itk1 = begin(k1);
								for(itk1;!atEnd(itk1);++itk1)
								{
									removeVertex(gOut,*itk1); //remove "old" neighbor
								      
								 }
								
							}else{goto weiter;}




						}else{goto weiter;}

					
						
					}
				
					
					
					weiter:
					//remove neighbor to avoid redundant analysis
					TOutEdgeIterator e5(gOut, neighbor);
					while(!atEnd(e5))
					{
						removeVertex(gOut, targetVertex(e5));
						goNext(e5);
					}
					removeVertex(gOut, neighbor);
					
					
					//iterate over containers that comprise position information
					TStrIntIterator itid = begin(seqIdAndPosA);
					TStrIntIterator iterPosE = begin(posE);
				
						for(itid,iterPosE;!atEnd(itid),!atEnd(iterPosE);itid+=2,++iterPosE)
						{
							TStrIntIterator vertexBegin = itid;
							goNext(vertexBegin);
							
							DnaString mergedVertex = infix(dataset[*itid],*vertexBegin,*vertexBegin+*iterPosE);
					
							if(length(mergedVertex)>=(int)floor(minL*0.7))  //*0.7 added as a preliminary solution for the "gap-problem" 
							  appendValue(motifs, mergedVertex);
							  appendValue(motifid, *itid);
							  appendValue(motifposS, *vertexBegin);
							  appendValue(motifposE, *iterPosE);
						}
					
							
				 }     
				 
			
		     }


		}
	
		
	
	} 
	

	if(length(motifs)==0){

		std::cout<<"No motifs were found!"<<std::endl;
	}
	
	else{

	//if more than one motif is found
	if(length(motifs) % length(dataset) == 0){ 

    int numMotif = length(motifs) /length(dataset); 

	int groupx = 0;

	for(int i = 0; i<numMotif; ++i){


	StringSet<String<Dna> > motifgroup;

    for(int y = groupx; y<length(dataset)+groupx; ++y){

		appendValue(motifgroup, motifs[y]);
		
	}

	//determination of consensus sequence(s)
	consensusSequence(motifgroup);

	groupx = groupx + length(dataset);

	clear(motifgroup);
	
	}

	
	//Key to the outputted consensus sequence
	std::cout<<"\n";
	std::cout<<"\n";
	std::cout<<"KEY"<<std::endl;
	std::cout<<"\n";
	std::cout<<"AC = "<<"Q"<<"\t"<<"AG = "<<"R"<<"\t"<<"AT = "<<"S"<<std::endl;
	std::cout<<"\n";
	std::cout<<"CG = "<<"U"<<"\t"<<"CT = "<<"M"<<"\t"<<"GT = "<<"V"<<std::endl;
	std::cout<<"\n";
	std::cout<<"ACG = "<<"W"<<"\t"<<"ACT = "<<"X"<<"\t"<<"AGT = "<<"Y"<<"\t"<<"CGT = "<<"P"<<std::endl;
	std::cout<<"\n";
	std::cout<<"ACGT = "<<"N"<<std::endl;
	
	}
	else{
	
		//rarely the set of found motifs is uneven and hence there is no determination of
		//a consensus sequence as it is hard to say, which of the entries belong together
		std::cout<<"Length of motifset uneven! No determination of a consensus sequence!"<<std::endl;

		for(int y = 0; y<length(motifs); ++y){

        
		std::cout<<motifs[y]<<std::endl;
		}
	
	}
	
	std::cout<<"****************************************"<<std::endl;

	//output the position information with found motif instances
	getIdMotifsetPos(motifid, motifposS, motifposE, motifs );

	}
	



	}


template <class TFloat, class TInteger>
TInteger round(TFloat& f) {
   return (TInteger)(f + 0.5);
}


//function to determine the consensus sequence in a given set of motifs
void
consensusSequence(StringSet<String<Dna> > motifs){

	//determine the average length of motif out of the found instances
	 
	 float consensusLength = 0;
	 
	 for(int h = 0; h<length(motifs);++h){

		 consensusLength += length(motifs[h]);

	  }

	    consensusLength = consensusLength/length(motifs);

        int consL = round<float,int>(consensusLength);

	  	   
		//Initialization of frequency matrix

	    const int rowsize = 4;  //four rows for A,C,G,T
		const int colsize = consL; 
		
		std::vector< std::vector<double> > matrix;
		matrix.resize( rowsize , std::vector<double>( colsize , 0 ) ); //all matrix entries set to 0

		
	
		//increment matrix entries according to frequency of bases per column

	    for(int j = 0; j<length(motifs);++j){

		  String<char> m = motifs[j];
		  for(int i = 0; i<colsize;++i){
            
			  switch (m[i])
              {
			    case 'A': ++matrix[0][i]; // row 0 -> A
			    break;
			    case 'C': ++matrix[1][i]; // row 1 -> C
                break;
			    case 'G': ++matrix[2][i]; // row 2 -> G
                break;
			    case 'T': ++matrix[3][i]; // row 3 -> T
                break;
			  
			  }
			   
			  
			  
		    }
		  
	      }

		
	 //matrix entries divided by number of sequences
	  for(int r=0;r<4;++r){ 

		 for(int c=0; c<colsize;++c){  
		  
			 matrix[r][c]= (matrix[r][c])/length(motifs); 
		   
		 }
	   }


	  //extraction of most frequent base in column
	  //entries with identical values are considered as well

	  int base1 = 0;
	  int base2 = 0;
	  int join1 = 0;
	  int join2 = 0;
	  int add1 = 0;
	  int add2 = 0;

	  String<int> bases;


	  for(int c=0; c<colsize;++c){

	    if((matrix[0][c]>matrix[1][c]) && (matrix[0][c]!=matrix[1][c])){  //A or C
		    base1 = 0;} else{base1 = 1;}
	  
		  
        if((matrix[2][c]>matrix[3][c]) && (matrix[2][c]!=matrix[3][c])){  //G or T
		    base2= 2;}else{base2 = 3;}

	  
		if(matrix[0][c] == matrix[1][c]){   // A and C
			  base1 = 0;
			  join1= 10;}else{join1 = 0;}

      if(matrix[2][c] == matrix[3][c]){   // G and T
		  base2 = 2;
		  join2 = 23;}else{join2 = 0;}

	  

      if((matrix[base1][c]>matrix[base2][c]) && (matrix[base1][c] != matrix[base2][c])){ 
		 
		  if((base1!=1) && (join1 > 0)){  //AC
			  add1= 10;} 
		  else if(base1 == 1){add1 = 2;} //C
		  else{add1 = 1;}
	  }
	 else if((matrix[base1][c]<matrix[base2][c]) && (matrix[base1][c] != matrix[base2][c])){
		  
		  if((base2!=3) && (join2>0)){    //GT
			  add1= 34;}    
		  else if(base2 == 3){add1=4;}  //T
		  else if((base2 == 2) && (join2 == 0)){add1 = 3;} //G

	  }
	 else{add1= 0;}

	 


	  if(matrix[base1][c] == matrix[base2][c]){

		  if((base1 == 1) && (base2== 3)){ //CT
			  add2 = 24;}
		  else if(((base1 == 0) && (join1>0)) &&(base2==2 && join2>0)){ //ACGT
			  add2 = 1234;}
		  else if((base1 == 0 && join1 == 0)&&(base2==2 && join2>0)){ //AGT
			  add2 = 134;}
		  else if((base1 == 0 && join1 == 0) && ( base2==2 && join2 == 0)){ //AG
			  add2 = 13;}
		  else if((base1 == 0 && join1>0) && (base2==2 && join2 ==0)){ //ACG
			  add2 = 123;}
		  else if((base1 == 0 && join1 == 0) && (base2 == 3)){ //AT
			  add2 = 14;}
		  else if(base1 == 1 && (base2 ==2 && join2 == 0)){ //CG
			  add2 = 23;}
		  else if((base1 == 0 && join1> 0) && (base2 == 3)){ //ACT
			  add2 = 124;}
		  else if((base1 == 1) && (base2 == 2 && join2>0)){ //CGT
			  add2 = 234;}
		  else{add2 = 0;}

	  }
		  

      
	   if(add1>0){

		   appendValue(bases, add1);
		   
	   
	   }else{appendValue(bases, add2);}





	  }


	   typedef Iterator<String<int>, Rooted>::Type TStrIntIterator;
	   TStrIntIterator itba = begin(bases);
	   String<char> consensusSeq;

	   while(!atEnd(itba)){
	 

       //convert back to actual bases
	   switch (*itba)
            {
			  case  1: appendValue(consensusSeq, 'A');
			  break;
			  case  2: appendValue(consensusSeq, 'C');
              break;
			  case  3: appendValue(consensusSeq, 'G');
              break;
			  case  4: appendValue(consensusSeq, 'T');
              break;
			  case  10: appendValue(consensusSeq, 'Q'); //AC
			  break;
			  case  34: appendValue(consensusSeq, 'V'); //GT
              break;
			  case  24: appendValue(consensusSeq, 'M'); //CT
              break;
			  case  1234: appendValue(consensusSeq, 'N'); //ACGT
              break;
			  case  134: appendValue(consensusSeq, 'Y'); //AGT
			  break;
			  case  13: appendValue(consensusSeq, 'R'); //AG
              break;
			  case  123: appendValue(consensusSeq, 'W'); //ACG
              break;
			  case  14: appendValue(consensusSeq, 'S'); //AT
              break;
			  case  23: appendValue(consensusSeq, 'U');  //CG
			  break;
			  case  124: appendValue(consensusSeq, 'Q'); //ACT
              break;
			  case  234: appendValue(consensusSeq, 'P'); //CGT
              break;
			 
			  
			}
			 


	  goNext(itba);
	 }


     //output the resulting consensus sequence
     std::cout<<"\n";
	 std::cout<<"Consensus sequence:  "<< consensusSeq<<"  "<<"length: "<<length(consensusSeq)<<std::endl;
	 std::cout<<"\n";



	 
}



//function to obtain position information and each found motif instance
void
getIdMotifsetPos(String<int>motifid, String<int>motifPosStart, String<int>motifPosEnd, StringSet<String<Dna> > motifs){


	typedef Iterator<StringSet<String<Dna>>, Rooted>::Type TStrSetIter;
	TStrSetIter itmotiv = begin(motifs);

	typedef Iterator<String<int>, Rooted>::Type TStrIntIterator;

	TStrIntIterator itsId = begin(motifid);
	TStrIntIterator itposStart = begin(motifPosStart);
    TStrIntIterator itposEnd = begin(motifPosEnd);

	//iterate over containers where the position information and motifs were stored
	while(!atEnd(itsId),!atEnd(itposStart),!atEnd(itposEnd), !atEnd(itmotiv)){

		//output motif instances and their locations in the set of sequences
		std::cout<<"seqId: "<<*itsId<<"\t"<<*itmotiv<<"\t"<<"location:"<<*itposStart<<" - "<<*itposStart+*itposEnd<<std::endl;


		goNext(itsId);
		goNext(itposStart);
		goNext(itposEnd);
		goNext(itmotiv);

	}


	

}



#endif  // #ifndef SANDBOX_SHALINI6_APPS_MOTIFSEARCH_MOTIFHEADER_H_
