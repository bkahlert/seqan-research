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

using namespace seqan;


//Sequenzen aus Fasta-Datei in StringSet einlesen

template<typename TStringSet, typename TId>
void readSequences(CharString &filename, TStringSet &seqs, StringSet<TId> &seqIDs){

	std::cout<<"read sequences..........."<<std::endl;
	typedef String<Dna> TSequence;
	typedef CharString TId;
	typedef StringSet<TSequence> TStringSet;

	MultiSeqFile multiSeqFile;
	if (!open(multiSeqFile.concat, toCString(filename), OPEN_RDONLY))
	    {
	        std::cerr << "Failed to open " << filename << " file." << std::endl;
       
	    }

    AutoSeqFormat format;
	guessFormat(multiSeqFile.concat, format);
	split(multiSeqFile, format);

	unsigned seqCount = length(multiSeqFile);

	reserve(seqs, seqCount, Exact());
	reserve(seqIDs, seqCount, Exact());

	TSequence seq;
	TId seqid;

	for (unsigned i = 0; i < seqCount; ++i){

		assignSeq(seq, multiSeqFile[i], format);   
		assignSeqId(seqid, multiSeqFile[i], format);   

		appendValue(seqs, seq, Generous());
		appendValue(seqIDs, seqid, Generous());
	}


}


//STELLAR-Aufruf

template<typename TSequence, typename TStringSet, typename TMatch>
void callStellar(TSequence &query, TStringSet &remainingSeqs, StringSet<QueryMatches<TMatch > > &matches){

	typedef String<Dna> TSequence;
	typedef CharString TId;
	typedef StringSet<TSequence, Dependent<> > TStringSet;

	//Stellarparameter
	
	//motivlänge in fasta dateien = 30

	double epsilon = 0.2; //3seq0err.fa = 0.01 | 4seq0err.fa = 0.01 | 4seq1err.fa = 0.04 |4seq5err.fa = 0.2| seqF22.fa = 0.01
	int minLength = 25;  // 3seq0err.fa = 20 | 4seq0err.fa = 20 | 4seq1err.fa = 20 |4seq5err.fa = 25
	unsigned xDrop = 5;  // bei allen fasta dateien 5
	
	//Berechnung von smin (s. stellar.cpp Z.320-324)
	int errMinLen = (int) floor(epsilon * minLength);
	int n = (int) ceil((errMinLen + 1) / epsilon);
	int errN = (int) floor(epsilon * n);
	unsigned smin = (unsigned) _min(ceil((double)(minLength-errMinLen)/(errMinLen+1)),
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

	//resize(matches, length(remainingSeqs));

	//Stellaraufruf
	stellar(sFinder, sPattern, epsilon, minLength, xDrop, matches, AllLocal());

}


//matchrefinement

template<typename TMatch, typename TStringSet, typename TSequence>
void refinementAndAlign(StringSet<QueryMatches<TMatch > > &matches, TStringSet &seqs, 
						Graph<Alignment<StringSet<TSequence, Dependent<> > > > &passgraph){

   
	typedef String<Dna> TSequence;
	typedef StringSet<TSequence, Dependent<> > TStringSet;
	typedef Size<typename TMatch::TAlign>::Type TSize;
	typedef Align<TSequence,ArrayGaps> TAlign;      
	typedef Iterator<String<TMatch> >::Type TIterator;	 
	

	String<TAlign> alignments;
	std::cout<<"start:"<<std::endl;
	
	for (TSize j = 0; j < length(matches); j++) {
	
		QueryMatches<TMatch> &queryMatches = value(matches, j);

		TIterator it = begin(queryMatches.matches);
	    TIterator itEnd = end(queryMatches.matches);

	 
		while(it!=itEnd){

			TAlign align;
	
			appendValue(align.data_rows, (*it).row1);
			appendValue(align.data_rows, (*it).row2);
			//std::cout<<"id1: "<<getObjectId(source((*it).row1))<<", id2: "<<getObjectId(source((*it).row2))<<std::endl;
			//std::cout<<align<<std::endl;
			appendValue(alignments, align);

			++it;
	     }

	  }
	
	
	
	typedef Graph<Alignment<StringSet<TSequence, Dependent<> > > > TAlignGraph;
	TAlignGraph outgraph(seqs);

	matchRefinement(alignments,seqs,outgraph);

	//std::cout<<outgraph<<std::endl;
	
	clear(alignments);
	
	
	passgraph = outgraph;

	
}


//progressive alignment

template<typename TStringSet, typename TSequence, typename TId>
void progressiveAlignment(Graph<Alignment<StringSet<TSequence, Dependent<> > > > &passgraph, TStringSet &allSeqs,
						  StringSet<TId> &seqIDs){

    std::cout<<"progressive alignment..........."<<std::endl;

	typedef Graph<Alignment<StringSet<TSequence, Dependent<> > > > TAlignGraph;

	tripletLibraryExtension(passgraph);

	typedef String<double> TDistanceMatrix;
	TDistanceMatrix distanceMatrix;
	getDistanceMatrix(passgraph, distanceMatrix,KmerDistance());
	
	typedef Graph<Tree<double> > TGuideTree;
    TGuideTree  guideTree;
    upgmaTree(distanceMatrix, guideTree);

	TAlignGraph gOut(allSeqs);
    progressiveAlignment(passgraph, guideTree, gOut);
	
	//std::cout<<gOut<<std::endl;
	
	
	
	////////////////////////////////////////////////////////////////////
	//                  MOTIF - EXTRACTION
	////////////////////////////////////////////////////////////////////
	 std::cout<<"motif extraction..........."<<std::endl;
	typedef typename VertexDescriptor<TAlignGraph>::Type TVertexDescriptor;
	typedef typename Iterator<TAlignGraph, OutEdgeIterator>::Type TOutEdgeIterator;
	typedef typename Iterator<TAlignGraph, VertexIterator>::Type TVertexIterator;
	typedef typename Iterator<String<TVertexDescriptor>, Rooted>::Type TStringIterator;
	typedef typename Iterator<String<int>, Rooted>::Type TStrIntIterator;
	
	
	TVertexDescriptor nilVertex = getNil<TVertexDescriptor>();
	
	
	std::cout<<"original graph:"<<std::endl;
	std::cout<<"# vertices: "<<numVertices(gOut)<<std::endl;
	std::cout<<"# edges: "<<numEdges(gOut)<<std::endl;
	std::cout<<"............................"<<std::endl;
	

	//KONSEKUTIVE KNOTEN ZUSAMMENFÜHREN

	
	//int r = 0;
	//gehe jeden knoten im graphen durch
	for(TVertexIterator itVertex(gOut);!atEnd(itVertex);++itVertex) 
	{   
		
		if(outDegree(gOut,*itVertex)>=1) //es muss mindestens eine Kante ausgehen
		{
			
		   //suche nach dem direkten nachbarknoten in einer sequence
		   TVertexDescriptor neighbor = findVertex(gOut, sequenceId(gOut, *itVertex),fragmentBegin(gOut,*itVertex)+fragmentLength(gOut,*itVertex));
		   if (neighbor == nilVertex)  continue; 
		  
		     //wenn vom aktuellen knoten und dem nachbarknoten die gleiche kantenanzahl ausgeht, dann ....
		     if(outDegree(gOut,*itVertex) == outDegree(gOut, neighbor))
		     {
			     TOutEdgeIterator e1(gOut, *itVertex);
			     TOutEdgeIterator e2(gOut, neighbor);

				 String<TVertexDescriptor> k;
				 String<TVertexDescriptor> n;

				 appendValue(k, *itVertex); //merke knoten und seinen nachbarknoten
				 appendValue(n, neighbor);

				 int countNeighbor = 0;

			     //iteriere über ausgehende kanten des aktuellen knoten und seinem nachbar
				 //und merke mir diese zielknoten
			     while(!atEnd(e1), !atEnd(e2))
			     {

				    TVertexDescriptor t1 = targetVertex(e1);
				    TVertexDescriptor t2 = targetVertex(e2);

					appendValue(k, t1);
					appendValue(n, t2);

					// wenn alle Knoten, also knoten, sein nachbarknoten sowie ihre jeweiligen zielknoten immer jeweils in derselben sequenz 
					// bzw nebeneinander liegen, dann wird countNeighbor hochgezählt
					if(sequenceId(gOut, t1) == sequenceId(gOut, t2) && (fragmentBegin(gOut,t1)+fragmentLength(gOut,t1))== fragmentBegin(gOut,t2))
						++countNeighbor;

					goNext(e1);
				    goNext(e2);

				 }

				 //std::cout<<"runde: "<<r<<std::endl;
				
				 //wenn Kantenausgangszahl == countNeighbor entspricht, d.h. wenn wirkliche alle Knoten, die obige Bedingung erfüllen
				 //dann wird die sequenzid, Anfangsposition vom 'linkesten' Knoten
				 //und  akutelle Länge des knoten+nachbarknoten gespeichert 
				 if(outDegree(gOut, *itVertex) == countNeighbor)
				 {
					 String<int> seqIdAndPosA; //seqid & posanfang
					 String<int> posE; //posEnde bzw totale Länge -> für addVertex später...

					 TStringIterator itk = begin(k);
					 TStringIterator itn = begin(n);

					 for(itk,itn;!atEnd(itk),!atEnd(itn);++itk,++itn)
					 {
							 appendValue(seqIdAndPosA, sequenceId(gOut, *itk));    //seqid
							 appendValue(seqIdAndPosA, fragmentBegin(gOut, *itk)); //knotenanfang
			
							 appendValue(posE, fragmentLength(gOut,*itk)+fragmentLength(gOut,*itn)); //knotenlänge
							 
						
						 
						     removeVertex(gOut,*itk); // der linkeste knoten wird nicht mehr benötigt
						 
					 }
	
					
					
					 
					//weitere nachbarknoten suchen
					while(findVertex(gOut, sequenceId(gOut, neighbor),fragmentBegin(gOut,neighbor)+fragmentLength(gOut,neighbor))!=nilVertex)
						  
					{
						
						TVertexDescriptor nextneigh = findVertex(gOut, sequenceId(gOut, neighbor),fragmentBegin(gOut,neighbor)+fragmentLength(gOut,neighbor)); 
						
						
						//gucken, ob gleiche Kantenanzahl ausgeht 
						//wenn ja merke wie oben alle knoten und ihre zielknoten
						if(outDegree(gOut,neighbor) == outDegree(gOut, nextneigh))
						{
							
							TOutEdgeIterator e3(gOut, neighbor);
							TOutEdgeIterator e4(gOut, nextneigh);

							String<TVertexDescriptor> k1;
							String<TVertexDescriptor> n1;

							appendValue(k1, neighbor);
							appendValue(n1, nextneigh);

							int countNeighbor2 = 0;

							

							
							//iteriere über ausgehende kanten des aktuellen knoten (neighbor) und seinem (neuen) nachbar
							while(!atEnd(e3), !atEnd(e4))
							{

								TVertexDescriptor t3 = targetVertex(e3);
								TVertexDescriptor t4 = targetVertex(e4);
								appendValue(k1, t3);
								appendValue(n1, t4);

								//wenn in gleicher sequenz und nebeneinander countNeighbor hochzählen
								if(sequenceId(gOut, t3) == sequenceId(gOut, t4) && (fragmentBegin(gOut,t3)+fragmentLength(gOut,t3))== fragmentBegin(gOut,t4))
									++countNeighbor2;

								goNext(e3);
								goNext(e4);
							}
							
							 
							
							TStringIterator itn1 = begin(n1);
							TStrIntIterator itPosE = begin(posE);

							//wenn wirklich alle nebeneinander
							if(outDegree(gOut, neighbor) == countNeighbor2)
							{
								    
									for(itn1, itPosE;!atEnd(itn1),!atEnd(itPosE);++itn1,++itPosE)
									{
										//alte Knotenlänge plus Knotenlänge vom neuen Nachbarknoten
										assignValue(posE,position(itPosE), *itPosE +fragmentLength(gOut, *itn1));  
										
									}
									
								
								
								neighbor = nextneigh; //neuer Nachbarknoten(nextneigh) ist nun aktueller Knoten für den ein weiterer Nachbar gesucht wird
								
								TStringIterator itk1 = begin(k1);
								for(itk1;!atEnd(itk1);++itk1)
								{
									removeVertex(gOut,*itk1); //lösche linken nachbarn vom aktuellen knoten( alter neighbor)
								      
								 }
								
							}else
							 {  /*
								TStringIterator itk1 = begin(k1);
								for(itk1;!atEnd(itk1);++itk1)
								{
									//lösche linken nachbarn vom aktuellen knoten( alter neighbor) und beende suche nach weiteren nachbarknoten
									removeVertex(gOut,*itk1); 
								}
								break;*/goto weiter;
								
							 }




						}else{goto weiter;} //wenn erst gar kein Nachbar gefunden, der dieselbe Kantenanzahl hat dann neighbor löschen, damit später addVertex funzt

					
						
					}
				
					
					
					weiter:
					TOutEdgeIterator e5(gOut, neighbor);
					while(!atEnd(e5))
					{
						removeVertex(gOut, targetVertex(e5));
						goNext(e5);
					}
					removeVertex(gOut, neighbor);
					


					//konsekutive knoten in graphen einfügen
					
					TStrIntIterator itid = begin(seqIdAndPosA);
					TStrIntIterator iterPosE = begin(posE);

					
					if(*iterPosE >=30) //30 ist MotifLänge
					{						
						for(itid,iterPosE;!atEnd(itid),!atEnd(iterPosE);itid+=2,++iterPosE)
						{
							TStrIntIterator vertexBegin = itid;
							goNext(vertexBegin);
							//std::cout<<*itid<<" : "<<*vertexBegin<<" - "<<*iterPosE<<std::endl;
							addVertex(gOut,*itid,*vertexBegin,*iterPosE); //iterPosE = length of new vertex

							
						}
					}
							
				 }     
				 
			
		     }


		}
	
    //++r;
	
	
	} 

	//std::cout<<gOut<<std::endl;
	std::cout<<"changed graph:"<<std::endl;
	std::cout<<"# vertices: "<<numVertices(gOut)<<std::endl;
	std::cout<<"# edges: "<<numEdges(gOut)<<std::endl;
	

	
	String<TVertexDescriptor > motive;
	// potentielle Motive heraussuchen
	for(TVertexIterator itMotive(gOut);!atEnd(itMotive);++itMotive) 
	{  
		//std::cout<<label(gOut, *itMotive)<<"  deg: "<<outDegree(gOut, *itMotive)<<std::endl;
		//nach konsekutive-knoten-suche zusammengeführte nachbarknoten zwar als knoten eingefügt, aber keine kanten eingfügt
		//deshalb sind potentielle knoten mit kantenanzahl = 0 und knotenlänge = motivlänge  
		//-> aber klappt manchmal nicht, da durch erlaubte mismatches/indels ein in der sequenz verstecktes 30-basenlanges motiv zu 32b wird??
		//oder knoten von denen sequenzanzahl-1 kanten ausgehen und Knotenlänge = motivlänge
		//-> klappt nicht immer, weil stellar matches und anschließendes matchrefinement nicht alle motive alignieren ??
		//deshalb fürs erste :  30<= Knotenlänge <=35    als Bedingung
		if((outDegree(gOut, *itMotive)==0 && (fragmentLength(gOut, *itMotive)>=30 && fragmentLength(gOut, *itMotive)<=35))
			|| (outDegree(gOut, *itMotive)>=1 && (fragmentLength(gOut, *itMotive)>=30 && fragmentLength(gOut, *itMotive)<=35))){ //outDegree(gOut, *itMotive)>= anzahl der sequenzen -1 -> bedeutet, dass das Motif in mehreren sequenzen vorkommt
			//std::cout<<"id: "<<sequenceId(gOut,*itMotive)<<" "<<label(gOut, *itMotive)<<" "<<"deg: "<<outDegree(gOut, *itMotive)<<" "<<fragmentBegin(gOut, *itMotive)<<"-"<<fragmentLength(gOut, *itMotive)<<std::endl;
			//appendValue(motive, label(gOut,*itMotive));
		appendValue(motive, *itMotive);}
		
	}

	
	typedef typename Iterator<String<TVertexDescriptor>, Rooted>::Type TDnaIterator;

	TDnaIterator itmotif = begin(motive);
	String<TVertexDescriptor > lmers;
	appendValue(lmers, *itmotif);
	
	
	for(itmotif ;!atEnd(itmotif);++itmotif)
	{

		int duplicate = 0;
		TDnaIterator itlmer = itmotif;
	
	
	}
	


		
	clear(passgraph);
	clear(guideTree);
}



int main(int argc, char const ** argv){


	if(argc!=2){

			std::cerr << "Error: Too many arguments"<<std::endl
					  << "USAGE: " << argv[0] << " FILE" << std::endl;

		    return 1;
	}


	typedef String<Dna> TSequence;
	typedef CharString TId;
	typedef StringSet<TSequence> TStringSet;
	
	
	CharString filename = argv[1];
	TStringSet seqs;
	StringSet<TId> seqIDs;


	//Sequenzen werden eingelesen
	readSequences(filename, seqs, seqIDs);

	
	//zum Speichern der Matches
	typedef StellarMatch<TSequence, TId> TMatch;
	StringSet<QueryMatches<TMatch > > matches;
	resize(matches, length(seqs));
	
	//alle Sequenzen durchgehen und jeweils eine(query) gegen den Rest(remainingSeqs) prüfen
	std::cout<<"call stellar..........."<<std::endl;
    for(unsigned r = 0; r<length(seqs)-1; ++r)
	{

		TSequence &query = seqs[r];
		StringSet<TSequence, Dependent<> > remainingSeqs;

		for(unsigned s = r+1; s<length(seqs);++s)
		{

			appendValue(remainingSeqs, seqs[s]);
		}

	
		//STELLAR-Aufruf mit eingegebenen Sequenzen
		callStellar(query,remainingSeqs, matches);

		
	 }

		//matchrefinement
	    std::cout<<"match refinement..........."<<std::endl;
		typedef Graph<Alignment<StringSet<TSequence, Dependent<> > > > TAlignGraph;
		TAlignGraph passgraph(seqs);
		clear(passgraph);

		refinementAndAlign(matches, seqs, passgraph); 
			
	    //progressives Alignement
		
		progressiveAlignment(passgraph, seqs, seqIDs);


    return 0;
}