#include <iostream>
#include <seqan/graph_types.h>
#include <seqan/graph_algorithms.h>
using namespace seqan;

int main (){

	//typedef unsigned int TCargo;
 //   typedef Graph<Directed<TCargo> > TGraph;
 //   typedef VertexDescriptor<TGraph>::Type TVertexDescriptor;
	//TGraph g;
	//TVertexDescriptor edges[]={1,0,0,4,2,1,4,1,5,1,6,2,3,2,2,3,7,3,5,4,6,5,5,6,7,6,7,7};
	//addEdges(g,edges,14);
	////FILE* strmWrite = fopen("graph.dot", "w");
 ////   write(strmWrite, g, DotDrawing());
 ////   fclose(strmWrite);
	////::std::cout << g << ::std::endl;

	//typedef String<char> TCityName;
	//typedef String<TCityName> TProperties;
	//TProperties cityNames;
	//resizeVertexMap(g,cityNames);
	//assignValue(cityNames,0,"a");
	//assignValue(cityNames,1,"b");
	//assignValue(cityNames,2,"c");
	//assignValue(cityNames,3,"d");
	//assignValue(cityNames,4,"e");
	//assignValue(cityNames,5,"f");
	//assignValue(cityNames,6,"g");
	//assignValue(cityNames,7,"h");
	///*
	//String<char> nameMap;
 //   char names[] = {'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h'};
 //   assignVertexMap(g,nameMap, names);

	//*/

	//::std::cout << g << ::std::endl;
	//::std::cout << cityNames << ::std::endl;

	//typedef Iterator<TGraph, DfsPreorder>::Type TDfsIterator;
	//TVertexDescriptor start;
	//for(int i =0;i<length(cityNames);++i){
	//start=i;
 //   TDfsIterator itV(g,start);
 //   for(;!atEnd(itV);goNext(itV)) {
 //       ::std::cout << value(itV) << ':' << getProperty(cityNames, value(itV)) << ::std::endl;
 //   }
	//::std::cout<<::std::endl;
	//}
	//String<TCargo> component;



	//stronglyConnectedComponents(g,component);
	//typedef Iterator<TGraph, VertexIterator>::Type TVertexIterator;
	//TVertexIterator TvI(g);
	//while(!atEnd(TvI)){

	//	::std::cout<<getProperty(cityNames,value(TvI))<<" ";
	//	::std::cout<<getProperty(component,value(TvI))<<::std::endl;

	//	++TvI;
	//}

	
	typedef Dna TAlphabet;
	typedef Size<TAlphabet>::Type TSize;
	typedef Graph<Hmm<TAlphabet,double,Default> > THmm;
	typedef VertexDescriptor<THmm>::Type TVertexDescriptor;
	typedef EdgeDescriptor<THmm>::Type TEdgeDescriptor;

	Dna dnaA = Dna('A');
	Dna dnaC = Dna('C');
	Dna dnaG = Dna('G');
	Dna dnaT = Dna('T');

	THmm hmm;

	TVertexDescriptor begState = addVertex(hmm);
	assignBeginState(hmm,begState);
	

	TVertexDescriptor exon =addVertex(hmm);
	assignEmissionProbability(hmm,exon,dnaA,0.25);
	assignEmissionProbability(hmm,exon,dnaC,0.25);
	assignEmissionProbability(hmm,exon,dnaT,0.25);
	assignEmissionProbability(hmm,exon,dnaG,0.25);
	

	TVertexDescriptor splice = addVertex(hmm);
	assignEmissionProbability(hmm,splice,dnaA,0.05);
	assignEmissionProbability(hmm,splice,dnaC,0.0);
	assignEmissionProbability(hmm,splice,dnaT,0.0);
	assignEmissionProbability(hmm,splice,dnaG,0.95);
	
	
	TVertexDescriptor intron =addVertex(hmm);
	assignEmissionProbability(hmm,intron,dnaA,0.4);
	assignEmissionProbability(hmm,intron,dnaC,0.1);
	assignEmissionProbability(hmm,intron,dnaT,0.4);
	assignEmissionProbability(hmm,intron,dnaG,0.1);
	

	
	TVertexDescriptor endState = addVertex(hmm);
	assignEndState(hmm,endState);
	
	/*assignTransitionProbability(hmm,begState,exon,1); --> benötigt vorher eine kante
	assignTransitionProbability(hmm,exon,exon,0.9);
	assignTransitionProbability(hmm,exon,splice,0.1);
	assignTransitionProbability(hmm,splice,intron,1);
	assignTransitionProbability(hmm,intron,intron,0.9);
	assignTransitionProbability(hmm,intron,endState,0.1);*/
	addEdge(hmm,begState,exon,1);
	addEdge(hmm,exon,exon,0.9);
	addEdge(hmm,exon,splice,0.1);
	addEdge(hmm,splice,intron,1);
	addEdge(hmm,intron,intron,0.9);
	addEdge(hmm,intron,endState,0.1);

	::std::cout<<hmm<<::std::endl;

	String<Dna> Seq;
	String<TVertexDescriptor> Path;
	Seq="CTTCATGTGAAAGCAGACGTAAGTCA";
	LogProb<> P =viterbiAlgorithm(hmm,Seq,Path);
	
	LogProb<> P1 =forwardAlgorithm(hmm,Seq);
	LogProb<> P2 =backwardAlgorithm(hmm,Seq);
	::std::cout<<P<<" "<<Path<<::std::endl;
	::std::cout<<P1<<" "<<P2;

	return 0;

}