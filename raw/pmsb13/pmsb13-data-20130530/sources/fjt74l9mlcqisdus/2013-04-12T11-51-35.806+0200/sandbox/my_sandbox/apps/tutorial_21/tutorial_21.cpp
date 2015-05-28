#include <iostream>
#include <seqan/graph_types.h>
#include <seqan/graph_algorithms.h>
#include <seqan/find_motif.h>
#include <seqan/basic.h>


using namespace seqan;

int main ()
{
	typedef unsigned int TCargo;
    typedef Graph<Hmm<Dna, LogProb<double>, Default> > TGraph;
    typedef VertexDescriptor<TGraph>::Type TVertexDescriptor;
	TGraph g;
	TVertexDescriptor exon = addVertex(g);
	TVertexDescriptor splice = addVertex(g);
	TVertexDescriptor intron = addVertex(g);
	TVertexDescriptor start = addVertex(g);
	TVertexDescriptor end = addVertex(g);
	assignBeginState(g,start);
	assignEndState(g,end);

	addEdge(g,exon,exon,0.9);
	addEdge(g,exon,splice,0.1);
	addEdge(g,splice,intron,1);
	addEdge(g,intron,intron,0.9);
	
	addEdge(g,start,exon,1);
	addEdge(g,intron,end,0.1);
	addEdge(g,end,end,1);
	
	emissionProbability(g,exon,(Dna)'A')=0.25;
	emissionProbability(g,exon,(Dna)'C')=0.25;
	emissionProbability(g,exon,(Dna)'G')=0.25;
	emissionProbability(g,exon,(Dna)'T')=0.25;
   	
	
	emissionProbability(g,splice,(Dna)'A')=0.05;
	emissionProbability(g,splice,(Dna)'C')=0.0;
	emissionProbability(g,splice,(Dna)'G')=0.95;
	emissionProbability(g,splice,(Dna)'T')=0.0;

	emissionProbability(g,intron,(Dna)'A')=0.4;
	emissionProbability(g,intron,(Dna)'C')=0.1;
	emissionProbability(g,intron,(Dna)'G')=0.1;
	emissionProbability(g,intron,(Dna)'T')=0.4;
	
	String<TVertexDescriptor> path;
	String<Dna> s = "CTTCATGTGAAAGCAGACGTAAGTCA";
	std::cout << viterbiAlgorithm(g,s,path)<<std::endl;
	
	//::std::cout << g << ::std::endl;
	return 0;
}