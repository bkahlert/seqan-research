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
	TVertexDescriptor s_t = addVertex(g);

	addEdge(g,s_t,exon,1);
	addEdge(g,exon,exon,0.9);
	addEdge(g,exon,splice,0.1);
	addEdge(g,splice,intron,0.9);
	addEdge(g,intron,intron,0.9);
	addEdge(g,intron,s_t,0.1);

	emissionProbability(g,exon,'A')=0.25;
	emissionProbability(g,exon,'C')=0.25;
	emissionProbability(g,exon,'G')=0.25;
	emissionProbability(g,exon,'T')=0.25;
   	
	
	emissionProbability(g,splice,'A')=0.05;
	emissionProbability(g,splice,'C')=0.0;
	emissionProbability(g,splice,'G')=0.95;
	emissionProbability(g,splice,'T')=0.0;

	emissionProbability(g,intron,'A')=0.4;
	emissionProbability(g,intron,'C')=0.1;
	emissionProbability(g,intron,'G')=0.1;
	emissionProbability(g,intron,'T')=0.4;
	::std::cout << g << ::std::endl;
	return 0;
}