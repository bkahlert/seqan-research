#include <iostream>
#include <seqan/graph_types.h>
#include <seqan/graph_algorithms.h>
#include <seqan/find_motif.h>
#include <seqan/basic.h>


using namespace seqan;

int main ()
{
	typedef unsigned int TCargo;
    typedef Graph<Undirected<TCargo> > TGraph;
    typedef VertexDescriptor<TGraph>::Type TVertexDescriptor;
	TGraph g;
	TVertexDescriptor exon = addVertex(g);
	TVertexDescriptor splice = addVertex(g);
	TVertexDescriptor intron = addVertex(g);
	TVertexDescriptor s_t = addVertex(g);

	addEdges(g, s_t, exon, LogProb<double> 1);
	addEdges(g, exon, exon, LogProb<double> 0.9);
	addEdges(g, exon, splice,  LogProb<double> 0.1);
	addEdges(g, splice, intron,  LogProb<double> 0.9);
	addEdges(g, intron, intron,  LogProb<double> 0.9);
	addEdges(g, intron, s_t,  LogProb<double> 0.1);
   	::std::cout << g << ::std::endl;
	return 0;
}