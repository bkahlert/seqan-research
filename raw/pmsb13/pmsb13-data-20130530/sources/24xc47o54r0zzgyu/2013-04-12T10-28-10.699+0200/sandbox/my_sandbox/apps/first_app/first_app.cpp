#include <iostream>
#include <seqan/graph_types.h>
#include <seqan/graph_algorithms.h>
#include <seqan/basic.h>

using namespace seqan;

int main ()
{
	typedef LogProb<> TCargo;
	typedef Graph<Directed<TCargo> > TGraph;
	typedef VertexDescriptor<TGraph>::Type TVertexDescriptor;

	TGraph g;

	TVertexDescriptor start = addVertex(g);
	TVertexDescriptor exon = addVertex(g);
	TVertexDescriptor intron = addVertex(g);
	TVertexDescriptor splice = addVertex(g);
	TVertexDescriptor end = addVertex(g);
	

	addEdge(g, start, exon, (TCargo)1.0);
	addEdge(g, exon, exon, (TCargo)0.9);
	addEdge(g, exon, splice, (TCargo)0.1);
	addEdge(g, splice, intron, (TCargo)1.0);
	addEdge(g, intron, intron, (TCargo)0.9);
	addEdge(g, intron, end, (TCargo)0.1);
	
	::std::cout << g << ::std::endl;

    
    return 0;
}