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

	addEdge(g,s_t,exon,1);
	addEdge(g,exon,exon,0.9);
	addEdge(g,exon,0.1);
	addEdge(g,splice,intron,0.9);
	addEdge(g,intron,intron,0.9);
	addEdge(g,intron,s_t,0.1);
   	::std::cout << g << ::std::endl;
	return 0;
}