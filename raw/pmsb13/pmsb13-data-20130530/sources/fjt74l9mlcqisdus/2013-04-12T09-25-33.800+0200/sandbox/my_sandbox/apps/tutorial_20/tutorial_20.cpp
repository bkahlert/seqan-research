#include <iostream>
#include <seqan/graph_types.h>
#include <seqan/graph_algorithms.h>
using namespace seq

int main ()
{
    typedef unsigned int TCargo;
    typedef Graph<Undirected<TCargo> > TGraph;
    typedef VertexDescriptor<TGraph>::Type TVertexDescriptor;
	TGraph g;
	TVertexDescriptor vertBerlin = addVertex(g);
    TVertexDescriptor vertHamburg = addVertex(g);
    TVertexDescriptor vertHannover = addVertex(g);
    TVertexDescriptor vertMainz = addVertex(g);
    TVertexDescriptor vertMuenchen = addVertex(g);
	 addEdge(g, vertBerlin, vertHamburg, 289);
    addEdge(g, vertBerlin, vertHannover, 286);
    addEdge(g, vertBerlin, vertMainz, 573);
    addEdge(g, vertBerlin, vertMuenchen, 586);
    addEdge(g, vertHannover, vertMuenchen, 572);
    addEdge(g, vertHamburg, vertMainz, 521);
	FILE* strmWrite = fopen("graph.dot", "w");
    write(strmWrite, g, DotDrawing());
    fclose(strmWrite);
	return 0;
}