#include <iostream>
#include <seqan/graph_types.h>
#include <seqan/graph_algorithms.h>

using namespace seqan;

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

	addEdge(g, vertBerlin, vertHamburg, (TCargo)289);
	addEdge(g, vertBerlin, vertHannover, (TCargo)286);
	addEdge(g, vertBerlin, vertHannover, (TCargo)573);
	addEdge(g, vertBerlin, vertHannover, (TCargo)586);
	addEdge(g, vertHannover, vertMuenchen, (TCargo)572);
	addEdge(g, vertHamburg, vertMainz, (TCargo)521);

	FILE* strmWrite = fopen("C:/Users/Hox/Desktop/graph.dot", "w");
    write(strmWrite, g, DotDrawing());
    fclose(strmWrite);

	::std::cout << g << ::std::endl;

    
    return 0;
}