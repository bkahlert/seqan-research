#include <iostream>
#include <seqan/graph_algorithms.h>
using namespace seqan;

int main() 
{
    typedef Graph<Undirected<> > TGraph;
    typedef VertexDescriptor<TGraph>::Type TVertexDescriptor;
    typedef EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
    typedef Size<TGraph>::Type TSize;
	typedef String<char> TCityName;
    typedef String<TCityName> TProperties;

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
 
	TProperties cityNames;
    resizeVertexMap(g, cityNames);
	assignProperty(cityNames, vertBerlin, "Berlin");
    assignProperty(cityNames, vertHamburg, "Hamburg");
    assignProperty(cityNames, vertMuenchen, "Munich");
    assignProperty(cityNames, vertMainz, "Mainz");
    assignProperty(cityNames, vertHannover, "Hannover");

	typedef Iterator<TGraph, VertexIterator>::Type TVertexIterator;
    TVertexIterator itV(g);
    for(;!atEnd(itV);goNext(itV)) {
        std::cout << value(itV) << ':' << getProperty(cityNames, value(itV)) << std::endl;
    }

	std::cout << g << std::endl;
}