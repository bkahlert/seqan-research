#include <iostream>
#include <seqan/graph_types.h>
#include <seqan/graph_algorithms.h>
using namespace seqan;
int main()
{
    typedef unsigned int TCargo;
    typedef Graph<Undirected<TCargo> > TGraph;
    typedef VertexDescriptor<TGraph>::Type TVertexDescriptor;
    TGraph g;
    typedef String<char> TCityName;
    typedef String<TCityName> TProperties;
    TProperties cityNames;

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
    resizeVertexMap(g, cityNames);
    String<char> shortcuts;
    typedef Iterator<TGraph, VertexIterator>::Type TVertexIterator;
    TVertexIterator itV(g);
    for (int i = 0; !atEnd(itV); goNext(itV),++i)
    {
        assignProperty(cityNames,value(itV),i+'a');
        std::cout << value(itV) << ':' << getProperty(cityNames, value(itV))
                  << std::endl;
    }
    TGraph a;
    TGraph b;
    TGraph c;
    depthFirstSearch(g,a,b,c);
    return 0;
}
