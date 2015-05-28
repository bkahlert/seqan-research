#include <iostream>
#include <seqan/graph_algorithms.h>
using namespace seqan;

int main() 
{
    typedef Graph<Directed<> > TGraph;
    typedef VertexDescriptor<TGraph>::Type TVertexDescriptor;
    typedef EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
    typedef Size<TGraph>::Type TSize;
	TSize numEdges = 14;
    TVertexDescriptor edges[] = {1,0, 0,4, 2,1, 4,1, 5,1, 6,2, 3,2, 2,3, 7,3, 5,4, 6,5, 5,6, 7,6, 7,7};
    TGraph g;
    addEdges(g, edges, numEdges);
    std::cout << g << std::endl;
}