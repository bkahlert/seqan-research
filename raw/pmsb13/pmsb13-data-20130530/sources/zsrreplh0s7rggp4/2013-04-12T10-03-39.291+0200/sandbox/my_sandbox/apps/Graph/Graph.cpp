#include <iostream>
#include <seqan/graph_types.h>
#include <seqan/graph_algorithms.h>
using namespace seqan;

int main ()
{
    typedef unsigned int TCargo;
    typedef Graph<Directed<TCargo> > TGraph;
    typedef VertexDescriptor<TGraph>::Type TVertexDescriptor;
    
    TGraph g;
    
    TVertexDescriptor vert0 = addVertex(g);
    TVertexDescriptor vert1 = addVertex(g);
    TVertexDescriptor vert2 = addVertex(g);
    TVertexDescriptor vert3 = addVertex(g);
    TVertexDescriptor vert4 = addVertex(g);
    TVertexDescriptor vert5 = addVertex(g);
    TVertexDescriptor vert6 = addVertex(g);
    TVertexDescriptor vert7 = addVertex(g);
    
    /*addEdge(g, vertBerlin, vertHamburg, 289);
    addEdge(g, vertBerlin, vertHannover, 286);
    addEdge(g, vertBerlin, vertMainz, 573);
    addEdge(g, vertBerlin, vertMuenchen, 586);
    addEdge(g, vertHannover, vertMuenchen, 572);
    addEdge(g, vertHamburg, vertMainz, 521);*/
    
    int arrayedges[]={1,0,0,4,2,1,4,1,5,1,6,2,3,2,2,3,7,3,5,4,6,5,5,6,7,6,7,7};
    
    addEdges(g,arrayedges,28);
    //FILE* strmWrite = fopen("graph.dot", "w");
    //write(strmWrite, g, DotDrawing());
    //fclose(strmWrite);
    ::std::cout << g << ::std::endl;
    
    return 0;
}