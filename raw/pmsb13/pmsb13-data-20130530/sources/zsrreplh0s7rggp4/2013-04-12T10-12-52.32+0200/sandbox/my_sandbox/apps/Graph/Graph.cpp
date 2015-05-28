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
    
    TVertexDescriptor arrayedges[]={vert1,vert0,vert0,vert4,vert2,vert1,vert4,vert1,vert5,vert1,vert6,vert2,vert3,vert2,vert2,vert3,vert7,vert3,vert5,vert4,vert6,vert5,vert5,vert6,vert7,vert6,vert7,vert7};
    
    addEdges(g,arrayedges,14);
    //FILE* strmWrite = fopen("graph.dot", "w");
    //write(strmWrite, g, DotDrawing());
    //fclose(strmWrite);
    ::std::cout << g << ::std::endl;
    
    return 0;
}