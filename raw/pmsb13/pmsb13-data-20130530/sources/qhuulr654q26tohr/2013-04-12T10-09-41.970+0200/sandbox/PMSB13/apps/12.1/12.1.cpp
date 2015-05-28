#include <iostream>
#include <seqan/graph_types.h>
#include <seqan/graph_algorithms.h>
using namespace seqan;
using namespace std;

int main ()
{
    typedef unsigned int TCargo;
    typedef Graph<Directed<TCargo> > TGraph;
    typedef VertexDescriptor<TGraph>::Type TVertexDescriptor;
    
    TGraph g;
    
    TVertexDescriptor v0 = addVertex(g);
    TVertexDescriptor v1 = addVertex(g);
    TVertexDescriptor v2 = addVertex(g);
    TVertexDescriptor v3 = addVertex(g);
    TVertexDescriptor v4 = addVertex(g);
    TVertexDescriptor v5 = addVertex(g);
    TVertexDescriptor v6 = addVertex(g);
    TVertexDescriptor v7 = addVertex(g);
    
    //String<TVertexDescriptor> edges;
    
    TVertexDescriptor edges[] = {v1,v0,v0,v4,v2,v1,v4,v1,v5,v1,v6,v2,v3,v2,v2,v3,v7,v3,v5,v4,v6,v5,v5,v6,v7,v6,v7,v7};
    TVertexDescriptor edges2[] = {v1,v0,v0,v4}
    
    addEdges(g,edges2,4);
    
    //="(1,0), (0,4), (2,1), (4,1), (5,1), (6,2), (3,2), (2,3), (7,3), (5,4), (6,5), (5,6), (7,6), (7,7)";
    
    
    FILE* strmWrite = fopen("graph.dot", "w");
    write(strmWrite, g, DotDrawing());
    fclose(strmWrite);
}