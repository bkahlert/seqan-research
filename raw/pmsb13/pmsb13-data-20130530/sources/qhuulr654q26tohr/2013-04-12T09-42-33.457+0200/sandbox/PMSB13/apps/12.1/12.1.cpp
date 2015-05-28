#include <iostream>
#include <seqan/graph_types.h>
#include <seqan/graph_algorithms.h>
using namespace seqan;
using namespace std;

int main ()
{
    typedef unsigned int TCargo;
    typedef Graph<Undirected<TCargo> > TGraph;
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
    
    TVertexDescriptor edges[] = {v1,v0,v0,v4,v2,v1};
    
    cout << edges.length() << endl;
    addEdges(g,edges,edges.length());
    
    //="(1,0), (0,4), (2,1), (4,1), (5,1), (6,2), (3,2), (2,3), (7,3), (5,4), (6,5), (5,6), (7,6), (7,7)";
    
    
    FILE* strmWrite = fopen("graph.dot", "w");
    write(strmWrite, g, DotDrawing());
    fclose(strmWrite);
}