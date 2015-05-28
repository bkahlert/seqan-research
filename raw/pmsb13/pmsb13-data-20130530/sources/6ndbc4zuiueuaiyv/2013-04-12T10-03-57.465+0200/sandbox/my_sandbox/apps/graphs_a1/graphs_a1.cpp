#include <iostream>
#include <seqan/graph_types.h>
#include <seqan/graph_algorithms.h>
using namespace seqan;
using namespace std;


int main ()
{
    typedef unsigned int TCargo;
    typedef unsigned int TSpec;
    typedef Graph<Directed<TCargo, TSpec> > TGraph;
    typedef VertexDescriptor<TGraph>::Type TVertexDescriptor;

    TGraph g;
    TVertexDescriptor edges[] = {1,0, 0,4, 2,1, 4,1, 5,1, 6,2, 3,2, 2,3, 7,3, 5,4, 6,5, 5,6, 7,6, 7,7};
    
    cout << sizeof(edges)/sizeof(TVertexDescriptor) << endl;
    addEdges(g, edges, 10);


    FILE* strmWrite = fopen("/home/stefan/Dokumente/subversion/seqan-trunk/sandbox/my_sandbox/apps/graphs_a1/graph.dot", "w");
    write(strmWrite, g, DotDrawing());
    fclose(strmWrite);

    ::std::cout << g << ::std::endl;
    
    return 0;
    
}
