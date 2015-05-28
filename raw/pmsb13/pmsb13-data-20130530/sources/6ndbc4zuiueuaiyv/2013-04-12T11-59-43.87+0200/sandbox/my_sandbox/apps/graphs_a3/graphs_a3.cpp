#include <iostream>
#include <seqan/graph_types.h>
#include <seqan/graph_algorithms.h>
#include <seqan/basic/basic_math.h>
using namespace seqan;
using namespace std;

int main()
{
    typedef String<char> TCityName;
    typedef String<TCityName> TProperties;
    TGraph g;
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

    return 0;
    
}

