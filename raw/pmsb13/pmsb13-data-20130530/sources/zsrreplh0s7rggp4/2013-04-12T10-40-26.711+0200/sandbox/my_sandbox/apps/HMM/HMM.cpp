#include <iostream>
#include <fstream>
#include <seqan/basic.h>
#include <seqan/graph_algorithms.h>

typedef LogProb<> TProbability;
typedef Dna TAlphabet;    
typedef Size<TAlphabet>::Type TSize;
typedef Graph<Hmm<TAlphabet, TProbability, Default> > THmm;
typedef VertexDescriptor<THmm>::Type TVertexDescriptor;
typedef EdgeDescriptor<THmm>::Type TEdgeDescriptor;  

Dna dnaA = Dna('A');
Dna dnaC = Dna('C');
Dna dnaG = Dna('G');
Dna dnaT = Dna('T');

THmm hmm;

TVertexDescriptor begState = addVertex(hmm);
assignBeginState(hmm, begState);
