#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/arg_parse.h>
#include <iostream>
#include <ios>
#include <fstream>
#include <climits>
#include <sstream>
#include <seqan/sequence.h>
#include <seqan/seq_io.h>


// --------------------------------------------------------------------------
// Class AppOptions
// --------------------------------------------------------------------------

// This struct stores the options from the command line.
//
// You might want to rename this to reflect the name of your app.

struct AppOptions
{
    // Verbosity level.  0 -- quiet, 1 -- normal, 2 -- verbose, 3 -- very verbose.
    int verbosity;

    // The first (and only) argument of the program is stored here.
    seqan::CharString inputFile;

    // maximal number of reads to analyse, exit as if the file contained only the first #maxReads 
    __uint64 maxReads;





    // dump a human readable representation into a given output stream
    void out(std::ostream & os);

    AppOptions() :
        verbosity(1),
        maxReads(ULONG_MAX)
    {}
};
