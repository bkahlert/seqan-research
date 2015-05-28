#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/arg_parse.h>
#include <iostream>
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
    seqan::CharString text;
    
    void stream(basic_ostream & os);

    AppOptions() :
        verbosity(1)
    {}
};
