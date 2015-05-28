// // ==========================================================================
// //                                workshop_app
// // ==========================================================================
// // Copyright (c) 2006-2012, Knut Reinert, FU Berlin
// // All rights reserved.
// //
// // Redistribution and use in source and binary forms, with or without
// // modification, are permitted provided that the following conditions are met:
// //
// //     * Redistributions of source code must retain the above copyright
// //       notice, this list of conditions and the following disclaimer.
// //     * Redistributions in binary form must reproduce the above copyright
// //       notice, this list of conditions and the following disclaimer in the
// //       documentation and/or other materials provided with the distribution.
// //     * Neither the name of Knut Reinert or the FU Berlin nor the names of
// //       its contributors may be used to endorse or promote products derived
// //       from this software without specific prior written permission.
// //
// // THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// // AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// // IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// // ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
// // FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// // DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// // SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// // CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// // LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// // OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// // DAMAGE.
// //
// // ==========================================================================
// // Author: Your Name <your.email@example.net>
// // ==========================================================================
// 
// #include <iostream>
// 
// #include <seqan/basic.h>
// #include <seqan/sequence.h>
// #include <seqan/seq_io.h>
// 
// #include <seqan/arg_parse.h>
// 
// // ==========================================================================
// // Classes
// // ==========================================================================
// 
// // --------------------------------------------------------------------------
// // Class AppOptions
// // --------------------------------------------------------------------------
// 
// // This struct stores the options from the command line.
// //
// // You might want to rename this to reflect the name of your app.
// 
// struct AppOptions
// {
//     // Verbosity level.  0 -- quiet, 1 -- normal, 2 -- verbose, 3 -- very verbose.
//     int verbosity;
// 
//     // The first (and only) argument of the program is stored here.
//     seqan::CharString text;
// 
//     AppOptions() :
//         verbosity(1)
//     {}
// };
// 
// // ==========================================================================
// // Functions
// // ==========================================================================
// 
// // --------------------------------------------------------------------------
// // Function parseCommandLine()
// // --------------------------------------------------------------------------
// 
// seqan::ArgumentParser::ParseResult
// parseCommandLine(AppOptions & options, int argc, char const ** argv)
// {
//     // Setup ArgumentParser.
//     seqan::ArgumentParser parser("workshop_app");
//     // Set short description, version, and date.
//     setShortDescription(parser, "Put a Short Description Here");
//     setVersion(parser, "0.1");
//     setDate(parser, "July 2012");
// 
//     // Define usage line and long description.
//     addUsageLine(parser, "[\\fIOPTIONS\\fP] \"\\fITEXT\\fP\"");
//     addDescription(parser, "This is the application skelleton and you should modify this string.");
// 
//     // We require one argument.
//     addArgument(parser, seqan::ArgParseArgument(seqan::ArgParseArgument::STRING, "TEXT"));
// 
//     addOption(parser, seqan::ArgParseOption("q", "quiet", "Set verbosity to a minimum."));
//     addOption(parser, seqan::ArgParseOption("v", "verbose", "Enable verbose output."));
//     addOption(parser, seqan::ArgParseOption("vv", "very-verbose", "Enable very verbose output."));
// 
//     // Add Examples Section.
//     addTextSection(parser, "Examples");
//     addListItem(parser, "\\fBworkshop_app\\fP \\fB-v\\fP \\fItext\\fP",
//                 "Call with \\fITEXT\\fP set to \"text\" with verbose output.");
// 
//     // Parse command line.
//     seqan::ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);
// 
//     // Only extract  options if the program will continue after parseCommandLine()
//     if (res != seqan::ArgumentParser::PARSE_OK)
//         return res;
// 
//     // Extract option values.
//     if (isSet(parser, "quiet"))
//         options.verbosity = 0;
//     if (isSet(parser, "verbose"))
//         options.verbosity = 2;
//     if (isSet(parser, "very-verbose"))
//         options.verbosity = 3;
//     seqan::getArgumentValue(options.text, parser, 0);
// 
//     return seqan::ArgumentParser::PARSE_OK;
// }
// 
// // --------------------------------------------------------------------------
// // Function main()
// // --------------------------------------------------------------------------
// 
// // Program entry point.
// 
// int main(int argc, char const ** argv)
// {
//     // Parse the command line.
//     seqan::ArgumentParser parser;
//     AppOptions options;
//     seqan::ArgumentParser::ParseResult res = parseCommandLine(options, argc, argv);
// 
//     // If there was an error parsing or built-in argument parser functionality
//     // was triggered then we exit the program.  The return code is 1 if there
//     // were errors and 0 if there were none.
//     if (res != seqan::ArgumentParser::PARSE_OK)
//         return res == seqan::ArgumentParser::PARSE_ERROR;
// 
//     std::cout << "EXAMPLE PROGRAM\n"
//               << "===============\n\n";
// 	
//     // Print the command line arguments back to the user.
//     if (options.verbosity > 0)
//     {
//         std::cout << "__OPTIONS____________________________________________________________________\n"
//                   << '\n'
//                   << "VERBOSITY\t" << options.verbosity << '\n'
//                   << "TEXT     \t" << options.text << "\n\n";
//     }
// 
// 	seqan::CharString id;
//     seqan::Dna5String seq;
//     seqan::SequenceStream seqStream("example.fa.gz");
//     if (!isGood(seqStream))
//     {
//         std::cerr << "ERROR: Could not open the file.\n";
//         return 1;
//     }
//     
// 	while(! atEnd(seqStream)){
// 		if (readRecord(id, seq, seqStream) != 0){
// 	        std::cerr << "ERROR: Could not read from example.fa!\n";
// 	        return 1;
// 	    }
// 		
// 		std::cout << id << '\t' << seq << '\n';
// 	}
// 	
// 	
// 	if (argc < 2)
// 	    {
// 	        std::cerr << "USAGE: basic_seq_io_example FILENAME\n";
// 	        return 1;
// 	    }
// 	    seqan::SequenceStream seqStream2(argv[1], seqan::SequenceStream::WRITE);
// 	    if (!isGood(seqStream2))
// 	    {
// 	        std::cerr << "ERROR: Could not open the file.\n";
// 	        return 1;
// 	    }
// 	    seqan::CharString id2 = "seq1";
// 	    seqan::Dna5String seq2 = "CGAATTTACCCAGGGATAAGT";
// 	    if (writeRecord(seqStream2, id2, seq2) != 0)
// 	    {
// 	        std::cerr << "ERROR: Could not write to file!\n";
// 	        return 1;
// 	    }
// 	
// 		
//     return 0;
// }



#include <iostream>
#include <seqan/store.h>
#include <seqan/arg_parse.h>
#include <seqan/misc/misc_interval_tree.h>
#include <seqan/parallel.h>
#include <seqan/file.h>

using namespace seqan;
// define used types
typedef FragmentStore<> TStore;
// define options
struct Options
{
    std::string annotationFileName;
    std::string alignmentFileName;
};
//
// 1. Parse command line and fill Options object
//
ArgumentParser::ParseResult parseOptions(Options & options, int argc, char const * argv[])
{
    ArgumentParser parser("gene_quant");
    setShortDescription(parser, "A simple gene quantification tool");
    setVersion(parser, "1.0");
    setDate(parser, "Sep 2012");
    addArgument(parser, ArgParseArgument(ArgParseArgument::INPUTFILE));
    addArgument(parser, ArgParseArgument(ArgParseArgument::INPUTFILE));
    addUsageLine(parser, "[\\fIOPTIONS\\fP] <\\fIANNOTATION FILE\\fP> <\\fIREAD ALIGNMENT FILE\\fP>");
    // Parse command line
    ArgumentParser::ParseResult res = parse(parser, argc, argv);
    if (res == ArgumentParser::PARSE_OK)
    {
        // Extract option values
        getArgumentValue(options.annotationFileName, parser, 0);
        getArgumentValue(options.alignmentFileName, parser, 1);
    }
    return res;
}
//
// 2. Load annotations and alignments from files
//
bool loadFiles(TStore & store, Options const & options)
{
    // INSERT YOUR CODE HERE ...
    //
	std::ifstream is;
	
	is.open(options.alignmentFileName.c_str());
	
	
    return true;
}
int main(int argc, char const * argv[])
{
    Options options;
    TStore store;
    ArgumentParser::ParseResult res = parseOptions(options, argc, argv);
    if (res != ArgumentParser::PARSE_OK)
        return res == ArgumentParser::PARSE_ERROR;
    if (!loadFiles(store, options))
        return 1;
    return 0;
}