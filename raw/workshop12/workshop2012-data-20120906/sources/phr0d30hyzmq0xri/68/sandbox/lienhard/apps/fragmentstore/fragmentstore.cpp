#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/arg_parse.h>
#include <fstream>
#include <iostream>
#include <seqan/file.h>
#include <seqan/store.h>
using namespace seqan;
// ==========================================================================
// Classes
// ==========================================================================

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
    seqan::CharString file;

    AppOptions() :
        verbosity(1)
    {}
};

// ==========================================================================
// Functions
// ==========================================================================

// --------------------------------------------------------------------------
// Function parseCommandLine()
// --------------------------------------------------------------------------

seqan::ArgumentParser::ParseResult
parseCommandLine(AppOptions & options, int argc, char const ** argv)
{
    // Setup ArgumentParser.
    seqan::ArgumentParser parser("fragmentstore");
    // Set short description, version, and date.
    setShortDescription(parser, "Fragment Store Example");
    setVersion(parser, "0.1");
    setDate(parser, "July 2012");

    // Define usage line and long description.
    addUsageLine(parser, "[\\fIOPTIONS\\fP] \"\\fIFILE\\fP\"");
    addDescription(parser, "This is a fragmentstore example.");

    // We require one argument.
    addArgument(parser, seqan::ArgParseArgument(seqan::ArgParseArgument::STRING, "FILE"));

    addOption(parser, seqan::ArgParseOption("q", "quiet", "Set verbosity to a minimum."));
    addOption(parser, seqan::ArgParseOption("v", "verbose", "Enable verbose output."));
    addOption(parser, seqan::ArgParseOption("vv", "very-verbose", "Enable very verbose output."));

    // Add Examples Section.
    addTextSection(parser, "Examples");
    addListItem(parser, "\\fBfragmentstore\\fP \\fB-v\\fP \\fIfile.txt\\fP",
                "Call with \\fIFILE\\fP set to \"file.txt\" with verbose output.");

    // Parse command line.
    seqan::ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);

    // Only extract  options if the program will continue after parseCommandLine()
    if (res != seqan::ArgumentParser::PARSE_OK)
        return res;

    // Extract option values.
    if (isSet(parser, "quiet"))
        options.verbosity = 0;
    if (isSet(parser, "verbose"))
        options.verbosity = 2;
    if (isSet(parser, "very-verbose"))
        options.verbosity = 3;
    seqan::getArgumentValue(options.file, parser, 0);

    return seqan::ArgumentParser::PARSE_OK;
}

// --------------------------------------------------------------------------
// Function main()
// --------------------------------------------------------------------------

// Program entry point.

int main(int argc, char const ** argv)
{
    // Parse the command line.
    seqan::ArgumentParser parser;
    AppOptions options;
    seqan::ArgumentParser::ParseResult res = parseCommandLine(options, argc, argv);

    // If there was an error parsing or built-in argument parser functionality
    // was triggered then we exit the program.  The return code is 1 if there
    // were errors and 0 if there were none.
    if (res != seqan::ArgumentParser::PARSE_OK)
        return res == seqan::ArgumentParser::PARSE_ERROR;

    std::cout << "Fragment Store EXAMPLE PROGRAM\n"
              << "===============\n\n";
    
    // Print the command line arguments back to the user.

    if (options.verbosity > 1)
    {
        std::cout << "__OPTIONS____________________________________________________________________\n"
                  << '\n'
                  << "VERBOSITY\t" << options.verbosity << '\n'
                  << "FILE     \t" << options.file << "\n\n";
    }

    FragmentStore<> store;
    std::ifstream file(toCString(options.file), std::ios_base::in | std::ios_base::binary);
    if(! file.good() ){
      std::cerr << "Cannot open file " << options.file << "for reading\n";
      return 1;
    }
    read(file, store, Gtf());
    // Create AnnotationTree iterator
    Iterator<FragmentStore<>, AnnotationTree<> >::Type it;
    it = begin(store, AnnotationTree<>());
    // Move iterator one node down
    while(goDown(it));
    unsigned counter=0;
    if(strcmp(toCString(getType(it)),"exon")){++counter;}
    while(goRight(it)){
      if(strcmp(toCString(getType(it)),"exon")){++counter;}
    }
    std::cout << "In mRNA " << getParentName(it) << " gibt es " << counter <<" Exons\n";
    unsigned sum=counter;
    while(! atEnd(it)){
      goNext(it);    
      while(goDown(it));
      counter=0;
      if(strcmp(toCString(getType(it)),"exon"))
        ++counter;
      else stdcout << getType(it) <<"\n";
      while(goRight(it))
        if(strcmp(toCString(getType(it)),"exon"))
	  ++counter;
	else stdcout << getType(it) <<"\n";
      
      std::cout << "In mRNA " << getParentName(it) << " gibt es " << counter <<" Exons\n";
      sum+=counter;
    }
    std::cout << "Insgesammt: " << sum <<" Exons\n";
   


    return 0;
}