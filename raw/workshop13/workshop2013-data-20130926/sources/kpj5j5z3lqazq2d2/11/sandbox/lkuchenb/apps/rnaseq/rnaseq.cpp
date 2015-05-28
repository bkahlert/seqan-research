#include <iostream>
#include <seqan/store.h>
#include <seqan/arg_parse.h>
#include <seqan/misc/misc_interval_tree.h>
#include <seqan/parallel.h>

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
    std::ifstream aligStream(toCString(options.alignmentFileName));
    std::ifstream annoStream(toCString(options.annotationFileName));
	    
    if (!alignSequence.good()) {
	std::cerr << "Could not open the specified SAM file " << options.annotationFileName << std::endl;
	return false;
    }
    read(aligStream, store, Sam());
    
    if (!annoStream.good()) {
	std::cerr << "Could not open the specified GFF file " << options.annotationFileName << std::endl;
	return false;
    }
    read(annoStream, store, Gtf())

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