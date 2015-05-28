#include <iostream>
#include <seqan/store.h>
#include <seqan/arg_parse.h>
#include <seqan/misc/misc_interval_tree.h>
#include <seqan/parallel.h>
#include <seqan/sequence.h>
#include <seqan/seq_io.h>
#include <fstream>

using namespace seqan;

typedef FragmentStore<> TStore;

struct Options
{
    std::string annotationFileName;
    std::string alignmentFileName;
};

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
	/*
    std::ifstream alignmentFile(options.alignmentFileName.c_str());
    if (!alignmentFile.good())
    {
        std::cerr << "Couldn't open alignment file " << options.alignmentFileName << std::endl;
        return false;
    }
    std::cerr << "Loading read alignments ..... " << std::flush;
    read(alignmentFile, store, Sam());
    std::cerr << "[" << length(store.alignedReadStore) << "]" << std::endl;
    // load annotations
    std::ifstream annotationFile(options.annotationFileName.c_str());
    if (!annotationFile.good())
    {
        std::cerr << "Couldn't open annotation file" << options.annotationFileName << std::endl;
        return false;
    }
    std::cerr << "Loading genome annotation ... " << std::flush;
    read(annotationFile, store, Gtf());
    std::cerr << "[" << length(store.annotationStore) << "]" << std::endl;
    return true;*/
    
	std::ifstream file("/Informatik/Development/assignment_annotations.gtf", std::ios_base::in | std::ios_base::binary);
	
	read(file, store, Gtf());
    // Create AnnotationTree iterator
	
    std::ifstream annotationFile(options.annotationFileName.c_str());

    Iterator<FragmentStore<>, AnnotationTree<> >::Type it;
    it = begin(store, AnnotationTree<>());
    // Move iterator one node down
    goDown(it); 

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