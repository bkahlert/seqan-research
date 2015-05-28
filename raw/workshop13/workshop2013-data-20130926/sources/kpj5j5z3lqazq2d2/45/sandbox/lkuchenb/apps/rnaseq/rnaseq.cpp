#include <iostream>
#include <seqan/store.h>
#include <seqan/arg_parse.h>
#include <seqan/misc/misc_interval_tree.h>
#include <seqan/parallel.h>

using namespace seqan;
// define used types
typedef FragmentStore<> TStore;
typedef Value<TStore::TAnnotationStore>::Type TAnnotation;
typedef TAnnotation::TId TId;
typedef TAnnotation::TId TPos;
typedef IntervalAndCargo<TPos, TId> TInterval;
typedef IntervalTree<TPos, TId>                 TIntervalTree;
typedef Value<TStore::TAlignedReadStore>::Type TAlignedRead;

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
    std::ifstream alignStream(toCString(options.alignmentFileName));
    std::ifstream annoStream(toCString(options.annotationFileName));
	    
    if (!alignStream.good()) {
	std::cerr << "Could not open the specified SAM file " << options.annotationFileName << std::endl;
	return false;
    }
    read(alignStream, store, Sam());
    std::cout << "read " << length(store.alignedReadStore) << " aligned reads" << std::endl;
    
    if (!annoStream.good()) {
	std::cerr << "Could not open the specified GFF file " << options.annotationFileName << std::endl;
	return false;
    }
    read(annoStream, store, Gtf());
    std::cout << "read " << length(store.annotationStore) << " annotations" << std::endl;

    return true;
}

//
// 3. Extract intervals from gene annotations (grouped by contigId)
//
void extractGeneIntervals(String<String<TInterval> > & intervals, TStore const & store)
{
    Iterator<TStore const, AnnotationTree<> >::Type it = begin(store,  AnnotationTree<>());
    if  (!goDown(it)) {
	return;
    }
    resize(intervals, length(store.contigStore));

    do {
	TPos beginPos = getAnnotation(it).beginPos;
	TPos endPos   = getAnnotation(it).endPos;
	TId  contigId = getAnnotation(it).contigId;
	if (endPos < beginPos)
	    std::swap(beginPos, endPos);
	appendValue(intervals[contigId], TInterval(beginPos, endPos, value(it))); 
    } while (goRight(it));

}

//
// 4. Construct interval trees
//
void constructIntervalTrees(String<TIntervalTree> & intervalTrees, String<String<TInterval> > const & intervals)
{
    resize(intervalTrees, length(intervals));
    for (unsigned i=0; i<length(intervals); ++i) {
	createIntervalTree(intervalTrees[i], intervals[i] );
    }
}

//
// 5. Count reads per gene
//
void countReadsPerGene(String<unsigned> & readsPerGene, String<TIntervalTree> const & intervalTrees, TStore const & store)
{
    typedef Iterator<TStore::TAlignedReadStore const, Rooted>::Type TAlignedReadStoreIterator;

    resize(readsPerGene, length(store.annotationStore), 0);

    for (TAlignedReadStoreIterator it = begin(store.alignedReadStore); !atEnd(it); goNext(it)) {
	TPos beginPos = it->beginPosition;
    }
}

int main(int argc, char const * argv[])
{
    Options options;
    TStore store;
    String<String<TInterval> > intervals;
    String<TIntervalTree> intervalTrees;
    String<unsigned> readsPerGene;
    

    ArgumentParser::ParseResult res = parseOptions(options, argc, argv);
    if (res != ArgumentParser::PARSE_OK)
        return res == ArgumentParser::PARSE_ERROR;
    if (!loadFiles(store, options))
        return 1;

    extractGeneIntervals(intervals, store);
    constructIntervalTrees(intervalTrees, intervals);
    countReadsPerGene(readsPerGene, intervalTrees, store);

    return 0;
}