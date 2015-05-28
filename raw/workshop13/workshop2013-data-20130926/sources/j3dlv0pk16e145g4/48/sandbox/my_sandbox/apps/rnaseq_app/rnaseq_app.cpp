#include <iostream>
#include <seqan/store.h>
#include <seqan/arg_parse.h>
#include <seqan/misc/misc_interval_tree.h>
#include <seqan/parallel.h>

using namespace seqan;
// define used types
typedef FragmentStore<> TStore;
typedef Value<TStore::TAnnotationStore>::Type   TAnnotation;
typedef TAnnotation::TId                        TId;
typedef TAnnotation::TId                        TPos;
typedef IntervalAndCargo<TPos, TId>             TInterval;
typedef IntervalTree<TPos, TId>                 TIntervalTree;

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
    std::ifstream alnFile(toCString(options.alignmentFileName), std::ifstream::in);
    if (!alnFile.good())
    {
        std::cerr << "ERROR: Could not open alignment file." << std::endl;
        return false;
    }
    read(alnFile, store, Sam());
    alnFile.close();
    std::ifstream annotationFile(toCString(options.annotationFileName), std::ifstream::in);
    if (!annotationFile.good())
    {
        std::cerr << "ERROR: Could not open annotation file." << std::endl;
        return false;
    }
    read(annotationFile, store, Gtf());
    annotationFile.close();
    return true;
}

//
// 3. Extract intervals from gene annotations (grouped by contigId)
//
void extractGeneIntervals(String<String<TInterval> > & intervals, TStore const & store)
{
    // extract intervals from gene annotations (grouped by contigId)
    resize(intervals, length(store.contigStore));
    Iterator<TStore const, AnnotationTree<> >::Type it = begin(store, AnnotationTree<>());
    if (!goDown(it))
        return;
    do
    {
        SEQAN_ASSERT_EQ(getType(it), "gene");
        TPos beginPos = getAnnotation(it).beginPos;
        TPos endPos = getAnnotation(it).endPos;
        TId contigId = getAnnotation(it).contigId;
        if (beginPos > endPos)
            std::swap(beginPos, endPos);
        // insert forward-strand interval of the gene and its annotation id
        appendValue(intervals[contigId], TInterval(beginPos, endPos, value(it)));
    } while (goRight(it));
}

//
// 4. Construct interval trees
//
void constructIntervalTrees(String<TIntervalTree> & intervalTrees, String<String<TInterval> > const & intervals)
{
    Iterator<String<String<TInterval> > >::Type iit = begin(intervals);
    do
    {
        appendValue(intervalTrees, TIntervalTree(getValue(iit)));
    } while (iit.goNext())
}


int main(int argc, char const * argv[])
{
    Options options;
    TStore store;
    String<String<TInterval> > intervals;
    String<TIntervalTree> intervalTrees;
    
    ArgumentParser::ParseResult res = parseOptions(options, argc, argv);
    if (res != ArgumentParser::PARSE_OK)
        return res == ArgumentParser::PARSE_ERROR;
    if (!loadFiles(store, options))
        return 1;
    extractGeneIntervals(intervals, store);
    constructIntervalTrees(intervalTrees, intervals);
    return 0;
}