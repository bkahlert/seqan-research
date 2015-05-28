#include <seqan/basic.h>
#include <seqan/sequence.h>

#include <seqan/arg_parse.h>

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
    seqan::ArgumentParser parser("io_fasta");
    // Set short description, version, and date.
    setShortDescription(parser, "reads fastA file");
    setVersion(parser, "0.1");
    setDate(parser, "July 2012");

    // Define usage line and long description.
    addUsageLine(parser, "[\\fIOPTIONS\\fP] <\\fIfastA-file\\fP>");
    addDescription(parser, "IO tutorial example");

    // We require one argument.
    addArgument(parser, seqan::ArgParseArgument(seqan::ArgParseArgument::STRING, "TEXT"));

    addOption(parser, seqan::ArgParseOption("q", "quiet", "Set verbosity to a minimum."));
    addOption(parser, seqan::ArgParseOption("v", "verbose", "Enable verbose output."));
    addOption(parser, seqan::ArgParseOption("vv", "very-verbose", "Enable very verbose output."));

    // Add Examples Section.
    addTextSection(parser, "Examples");
    addListItem(parser, "\\fBio_fasta\\fP \\fB-v\\fP \\fItest.fa\\fP",
                "read \\fItest.fa\\fP set to \"test.fa\" with verbose output.");

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

    std::cout << "EXAMPLE PROGRAM\n"
              << "===============\n\n";
    
    // Print the command line arguments back to the user.
    if (options.verbosity > 1)
    {
        std::cout << "__OPTIONS____________________________________________________________________\n"
                  << '\n'
                  << "VERBOSITY\t" << options.verbosity << '\n'
                  << "file.fasta\t" << options.file << "\n\n";
    }
    seqan::CharString id;
    seqan::Dna5String seq;
    seqan::SequenceStream seqStream(options.file);
    readRecord(id, seq, seqStream);
    std::cout << id << '\t' << seq << '\n';
    return 0;

    return 0;
}