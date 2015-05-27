#include <seqan/stream.h>
#include <seqan/index.h>
#include <seqan/align.h>
#include <seqan/seeds2.h>

using namespace seqan;

// set parameters

unsigned const minimalLength = 12;
float const errorRate = 0.05;
unsigned const qGramLength = 8;
CharString const filename = "ministellar.out";

typedef int TScoreValue;
TScoreValue const matchScore = 1;
TScoreValue const mismatchScore = -2;
TScoreValue const gapScore = -2;
TScoreValue const minScore = 9;
TScoreValue const xDrop = 3;

template <typename TSequence>
inline bool loadFasta(StringSet<CharString>& ids, StringSet<TSequence>& sequences, char const* filename) {
    std::ifstream streamIn(filename, std::ios::binary);
    if (not streamIn)
        return false;

    RecordReader<std::ifstream, SinglePass<> > reader(streamIn);
    if (read2(ids, sequences, reader, Fasta()) != 0) {
        std::cerr << "Could not read file " << filename << std::endl;
        return false;
    }

    return true;
}

template<typename TSeed, typename TSeq>
void writeSeed(TSeed & seed, TSeq const & seq0, TSeq const & seq1) {
    std::cout << "Seed in seq1 from position " << getBeginDim0(seed);
    std::cout << " to " << getEndDim0(seed) << ": ";
    std::cout << infix(seq0, getBeginDim0(seed), getEndDim0(seed)) << std::endl;
    std::cout << "Seed in seq2 from position " << getBeginDim1(seed);
    std::cout << " to " << getEndDim1(seed) << ": ";
    std::cout << infix(seq1, getBeginDim1(seed), getEndDim1(seed)) << std::endl;
}

int main(int argc, char const ** argv) {
    if (argc != 3) {
        std::cerr << "ERROR: Invalid argument count!" << std::endl;
        std::cerr << "USAGE: ministallar IN1.fa IN2.fa" << std::endl;
        return 1;
    }

    Score<int, Simple> score(matchScore, mismatchScore, gapScore);

    // read sequences from fasta files
    typedef Dna5String TSequence;
    typedef Value<TSequence>::Type TValue;
    StringSet<CharString> ids1, ids2;
    StringSet<TSequence> seqs1, seqs2;

    if (not loadFasta(ids1, seqs1, argv[1]) or not loadFasta(ids2, seqs2, argv[2]))
        return 1;

    // prepare output file
    std::ofstream outFile(toCString(filename));

    // define finder and pattern
    // TODO: Define types for finder and pattern using SeqAn's find interface from the index
    //       module. Both should be specialized for the approproiate swift filter algorithm.
    //       The pattern will need the definition of a q-gram-index. Use q-grams of length 8
    //       for now, and open addressing.
    // HINT: The tutorial on Pattern Matching describes the find interface. At the end of that
    //       tutorial page, you will find a link to a Swift HowTo page.
    
    typedef Finder<TSequence, Swift<SwiftLocal> > TFinder;
    typedef Index<TSequence, IndexQGram<Shape<TValue, UngappedShape<8> >, OpenAddressing> > TIndex;
    typedef Pattern<TIndex, Swift<SwiftLocal> > TPattern;

    // TODO: Define variables of your pattern and finder types, and initialize them with the
    //       first sequences in the sequence sets seqs1 and seqs2.

    TIndex index(seqs2[0]);
    TPattern pattern(index);
    TFinder finder(seqs1[0]);

    // TODO: Repeat the following steps for all hits of the finder.

    // call the function find to obtain a swift hit
    // TODO: uncomment the following line

    while (find(finder, pattern, errorRate, static_cast<int>(minimalLength))) {
        // extract infixes from hit
        // TODO: uncomment the following three lines

        typedef Infix<TSequence>::Type TInfix;
        TInfix infix1 = infix(finder, seqs1[0]);
        TInfix infix2 = infix(pattern, seqs2[0]);

        // find best local alignment of infixes, and check for minimal score
        // TODO: Define an align object and initialize it with the infixes. Find the best local
        //       alignment of the infixes. Use the above defined scoring scheme. Make sure that
        //       it reaches the minimal score minScore.

        Align<TSequence> align;
        reserve(rows(align), 2, Exact());
        appendValue(rows(align), infix1);
        appendValue(rows(align), infix2);
        if (localAlignment(align, score) < minScore)
            continue;

        // create a seed for the local alignment, and conduct gapped X-drop extension
        // TODO: Define a seed on the original sequences but for the subsequences of the local
        //       alignment. Extend this seed in both directions using gapped X-drop extension.
        //       Check that the extended seed has a length of at least minimalLength.
        // HINT: Have a look at the Seed-and-Extend tutorial.

        typedef Seed<Simple> TSeed;
        TSeed seed(
            clippedBeginPosition(row(align, 0)) + beginPosition(infix1),
            clippedBeginPosition(row(align, 1)) + beginPosition(infix2),
            clippedEndPosition(row(align, 0)) + beginPosition(infix1),
            clippedEndPosition(row(align, 1)) + beginPosition(infix2));
        std::cout << "Before extension: " << std::endl;
        writeSeed(seed, seqs1[0], seqs2[0]);
        extendSeed(seed, seqs2[0], seqs1[0], EXTEND_BOTH, score, xDrop, GappedXDrop());
        std::cout << "After extension: " << std::endl;
        writeSeed(seed, seqs1[0], seqs2[0]);

        // find best global alignment of extended seed
        // TODO: Compute the best global alignment of extended seed and its score.
        // HINT: Create an align object on infixes of the sequences.

        Align<TSequence> extendedAlign;
        reserve(rows(extendedAlign), 2, Exact());
        appendValue(rows(extendedAlign), infix(seqs1[0], getBeginDim0(seed), getEndDim0(seed)));
        appendValue(rows(extendedAlign), infix(seqs2[0], getBeginDim1(seed), getEndDim1(seed)));
        globalAlignment(extendedAlign, score);

        // TODO: Output the alignment as a match to the output file.

        outFile << clippedBeginPosition(row(extendedAlign, 0)) + beginPosition(row(extendedAlign, 0))
                << "\t"
                << clippedEndPosition(row(extendedAlign, 0)) + beginPosition(row(extendedAlign, 0)) << std::endl;
        //outFile << getBeginDim0(seed) << "\t" << getEndDim0(seed) << std::endl;
        //outFile << extendedAlign << std::endl;
    } // while(find(...))
}