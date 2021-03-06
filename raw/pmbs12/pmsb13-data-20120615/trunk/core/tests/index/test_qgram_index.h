// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2010, Knut Reinert, FU Berlin
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ==========================================================================
// Author: David Weese <david.weese@fu-berlin.de>
// ==========================================================================

#ifndef TESTS_INDEX_TEST_QGRAM_INDEX_H
#define TESTS_INDEX_TEST_QGRAM_INDEX_H


//////////////////////////////////////////////////////////////////////////////

namespace SEQAN_NAMESPACE_MAIN
{

SEQAN_DEFINE_TEST(testGappedShapes)
{
	String<char> shape_string = "0010011011101";
	Shape<Dna,GenericShape> shape1;
	stringToShape(shape1, shape_string);
	Shape<Dna,GenericShape> shape2 = Shape<Dna,GenericShape>(shape1);

	SEQAN_ASSERT(shape1.weight == shape2.weight);
	SEQAN_ASSERT(shape1.span == shape2.span);
	SEQAN_ASSERT(shape1.diffs == shape2.diffs);
	SEQAN_ASSERT(length(shape1) == length(shape2));
	SEQAN_ASSERT(weight(shape1) == weight(shape2));
/*
	Shape<Dna,GenericShape> shape3 = Shape<Dna,GenericShape>(5, 13);
	shape3[0]=2;
	shape3[1]=2;
	shape3[2]=1;
	shape3[3]=1;
	shape3[4]=2;
	shape3[5]=1;
	shape3[6]=1;
	shape3[7]=2;
	for(int i = 0; i < 8; ++i)
        SEQAN_ASSERT(shape1[i] == shape3[i]);
*/
}


SEQAN_DEFINE_TEST(testUngappedShapes)
{
	Shape<Dna,SimpleShape> shape1;
	resize(shape1, 4);
	Shape<Dna,SimpleShape> shape2 = Shape<Dna,SimpleShape>(shape1);

	SEQAN_ASSERT(shape1.span == shape2.span);
	SEQAN_ASSERT(shape1.leftFactor == shape2.leftFactor);
	SEQAN_ASSERT(length(shape1) == length(shape2));
	SEQAN_ASSERT(weight(shape1) == weight(shape2));
	

	Shape<Dna,SimpleShape> shape3 = Shape<Dna,SimpleShape>(4);
	SEQAN_ASSERT(shape3.leftFactor == 64);
	SEQAN_ASSERT(shape1.leftFactor == shape3.leftFactor);


}

template <typename TIndex>
void testStepSize()
{
    TIndex index("CATGATTACATA");
    setStepSize(index,2);
    hash(indexShape(index), "CAT");
    String<typename Position<DnaString>::Type> occs;
    occs = getOccurrences(index, indexShape(index));
    SEQAN_ASSERT_EQ(length(occs), 2u);
    SEQAN_ASSERT_EQ(occs[0], 0u);
    SEQAN_ASSERT_EQ(occs[1], 8u);
}

SEQAN_DEFINE_TEST(testStepSize)
{
    typedef Index<DnaString, IndexQGram< UngappedShape<3> > > TIndex1;
    typedef Index<DnaString, IndexQGram< UngappedShape<3>, OpenAddressing > > TIndex2;

    testStepSize<TIndex1>();
    testStepSize<TIndex2>();
}

/*
void testQGramIndexSchnell()
{
	clock_t start, finish;
	double duration;

	String<Dna> text;
	fstream strm_t;
	strm_t.open(TEST_PATH "fasta.txt", ios_base::in);
	read(strm_t, text, Fasta());
	String<Dna> next;
	resize(next,length(text));
	for(int i = 0; i < 1; ++i)							// datei zu ende?
	{
		arrayCopyForward(begin(text),end(text),begin(next));
		append(text, next);
	}
	strm_t.close();

	String<char> shape_string = "1000100100001110011";
	int q = length(shape_string);
	Shape<Dna,GenericShape> shape;
	stringToShape(shape, shape_string);

	typedef Position<String<Dna> >::Type TPosition;
	String<TPosition> index;
	resize(index, length(text) - q + 2);	
	
	String<TPosition> pos;	
	int pos_size = _intPow((unsigned)ValueSize<Dna>::VALUE, weight(shape));
	pos_size += 1;	
	resize(pos, pos_size);

	start = clock();
	Nothing nothing;
	createQGramIndex(index, pos, nothing, text, shape, 1);
	finish = clock();
	duration = (double)(finish - start) / CLOCKS_PER_SEC;
	//std::cout << "\nQGramIndex bauen dauert: " << duration << " Sekunden.\n\n";
	
	
}
*/
/*
void testGappedQGramIndex()
{
	String<Dna> text = "CTGAACCCTAAACCCT";
	String<char> shape_string = "101";
	int q = length(shape_string);
	Shape<Dna,GenericShape> shape;
	stringToShape(shape, shape_string);

	typedef Position<String<Dna> >::Type TPosition;
	String<TPosition> index;
	resize(index, length(text) - q + 2);
	
	String<TPosition> pos;
    int pos_size = _intPow((unsigned)ValueSize<Dna>::VALUE, weight(shape));
	pos_size += 1;	
	resize(pos, pos_size);

	Nothing nothing;
	createQGramIndex(index, pos, nothing, text, shape, 1);
	
	SEQAN_ASSERT(pos[0] == 0);
	SEQAN_ASSERT(pos[1] == 1);
	SEQAN_ASSERT(pos[2] == 5);
	SEQAN_ASSERT(pos[3] == 5);
	SEQAN_ASSERT(pos[4] == 5);
	SEQAN_ASSERT(pos[5] == 6);
	SEQAN_ASSERT(pos[6] == 8);
	SEQAN_ASSERT(pos[7] == 9);
	SEQAN_ASSERT(pos[8] == 11);
	SEQAN_ASSERT(pos[9] == 12);
	SEQAN_ASSERT(pos[10] == 12);
	SEQAN_ASSERT(pos[11] == 12);
	SEQAN_ASSERT(pos[12] == 12);
	SEQAN_ASSERT(pos[13] == 14);
	SEQAN_ASSERT(pos[14] == 14);
	SEQAN_ASSERT(pos[15] == 14);
	SEQAN_ASSERT(pos[16] == 14);

	SEQAN_ASSERT(index[0] == 9);
	SEQAN_ASSERT(index[1] == 11);
	SEQAN_ASSERT(index[2] == 10);
	SEQAN_ASSERT(index[3] == 4);
	SEQAN_ASSERT(index[4] == 3);
	SEQAN_ASSERT(index[5] == 7);
	SEQAN_ASSERT(index[6] == 12);
	SEQAN_ASSERT(index[7] == 5);
	SEQAN_ASSERT(index[8] == 0);
	SEQAN_ASSERT(index[9] == 13);
	SEQAN_ASSERT(index[10] == 6);
	SEQAN_ASSERT(index[11] == 2);
	SEQAN_ASSERT(index[12] == 8);
	SEQAN_ASSERT(index[13] == 1);
	
}
*/
SEQAN_DEFINE_TEST(testUngappedQGramIndex)
{
	String<Dna> text = "CTGAACCCTAAACCCT";
	int q = 2;
	Shape<Dna,SimpleShape> shape;
	resize(shape, q);

	typedef Position<String<Dna> >::Type TPosition;
	String<TPosition> index;
	resize(index, length(text) - q + 1);
	
	String<TPosition> pos;
    int pos_size = _intPow((unsigned)ValueSize<Dna>::VALUE, q);
	pos_size += 1;	
	resize(pos, pos_size);

	Nothing nothing;
	createQGramIndex(index, pos, nothing, text, shape, 1);
	
	
	SEQAN_ASSERT(pos[0] == 0);
	SEQAN_ASSERT(pos[1] == 3);
	SEQAN_ASSERT(pos[2] == 5);
	SEQAN_ASSERT(pos[3] == 5);
	SEQAN_ASSERT(pos[4] == 5);
	SEQAN_ASSERT(pos[5] == 5);
	SEQAN_ASSERT(pos[6] == 9);
	SEQAN_ASSERT(pos[7] == 9);
	SEQAN_ASSERT(pos[8] == 12);
	SEQAN_ASSERT(pos[9] == 13);
	SEQAN_ASSERT(pos[10] == 13);
	SEQAN_ASSERT(pos[11] == 13);
	SEQAN_ASSERT(pos[12] == 13);
	SEQAN_ASSERT(pos[13] == 14);
	SEQAN_ASSERT(pos[14] == 14);
	SEQAN_ASSERT(pos[15] == 15);

	SEQAN_ASSERT(index[0] == 3);
	SEQAN_ASSERT(index[1] == 9);
	SEQAN_ASSERT(index[2] == 10);
	SEQAN_ASSERT(index[3] == 4);
	SEQAN_ASSERT(index[4] == 11);
	SEQAN_ASSERT(index[5] == 5);
	SEQAN_ASSERT(index[6] == 6);
	SEQAN_ASSERT(index[7] == 12);
	SEQAN_ASSERT(index[8] == 13);
	SEQAN_ASSERT(index[9] == 0);
	SEQAN_ASSERT(index[10] == 7);
	SEQAN_ASSERT(index[11] == 14);
	SEQAN_ASSERT(index[12] == 2);
	SEQAN_ASSERT(index[13] == 8);
	SEQAN_ASSERT(index[14] == 1);
}

//////////////////////////////////////////////////////////////////////////////

SEQAN_DEFINE_TEST(testQGramFind)
{
	typedef Index<String<char>, IndexQGram<UngappedShape<2> > > TQGramIndex;
	TQGramIndex idx("to be or not to be");
	Finder<TQGramIndex> finder(idx);

	SEQAN_ASSERT(find(finder, "be"));
	SEQAN_ASSERT(position(finder) == 3);
	SEQAN_ASSERT(find(finder, "be"));
	SEQAN_ASSERT(position(finder) == 16);
	SEQAN_ASSERT(!find(finder, "be"));
/*
	while (find(finder, "be"))
	{
		std::cout << position(finder) << "\n";
	}
*/
}

//////////////////////////////////////////////////////////////////////////////


} //namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
