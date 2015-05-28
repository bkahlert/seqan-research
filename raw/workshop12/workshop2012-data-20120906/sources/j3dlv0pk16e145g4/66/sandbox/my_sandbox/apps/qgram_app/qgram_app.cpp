int main ()
{
    typedef Index<DnaString, IndexQGram< OneGappedShape<4> > > TIndex;
    TIndex index("CATGATTACATA");
    stringToShape(indexShape(index), "1101");
    hash(indexShape(index), "ATA");
    for (unsigned i = 0; i < length(getOccurrences(index, indexShape(index))); ++i)
        std::cout << getOccurrences(index, indexShape(index))[i] << std::endl;
    return 0;
}