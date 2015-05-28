// BLASTX IMPLEMENTIERUNG VON ANNKATRIN BRESSIN UND MARJAN FAIZI
// SOFTWAREPROJEKT VOM 2.4. - 29.5.2012
// VERSION VOM 

int main(){

typedef String<char> TSequence;                 // sequence type
typedef Align<TSequence, ArrayGaps> TAlign;     // align type

TSequence seq1 = "GARFIELDTHECAT";
TSequence seq2 = "GARFIELDTHEBIGCAT";

TAlign align;
resize(rows(align),2);
assignSource(row(align,0),seq1);
assignSource(row(align,1),seq2);
int score = globalAlignment(align,Score<int,Simple>(1,-1,-1),MyerBitVector());

cout << "Score: " << score <<endl;
cout << align <<endl;

    return 0;
}
