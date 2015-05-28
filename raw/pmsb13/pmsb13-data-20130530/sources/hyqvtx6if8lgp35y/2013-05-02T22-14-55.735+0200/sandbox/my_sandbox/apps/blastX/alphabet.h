/**Datei alphabet.h
*In dieser Datei befindet sich die Deklaration der Klasse Alphabet
*/
#include<seqan/basic.h>
#include<seqan/sequence>

using namespace seqan;

#indef ALPHABET_H
#define ALPHABET_H

class Alphabet
{
	private:
		//Gibt Anzahl an verschiedenen Gruppen beim Clustern an
		int groupNumber;
		//k naechste Nachbarn beim Clustern
		int kNearest;
		//Bereits vorhandenes Alphabet nach Bindungskraeften gruppiert
		String<AminoAcid> alphabet;
		append(alphabet,4); // Alanin			0
		append(alphabet,3); // Arginin			1
		append(alphabet,5); // Asparagin		2
		append(alphabet,2); // Asparaginsäure	3
		append(alphabet,1); // Cystein			4
		append(alphabet,5);	// Glutamin			5
		append(alphabet,2); // Glutaminsäure	6
		append(alphabet,4); // Glycin			7
		append(alphabet,3); // Histidin			8
		append(alphabet,4); // Isoleucin		9
		append(alphabet,4); // Leucin			10
		append(alphabet,3); // Lysin			11
		append(alphabet,1); // Methionin		12
		append(alphabet,4); // Phenylalanin		13
		append(alphabet,4);	// Prolin			14
		append(alphabet,5); // Serin			15
		append(alphabet,5); // Threonin			16
		append(alphabet,5); // Tryptophan		17
		append(alphabet,5); // Tyrosin			18
		append(alphabet,4); // Valin			19
		append(alphabet,6);	// STOP				20

		//Speichert neu erstelltes Alphabet
		String<AminoAcid> newAlphabet;


	public:
		//Konstruktor zum generieren eines Alphabet-Objektes
		Alphabet();
		//Ein Destruktor zum zerstoeren eines Alphabet-Objektes
		~Alphabet();
		//Erstellt ein Alphabet nach dem k-nearest clustering 
		/*String<AminoAcid> createNewAlphabet(int groupNumber)
		{
			for (int i = 0; length(blossum62); ++i)
			{
				for (int j = 0;length(blossum62); ++j)
				{
					if (i == j)
						blossum62[i][j] = -10;
				}
				maxWert = max(blossum62)
				if (maxWert == blossum[i][j])
				{
					append(alphabet[i],0);
					append(alphabet[j],0);
				}
			}
			return alphabet;
		}*/
		
		//Gibt bereits vorhandenes Alphabet aus
		String<AminoAcid> getAlphabet()
		{
			return alphabet;
		}

};
#endif
