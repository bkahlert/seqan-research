// BLASTX IMPLEMENTIERUNG VON ANNKATRIN BRESSIN UND MARJAN FAIZI
// SOFTWAREPROJEKT VOM 2.4. - 29.5.2012
// VERSION VOM 04.MAI.2013

/**@brief output.cpp
*Beinhaltet Funktionen, die die Ergebnisse in eine Textdatei schreiben
*/

#include "own_functions.h"





void WRITE_TO_FILE(Match_found & seed_found, StringSet<String<char> > & proteinID, StringSet<String<char> > & readID,String<unsigned> numb_alp, unsigned & alp){
	/**@brief write_to_file schreibt die gefundenen matches in eine Textdatei
	*@param seed_found Instanz der Klasse Match_found enthaelt die Informationen der gefundenen matches 
	*@param proteinID IDs der einzelnen Proteinsequenzen
	*@param readID IDs der einzelnen Reads
	*@param numb_alp Speichert die Zahlen der ausgesuchten Alphabete
	*@param alp Aktuelle Position des benutzten Alphabets
	*@return Hat kein Rueckgabewert
	*/
	//Aktuelle Zeit speichern
	/*
	time_t rawtime;
	struct tm * timeinfo;
	time (&rawtime);
	timeinfo = localtime (&rawtime);
	string matrix[8] = {"Blosum30","Blosum45","Blosum80","Pam120","Pam200","Pam250","Pam40","Vtml200"};
	string file = "output" + matrix[numb_alp[alp]-49] + ".txt";
	ofstream outfile;
	outfile.open(toCString(file));
    if (outfile.is_open()){
		outfile <<"#Tool: BlastX Version 1"<<endl;
		outfile <<"#Author: Annkatrin Bressin und Marjan Faizi"<<endl;
		outfile <<"#Scorematrix: "<<matrix[numb_alp[alp]-49]<<endl;		
		outfile <<"#Date: "<<asctime(timeinfo)<<endl; 
		outfile << "read_id\tbegin_read\tend_read\tprotein_id\tbegin_protein\tend_protein\tframe\tside"<<endl;  
		for (unsigned position = 0; position < length(seed_found.position_read); ++position){
			unsigned read = (seed_found.position_read[position]/6);
			unsigned protein = seed_found.position_protein[position];
			unsigned frame = seed_found.position_read[position]%6;
			unsigned side = frame/3;
			outfile <<readID[read]
			<<"\t"<<seed_found.begin_read[position]
			<<"\t"<<seed_found.end_read[position]
			<<"\t"<<proteinID[protein]
			<<"\t"<<seed_found.begin_protein[position]
			<<"\t"<<seed_found.end_protein[position]
			<<"\t"<<frame
			<<"\t"<<side<<endl;
		}
		outfile.close();
	}
	*/
}


//void WRITE_TO_FILE(Match_found & verify_found, StringSet<String<char> > & proteinID, StringSet<String<char> > & readID, StringSet<String<Dna> > & reads, String<unsigned> numb_alp, unsigned & alp){
	/**@brief write_to_file schreibt die verifizierten matches in eine Textdatei
	*@param verify_found Instanz der Klasse Match_found enthaelt die Informationen der verifizierten matches 
	*@param proteinID IDs der einzelnen Proteinsequenzen
	*@param readID IDs der einzelnen Reads
	*@param reads Original Reads(nicht uebersetzt)
	*@param numb_alp Speichert die Zahlen der ausgesuchten Alphabete
	*@param alp Aktuelle Position des benutzten Alphabets
	*@return Hat kein Rueckgabewert
	*/
	//Aktuelle Zeit speichern
/*	
time_t rawtime;
	struct tm * timeinfo;
	time (&rawtime);
	timeinfo = localtime (&rawtime);
	string matrix[8] = {"Blosum30","Blosum45","Blosum80","Pam120","Pam200","Pam250","Pam40","Vtml200"};
	string file = "output" + matrix[numb_alp[alp]-49] + "verify.txt";
	ofstream outfile;
	outfile.open(toCString(file));
    if (outfile.is_open()){
		outfile <<"#Tool: BlastX Version 1"<<endl;
		outfile <<"#Author: Annkatrin Bressin und Marjan Faizi"<<endl;
		outfile <<"#Scorematrix: "<<matrix[numb_alp[alp]-49]<<endl;		
		outfile <<"#Date: "<<asctime(timeinfo)<<endl;
		outfile << "read_id\tread_sequence\tprotein_id\tbegin_protein\tend_protein\tframe\tside\tscore"<<endl;  
		for (unsigned position = 0; position < length(verify_found.position_read); ++position){
			unsigned read = verify_found.position_read[position]/6;
			unsigned protein = verify_found.position_protein[position];
			unsigned frame = verify_found.position_read[position]%6;
			unsigned side = frame/3;
			outfile << readID[read]
			<<"\t"<<reads[read]
			<<"\t"<<proteinID[protein]
			<<"\t"<<verify_found.begin_protein[position]
			<<"\t"<<verify_found.end_protein[position]
			<<"\t"<<frame
			<<"\t"<<side
			<<"\t"<<verify_found.score[position]<<endl;
		}
		outfile.close();
	}
}
*/
