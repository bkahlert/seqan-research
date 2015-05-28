
/*********************************************************
*   Grundlagen der Programmierung / Informatik 2
*   BESCHREIBUNG:	Toleranzmaße nach DIN ISO 2768-1
*					Vorlage für die praktische Aufgabe
*   Autor:          Udo Mustermann
*   Datum           12.4.2012
**********************************************************/

/*
Werte gemäß DIN ISO 2768-1

Nennmaßbereich          Toleranzklassen, Toleranzen in mm
in mm
                        f	m	c	v
0.50    3.00		0.05	0.10	0.15	-
3.00    6.00		0.05	0.10	0.20	0.50
6.00    30	        0.10	0.20	0.50	1.00
30      120		0.15	0.30	0.80	1.50
120     400		0.20	0.50	1.20	2.50
400     1000		0.30	0.80	2.00	4.00
1000    2000		0.50	1.20	3.00	6.00
2000    4000		-	2.00	4.00	8.00
*/


#include <stdio.h>

#ifdef _WIN32
#include <windows.h>    // notwendig für SetConsoleOutputCP(), siehe unten
#endif

int main(void)
{
	int row=0;
	int column=0;
	double DIN_Werte [6][8] = 	{{0.5,3.0,6.0,30.0,120.0,400.0,1000.0,2000.0},
					{3.0,6.0,30.0,120.0,400.0,1000.0,2000.0,4000.0},
					{0.05,0.05,0.10,0.15,0.20,0.30,0.50,0.0},
					{0.10,0.10,0.20,0.30,0.50,0.80,1.20,2.0},
					{0.15,0.20,0.50,0.80,1.20,2.0,3.0,4.0},
					{0.0,0.5,1.0,1.5,2.5,4.0,6.0,8.0}};


#ifdef _WIN32
    SetConsoleOutputCP(1252);     // Setzt unter Windows den ANSI-Zeichensatz auch in der
#endif                           // Konsole ()
    // falls unter Eigenschaften /Schriftart:
    // True-Type-Schriftart ausgewählt wurde
start:
    double laenge = 0;
    char toleranz = 'x';
    //printf("\nNN Berechnung des Höchstmaß und Mindestmaß für Längenmaße \n");
    //printf("ohne einzelne Toleranzeintragung nach DIN ISO 2768-1 \n");
    //printf("Dieses Programm wurde erstellt von Udo Mustermann\n\n");

        // Eingabe
        printf("\nGeben Sie die Länge in mm und die Toleranzklasse (f, m, c, v) ein: ");
        int ret = scanf("%lf %c", &laenge, &toleranz);
	if(ret!=2)
	{
		printf("Falsche Eingabe.\n\n");
		goto start;
	}        
	printf("\nEingegebene Länge: %.2f, Toleranzklasse: %c\n", laenge, toleranz);

	for(int i=0;i<8;i++)
	{
		if(DIN_Werte[0][i]<=laenge && DIN_Werte[1][i]>=laenge)
		{
			row=i;
			break;
		}
		else if(DIN_Werte[0][0]>laenge || DIN_Werte[1][7]<laenge)
		{
			printf("Eingegebene Länge entspricht nicht den DIN ISO 2768-1 Anforderungen.\nBitte um erneute Eingabe\n\n");
			goto start;
		}
	}
		
	if((row==0 && toleranz=='v')||(row==7 && toleranz=='f'))
	{
			printf("Eingegebene Toleranzklasse laut DIN ISO 2768-1 nicht für eingegebene Länge vorhanden..\nBitte um erneute Eingabe\n\n");
			goto start;		
	} 

        if(toleranz == 'f')
        {
            column=2;
        }
        else if(toleranz == 'm') 
        {
            column=3;
        }
        else if(toleranz == 'm') 
        {
            column=4;
        }
        else if(toleranz == 'm') 
        {
            column=5;
        }
	else
	{
		printf("Kein Toleranzmaß gefunden.\nBitte um erneute Eingabe\n\n");
		goto start;
	}

	printf("\n\nRESULTAT:\nMindestmaß:  %.2f, Höchstmaß: %.2f\n", laenge - DIN_Werte[column][row], laenge + DIN_Werte[column][row] );
	getchar();
	getchar();

}

