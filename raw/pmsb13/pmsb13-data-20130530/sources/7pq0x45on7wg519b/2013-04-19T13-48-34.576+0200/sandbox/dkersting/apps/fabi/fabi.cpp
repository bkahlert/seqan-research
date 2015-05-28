
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
    double laenge = 0;
    char toleranz = 'x';
	int row;
	int column;
	float DIN_Werte [6][8] = {{0.5,3.0,6.0,30.0,120.0,400.0,1000.0,2000.0},
{3.0,6.0,30.0,120.0,400.0,1000.0,2000.0,4000.0},
{0.05,0.05,0.10,0.15,0.20,0.30,0.50,NULL},
{0.10,0.10,0.20,0.30,0.50,0.80,1.20,2.0},
{0.15,0.20,0.50,0.80,1.20,2.0,3.0,4.0},
{NULL,0.5,1.0,1.5,2.5,4.0,6.0,8.0}}
;
//	DIN_Werte[0][]={0.5,3.0,6.0,30.0,120.0,400.0,1000.0,2000.0};
//	DIN_Werte[1][]={3.0,6.0,30.0,120.0,400.0,1000.0,2000.0,4000.0};
//	DIN_Werte[2][]={0.05,0.05,0.10,0.15,0.20,0.30,0.50,NULL};
//	DIN_Werte[3][]={0.10,0.10,0.20,0.30,0.50,0.80,1.20,2.0};
//	DIN_Werte[4][]={0.15,0.20,0.50,0.80,1.20,2.0,3.0,4.0};
//	DIN_Werte[5][]={NULL,0.5,1.0,1.5,2.5,4.0,6.0,8.0};

#ifdef _WIN32
    SetConsoleOutputCP(1252);     // Setzt unter Windows den ANSI-Zeichensatz auch in der
#endif                           // Konsole ()
    // falls unter Eigenschaften /Schriftart:
    // True-Type-Schriftart ausgewählt wurde
start:
    //printf("\nNN Berechnung des Höchstmaß und Mindestmaß für Längenmaße \n");
    //printf("ohne einzelne Toleranzeintragung nach DIN ISO 2768-1 \n");
    //printf("Dieses Programm wurde erstellt von Udo Mustermann\n\n");

        // Eingabe
        printf("\nGeben Sie die Länge in mm und die Toleranzklasse (f, m, c, v) ein: ");
        scanf("%lf %c", &laenge, &toleranz);
        printf("\nEingegebene Länge: %f, Toleranzklasse: %c\n", laenge, toleranz);

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
		
	if(row==0 && toleranz=='v'||row==7 && toleranz=='f')
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

printf("\n\n#################\n%f\n%f\n###################\n\n",row,column);
	printf("\n\nRESULTAT:\nMindestmaß:  %f, Höchstmaß: %f\n", laenge - DIN_Werte[column][row], laenge + DIN_Werte[column][row] );
	getchar();
	getchar();

    return 0;
}

