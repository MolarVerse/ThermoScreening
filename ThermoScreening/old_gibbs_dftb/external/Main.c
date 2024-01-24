#include <stdlib.h>
#include <stdio.h>
#include <string.h>

void schoenflies(int natoms, int *attype, double *coord, char *symbol, double *paramar);
int getsymnum(char *pgroup);

int main(void)
{

    int natoms = 3;
    int attype[3] = {6, 8, 8};
    double coord[9] = {0.09584664, 3.97763572, 0.00000000, 1.35424664, 3.97763572, 0.00000000, -1.16255336, 3.97763572, 0.00000000};
    // double coord[9] = {0.3571, 6.9067, 0.2994, -0.0247, 6.3474, 1.3326, 0.7383, 7.465, -0.7321};
    char symbol[100] = {0};
    double paramar[11];

    paramar[0] = -1;     // Verbosity
    paramar[1] = 10;     // MaxAxisOrder
    paramar[2] = 100;    // MaxOptCycles
    paramar[3] = 1e-3;   // ToleranceSame
    paramar[4] = 5e-1;   // TolerancePrimary
    paramar[5] = 1e-4;   // ToleranceFinal
    paramar[6] = 5e-1;   // MaxOptStep
    paramar[7] = 1e-7;   // MinOptStep
    paramar[8] = 1e-7;   // GradientStep
    paramar[9] = 1e-8;   // OptChangeThreshold
    paramar[10] = 1e-10; // OptChangeThreshold

    schoenflies(natoms, attype, coord, symbol, paramar);

    printf("%s\n%d\n", symbol, getsymnum(symbol));

    return EXIT_SUCCESS;
}
int getsymnum(char *pgroup)
{
    int symnum = 1;
    if (strstr(pgroup, "ci") != NULL)
        symnum = 1;
    if (strstr(pgroup, "ci") != NULL)
        symnum = 1;
    if (strstr(pgroup, "cs") != NULL)
        symnum = 1;
    if (strstr(pgroup, "c2") != NULL)
        symnum = 2;
    if (strstr(pgroup, "c3") != NULL)
        symnum = 3;
    if (strstr(pgroup, "c4") != NULL)
        symnum = 4;
    if (strstr(pgroup, "c5") != NULL)
        symnum = 5;
    if (strstr(pgroup, "c6") != NULL)
        symnum = 6;
    if (strstr(pgroup, "c7") != NULL)
        symnum = 7;
    if (strstr(pgroup, "c8") != NULL)
        symnum = 8;
    if (strstr(pgroup, "c9") != NULL)
        symnum = 9;
    if (strstr(pgroup, "c10") != NULL)
        symnum = 10;
    if (strstr(pgroup, "c11") != NULL)
        symnum = 11;
    if (strstr(pgroup, "s4") != NULL)
        symnum = 2;
    if (strstr(pgroup, "s6") != NULL)
        symnum = 3;
    if (strstr(pgroup, "s8") != NULL)
        symnum = 4;
    if (strstr(pgroup, "d2") != NULL)
        symnum = 4;
    if (strstr(pgroup, "d3") != NULL)
        symnum = 6;
    if (strstr(pgroup, "d4") != NULL)
        symnum = 8;
    if (strstr(pgroup, "d5") != NULL)
        symnum = 10;
    if (strstr(pgroup, "d6") != NULL)
        symnum = 12;
    if (strstr(pgroup, "d7") != NULL)
        symnum = 14;
    if (strstr(pgroup, "d8") != NULL)
        symnum = 16;
    if (strstr(pgroup, "d9") != NULL)
        symnum = 18;
    if (strstr(pgroup, "d10") != NULL)
        symnum = 20;
    if (strstr(pgroup, "t") != NULL)
        symnum = 12;
    if (strstr(pgroup, "th") != NULL)
        symnum = 12;
    if (strstr(pgroup, "td") != NULL)
        symnum = 12;
    if (strstr(pgroup, "o") != NULL)
        symnum = 24;
    if (strstr(pgroup, "oh") != NULL)
        symnum = 24;
    if (strstr(pgroup, "ih") != NULL)
        symnum = 60;
    if (strstr(pgroup, "c") != NULL)
        symnum = 1; //.and.lin
    if (strstr(pgroup, "d") != NULL)
        symnum = 2; //.and.lin
    if (strstr(pgroup, "Ci") != NULL)
        symnum = 1;
    if (strstr(pgroup, "Cs") != NULL)
        symnum = 1;
    if (strstr(pgroup, "C2") != NULL)
        symnum = 2;
    if (strstr(pgroup, "C3") != NULL)
        symnum = 3;
    if (strstr(pgroup, "C4") != NULL)
        symnum = 4;
    if (strstr(pgroup, "C5") != NULL)
        symnum = 5;
    if (strstr(pgroup, "C6") != NULL)
        symnum = 6;
    if (strstr(pgroup, "C7") != NULL)
        symnum = 7;
    if (strstr(pgroup, "C8") != NULL)
        symnum = 8;
    if (strstr(pgroup, "C9") != NULL)
        symnum = 9;
    if (strstr(pgroup, "C10") != NULL)
        symnum = 10;
    if (strstr(pgroup, "C11") != NULL)
        symnum = 11;
    if (strstr(pgroup, "S4") != NULL)
        symnum = 2;
    if (strstr(pgroup, "S6") != NULL)
        symnum = 3;
    if (strstr(pgroup, "S8") != NULL)
        symnum = 4;
    if (strstr(pgroup, "D2") != NULL)
        symnum = 4;
    if (strstr(pgroup, "D3") != NULL)
        symnum = 6;
    if (strstr(pgroup, "D4") != NULL)
        symnum = 8;
    if (strstr(pgroup, "D5") != NULL)
        symnum = 10;
    if (strstr(pgroup, "D6") != NULL)
        symnum = 12;
    if (strstr(pgroup, "D7") != NULL)
        symnum = 14;
    if (strstr(pgroup, "D8") != NULL)
        symnum = 16;
    if (strstr(pgroup, "D9") != NULL)
        symnum = 18;
    if (strstr(pgroup, "D10") != NULL)
        symnum = 20;
    if (strstr(pgroup, "T") != NULL)
        symnum = 12;
    if (strstr(pgroup, "Th") != NULL)
        symnum = 12;
    if (strstr(pgroup, "Td") != NULL)
        symnum = 12;
    if (strstr(pgroup, "O") != NULL)
        symnum = 24;
    if (strstr(pgroup, "Oh") != NULL)
        symnum = 24;
    if (strstr(pgroup, "Ih") != NULL)
        symnum = 60;
    if (strstr(pgroup, "C") != NULL)
        symnum = 1; // .and.lin
    if (strstr(pgroup, "D") != NULL)
        symnum = 2; //.and.lin

    return symnum;
}