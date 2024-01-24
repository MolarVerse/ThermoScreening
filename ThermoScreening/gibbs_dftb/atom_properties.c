/* Mass of the elements. Taken from www.webelements.com */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>

struct table {double mass; char element[5]; };
struct list  {int periodic_number; char element[5]; };


double init_mass(char name[]){

    int i        = 0;
    int elements = 0;
    double mass  = 0;

    while (name[i])
    {
        name[i] = tolower(name[i]);
        i++;
    }


    struct table molar_mass[] =
    {
        {1.00784      , "h"    },
        {2.01410175   , "d"    },
        {3.0160495    , "t"    },
        {4.002602     , "he"   },
        {6.941        , "li"   },
        {9.012182     , "be"   },
        {10.811       , "b"    },
        {12.0107      , "c"    },
        {14.0067      , "n"    },
        {15.999       , "o"    },
        {18.998403    , "f"    },
        {20.1797      , "ne"   },

        {22.989769    , "na"   },
        {24.3050      , "mg"   },
        {26.981538    , "al"   },
        {28.0855      , "si"   },
        {30.973762    , "p"    },
        {32.065       , "s"    },
        {35.453       , "cl"   },
        {39.948       , "ar"   },

        {39.0983      , "k"    },
        {40.078       , "ca"   },
        {44.955910    , "sc"   },
        {47.880       , "ti"   },
        {50.9415      , "v"    },
        {51.9961      , "cr"   },
        {54.938049    , "mn"   },
        {55.845       , "fe"   },
        {58.933200    , "co"   },
        {58.6934      , "ni"   },
        {63.546       , "cu"   },
        {65.399       , "zn"   },
        {69.723       , "ga"   },
        {72.64        , "ge"   },
        {74.92160     , "as"   },
        {78.96        , "se"   },
        {79.904       , "br"   },
        {83.798       , "kr"   },

        {85.4678      , "rb"   },
        {87.62        , "sr"   },
        {88.90585     , "y"    },
        {91.224       , "zr"   },
        {92.90638     , "nb"   },
        {95.94        , "mo"   },
        {98.9063      , "tc"   },
        {101.07       , "ru"   },
        {102.9055     , "rh"   },
        {106.42       , "pd"   },
        {107.8682     , "ag"   },
        {112.411      , "cd"   },
        {114.818      , "in"   },
        {118.71       , "sn"   },
        {121.76       , "sb"   },
        {127.6        , "te"   },
        {126.90447    , "i"    },
        {131.293      , "xe"   },

        {132.90546    , "cs"   },
        {137.327      , "ba"   },
        {138.9055     , "la"   },
        {140.116      , "ce"   },
        {140.90765    , "pr"   },
        {144.24       , "nd"   },
        {146.9151     , "pm"   },
        {150.36       , "sm"   },
        {151.964      , "eu"   },
        {157.25       , "gd"   },
        {158.92534    , "tb"   },
        {162.5        , "dy"   },
        {164.93032    , "ho"   },
        {167.259      , "er"   },
        {168.93421    , "tm"   },
        {173.04       , "yb"   },
        {174.967      , "lu"   },
        {178.49       , "hf"   },
        {180.9479     , "ta"   },
        {183.84       , "w"    },
        {186.207      , "re"   },
        {190.23       , "os"   },
        {192.217      , "ir"   },
        {195.078      , "pt"   },
        {196.96655    , "au"   },
        {200.59       , "hg"   },
        {204.3833     , "tl"   },
        {207.2        , "pb"   },
        {208.98038    , "bi"   },
        {208.9824     , "po"   },
        {209.9871     , "at"   },
        {222.0176     , "rn"   },

        {223.0197     , "fr"   },
        {226.0254     , "ra"   },
        {227.0278     , "ac"   },
        {232.0381     , "th"   },
        {231.03588    , "pa"   },
        {238.0289     , "u"    },
        {237.0482     , "np"   },
        {244.0642     , "pu"   },
        {243.0614     , "am"   },
        {247.0703     , "cm"   },
        {247.0703     , "bk"   },
        {251.0796     , "cf"   },
        {252.0829     , "es"   },
        {257.0951     , "fm"   },
        {258.0986     , "md"   },
        {259.1009     , "no"   },
        {260.1053     , "lr"   },

        // turbomole dummy atom
        {999.00000    , "q"    },
        // gaussian dummy atom
        {999.00000    , "x"    },
        // Cavity particle
        {1000.00000   , "cav"  },
        // super-heavy particle
        {1000000.0    , "sup"  },
        // {not available, "db-uuo" },
        // light dummy
        {1.0          , "dum"  },
    };

    elements = sizeof(molar_mass)/sizeof(molar_mass[0]);

    for (i = 0; i < elements; i++)
    {
        if(strcasecmp(name,molar_mass[i].element) == 0)
        {
            mass = molar_mass[i].mass;
            break;
        }
    }

    return (mass);

}


int periodic_number(char name[])
{
    int i        = 0;
    int elements = 0;

    int number   = 0;

    while (name[i])
    {
        name[i] = tolower(name[i]);
        i++;
    }


    struct list periodic_number[] =
    {
        {   1   , "h"  },
        {   1   , "d"  },
        {   1   , "t"  },
        {   2   , "he" },

        {   3   , "li" },
        {   4   , "be" },
        {   5   , "b"  },
        {   6   , "c"  },
        {   7   , "n"  },
        {   8   , "o"  },
        {   9   , "f"  },
        {  10   , "ne" },

        {  11   , "na" },
        {  12   , "mg" },
        {  13   , "al" },
        {  14   , "si" },
        {  15   , "p"  },
        {  16   , "s"  },
        {  17   , "cl" },
        {  18   , "ar" },

        {  19   , "k"  },
        {  20   , "ca" },
        {  21   , "sc" },
        {  22   , "ti" },
        {  23   , "v"  },
        {  24   , "cr" },
        {  25   , "mn" },
        {  26   , "fe" },
        {  27   , "co" },
        {  28   , "ni" },
        {  29   , "cu" },
        {  30   , "zn" },
        {  31   , "ga" },
        {  32   , "ge" },
        {  33   , "as" },
        {  34   , "se" },
        {  35   , "br" },
        {  36   , "kr" },

        {  37   , "rb" },
        {  38   , "sr" },
        {  39   , "y"  },
        {  40   , "zr" },
        {  41   , "nb" },
        {  42   , "mo" },
        {  43   , "tc" },
        {  44   , "ru" },
        {  45   , "rh" },
        {  46   , "pd" },
        {  47   , "ag" },
        {  48   , "cd" },
        {  49   , "in" },
        {  50   , "sn" },
        {  51   , "sb" },
        {  52   , "te" },
        {  53   , "i"  },
        {  54   , "xe" },

        {  55   , "cs" },
        {  56   , "ba" },
        {  57   , "la" },
        {  58   , "ce" },
        {  59   , "pr" },
        {  60   , "nd" },
        {  61   , "pm" },
        {  62   , "sm" },
        {  63   , "eu" },
        {  64   , "gd" },
        {  65   , "tb" },
        {  66   , "dy" },
        {  67   , "ho" },
        {  68   , "er" },
        {  69   , "tm" },
        {  70   , "yb" },
        {  71   , "lu" },
        {  72   , "hf" },
        {  73   , "ta" },
        {  74   , "w"  },
        {  75   , "re" },
        {  76   , "os" },
        {  77   , "ir" },
        {  78   , "pt" },
        {  79   , "au" },
        {  80   , "hg" },
        {  81   , "tl" },
        {  82   , "pb" },
        {  83   , "bi" },
        {  84   , "po" },
        {  85   , "at" },
        {  86   , "rn" },

        {  87   , "fr" },
        {  88   , "ra" },
        {  89   , "ac" },
        {  90   , "th" },
        {  91   , "pa" },
        {  92   , "u"  },
        {  93   , "np" },
        {  94   , "pu" },
        {  95   , "am" },
        {  96   , "cm" },
        {  97   , "bk" },
        {  98   , "cf" },
        {  99   , "es" },
        { 100   , "fm" },
        { 101   , "md" },
        { 102   , "no" },
        { 103   , "lr" },

        // turbomole dummy atom
        { 999   , "q"  },
        // gaussian dummy atom
        { 999   , "x"  },
        // Cavity particle
        {1000   , "cav"},
        // super-heavy particle
        {1000000, "sup"},
        // light dummy
        {1      , "dum"},
        // {not available, "db-uuo" },
    };

    elements = sizeof(periodic_number)/sizeof(periodic_number[0]);

    for (i = 0; i < elements; i++)
    {
        if(strcasecmp(name,periodic_number[i].element) == 0)
        {
            number = periodic_number[i].periodic_number;
            break;
        }
    }

    return (number);
}
