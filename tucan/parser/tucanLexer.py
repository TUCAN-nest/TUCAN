# Generated from tucan.g4 by ANTLR 4.11.1
from antlr4 import *
from io import StringIO
import sys
if sys.version_info[1] > 5:
    from typing import TextIO
else:
    from typing.io import TextIO


def serializedATN():
    return [
        4,0,120,589,6,-1,2,0,7,0,2,1,7,1,2,2,7,2,2,3,7,3,2,4,7,4,2,5,7,5,
        2,6,7,6,2,7,7,7,2,8,7,8,2,9,7,9,2,10,7,10,2,11,7,11,2,12,7,12,2,
        13,7,13,2,14,7,14,2,15,7,15,2,16,7,16,2,17,7,17,2,18,7,18,2,19,7,
        19,2,20,7,20,2,21,7,21,2,22,7,22,2,23,7,23,2,24,7,24,2,25,7,25,2,
        26,7,26,2,27,7,27,2,28,7,28,2,29,7,29,2,30,7,30,2,31,7,31,2,32,7,
        32,2,33,7,33,2,34,7,34,2,35,7,35,2,36,7,36,2,37,7,37,2,38,7,38,2,
        39,7,39,2,40,7,40,2,41,7,41,2,42,7,42,2,43,7,43,2,44,7,44,2,45,7,
        45,2,46,7,46,2,47,7,47,2,48,7,48,2,49,7,49,2,50,7,50,2,51,7,51,2,
        52,7,52,2,53,7,53,2,54,7,54,2,55,7,55,2,56,7,56,2,57,7,57,2,58,7,
        58,2,59,7,59,2,60,7,60,2,61,7,61,2,62,7,62,2,63,7,63,2,64,7,64,2,
        65,7,65,2,66,7,66,2,67,7,67,2,68,7,68,2,69,7,69,2,70,7,70,2,71,7,
        71,2,72,7,72,2,73,7,73,2,74,7,74,2,75,7,75,2,76,7,76,2,77,7,77,2,
        78,7,78,2,79,7,79,2,80,7,80,2,81,7,81,2,82,7,82,2,83,7,83,2,84,7,
        84,2,85,7,85,2,86,7,86,2,87,7,87,2,88,7,88,2,89,7,89,2,90,7,90,2,
        91,7,91,2,92,7,92,2,93,7,93,2,94,7,94,2,95,7,95,2,96,7,96,2,97,7,
        97,2,98,7,98,2,99,7,99,2,100,7,100,2,101,7,101,2,102,7,102,2,103,
        7,103,2,104,7,104,2,105,7,105,2,106,7,106,2,107,7,107,2,108,7,108,
        2,109,7,109,2,110,7,110,2,111,7,111,2,112,7,112,2,113,7,113,2,114,
        7,114,2,115,7,115,2,116,7,116,2,117,7,117,2,118,7,118,2,119,7,119,
        1,0,1,0,1,1,1,1,1,1,1,2,1,2,1,2,1,3,1,3,1,3,1,4,1,4,1,5,1,5,1,6,
        1,6,1,7,1,7,1,8,1,8,1,9,1,9,1,9,1,10,1,10,1,10,1,11,1,11,1,11,1,
        12,1,12,1,12,1,13,1,13,1,13,1,14,1,14,1,15,1,15,1,16,1,16,1,16,1,
        17,1,17,1,17,1,18,1,18,1,19,1,19,1,19,1,20,1,20,1,20,1,21,1,21,1,
        21,1,22,1,22,1,23,1,23,1,23,1,24,1,24,1,24,1,25,1,25,1,25,1,26,1,
        26,1,26,1,27,1,27,1,27,1,28,1,28,1,28,1,29,1,29,1,29,1,30,1,30,1,
        30,1,31,1,31,1,31,1,32,1,32,1,32,1,33,1,33,1,33,1,34,1,34,1,34,1,
        35,1,35,1,35,1,36,1,36,1,36,1,37,1,37,1,37,1,38,1,38,1,39,1,39,1,
        39,1,40,1,40,1,40,1,41,1,41,1,41,1,42,1,42,1,42,1,43,1,43,1,43,1,
        44,1,44,1,44,1,45,1,45,1,45,1,46,1,46,1,46,1,47,1,47,1,47,1,48,1,
        48,1,48,1,49,1,49,1,49,1,50,1,50,1,50,1,51,1,51,1,51,1,52,1,52,1,
        53,1,53,1,53,1,54,1,54,1,54,1,55,1,55,1,55,1,56,1,56,1,56,1,57,1,
        57,1,57,1,58,1,58,1,58,1,59,1,59,1,59,1,60,1,60,1,60,1,61,1,61,1,
        61,1,62,1,62,1,62,1,63,1,63,1,63,1,64,1,64,1,64,1,65,1,65,1,65,1,
        66,1,66,1,66,1,67,1,67,1,67,1,68,1,68,1,68,1,69,1,69,1,69,1,70,1,
        70,1,70,1,71,1,71,1,71,1,72,1,72,1,72,1,73,1,73,1,74,1,74,1,74,1,
        75,1,75,1,75,1,76,1,76,1,76,1,77,1,77,1,77,1,78,1,78,1,78,1,79,1,
        79,1,79,1,80,1,80,1,80,1,81,1,81,1,81,1,82,1,82,1,82,1,83,1,83,1,
        83,1,84,1,84,1,84,1,85,1,85,1,85,1,86,1,86,1,86,1,87,1,87,1,87,1,
        88,1,88,1,88,1,89,1,89,1,89,1,90,1,90,1,90,1,91,1,91,1,92,1,92,1,
        92,1,93,1,93,1,93,1,94,1,94,1,94,1,95,1,95,1,95,1,96,1,96,1,96,1,
        97,1,97,1,97,1,98,1,98,1,98,1,99,1,99,1,99,1,100,1,100,1,100,1,101,
        1,101,1,101,1,102,1,102,1,102,1,103,1,103,1,103,1,104,1,104,1,104,
        1,105,1,105,1,105,1,106,1,106,1,106,1,107,1,107,1,107,1,108,1,108,
        1,108,1,109,1,109,1,109,1,110,1,110,1,110,1,111,1,111,1,111,1,112,
        1,112,1,112,1,113,1,113,1,113,1,114,1,114,1,114,1,115,1,115,1,115,
        1,116,1,116,1,116,1,117,1,117,1,117,1,118,1,118,1,119,1,119,4,119,
        586,8,119,11,119,12,119,587,0,0,120,1,1,3,2,5,3,7,4,9,5,11,6,13,
        7,15,8,17,9,19,10,21,11,23,12,25,13,27,14,29,15,31,16,33,17,35,18,
        37,19,39,20,41,21,43,22,45,23,47,24,49,25,51,26,53,27,55,28,57,29,
        59,30,61,31,63,32,65,33,67,34,69,35,71,36,73,37,75,38,77,39,79,40,
        81,41,83,42,85,43,87,44,89,45,91,46,93,47,95,48,97,49,99,50,101,
        51,103,52,105,53,107,54,109,55,111,56,113,57,115,58,117,59,119,60,
        121,61,123,62,125,63,127,64,129,65,131,66,133,67,135,68,137,69,139,
        70,141,71,143,72,145,73,147,74,149,75,151,76,153,77,155,78,157,79,
        159,80,161,81,163,82,165,83,167,84,169,85,171,86,173,87,175,88,177,
        89,179,90,181,91,183,92,185,93,187,94,189,95,191,96,193,97,195,98,
        197,99,199,100,201,101,203,102,205,103,207,104,209,105,211,106,213,
        107,215,108,217,109,219,110,221,111,223,112,225,113,227,114,229,
        115,231,116,233,117,235,118,237,119,239,120,1,0,3,1,0,50,57,1,0,
        49,57,1,0,48,57,589,0,1,1,0,0,0,0,3,1,0,0,0,0,5,1,0,0,0,0,7,1,0,
        0,0,0,9,1,0,0,0,0,11,1,0,0,0,0,13,1,0,0,0,0,15,1,0,0,0,0,17,1,0,
        0,0,0,19,1,0,0,0,0,21,1,0,0,0,0,23,1,0,0,0,0,25,1,0,0,0,0,27,1,0,
        0,0,0,29,1,0,0,0,0,31,1,0,0,0,0,33,1,0,0,0,0,35,1,0,0,0,0,37,1,0,
        0,0,0,39,1,0,0,0,0,41,1,0,0,0,0,43,1,0,0,0,0,45,1,0,0,0,0,47,1,0,
        0,0,0,49,1,0,0,0,0,51,1,0,0,0,0,53,1,0,0,0,0,55,1,0,0,0,0,57,1,0,
        0,0,0,59,1,0,0,0,0,61,1,0,0,0,0,63,1,0,0,0,0,65,1,0,0,0,0,67,1,0,
        0,0,0,69,1,0,0,0,0,71,1,0,0,0,0,73,1,0,0,0,0,75,1,0,0,0,0,77,1,0,
        0,0,0,79,1,0,0,0,0,81,1,0,0,0,0,83,1,0,0,0,0,85,1,0,0,0,0,87,1,0,
        0,0,0,89,1,0,0,0,0,91,1,0,0,0,0,93,1,0,0,0,0,95,1,0,0,0,0,97,1,0,
        0,0,0,99,1,0,0,0,0,101,1,0,0,0,0,103,1,0,0,0,0,105,1,0,0,0,0,107,
        1,0,0,0,0,109,1,0,0,0,0,111,1,0,0,0,0,113,1,0,0,0,0,115,1,0,0,0,
        0,117,1,0,0,0,0,119,1,0,0,0,0,121,1,0,0,0,0,123,1,0,0,0,0,125,1,
        0,0,0,0,127,1,0,0,0,0,129,1,0,0,0,0,131,1,0,0,0,0,133,1,0,0,0,0,
        135,1,0,0,0,0,137,1,0,0,0,0,139,1,0,0,0,0,141,1,0,0,0,0,143,1,0,
        0,0,0,145,1,0,0,0,0,147,1,0,0,0,0,149,1,0,0,0,0,151,1,0,0,0,0,153,
        1,0,0,0,0,155,1,0,0,0,0,157,1,0,0,0,0,159,1,0,0,0,0,161,1,0,0,0,
        0,163,1,0,0,0,0,165,1,0,0,0,0,167,1,0,0,0,0,169,1,0,0,0,0,171,1,
        0,0,0,0,173,1,0,0,0,0,175,1,0,0,0,0,177,1,0,0,0,0,179,1,0,0,0,0,
        181,1,0,0,0,0,183,1,0,0,0,0,185,1,0,0,0,0,187,1,0,0,0,0,189,1,0,
        0,0,0,191,1,0,0,0,0,193,1,0,0,0,0,195,1,0,0,0,0,197,1,0,0,0,0,199,
        1,0,0,0,0,201,1,0,0,0,0,203,1,0,0,0,0,205,1,0,0,0,0,207,1,0,0,0,
        0,209,1,0,0,0,0,211,1,0,0,0,0,213,1,0,0,0,0,215,1,0,0,0,0,217,1,
        0,0,0,0,219,1,0,0,0,0,221,1,0,0,0,0,223,1,0,0,0,0,225,1,0,0,0,0,
        227,1,0,0,0,0,229,1,0,0,0,0,231,1,0,0,0,0,233,1,0,0,0,0,235,1,0,
        0,0,0,237,1,0,0,0,0,239,1,0,0,0,1,241,1,0,0,0,3,243,1,0,0,0,5,246,
        1,0,0,0,7,249,1,0,0,0,9,252,1,0,0,0,11,254,1,0,0,0,13,256,1,0,0,
        0,15,258,1,0,0,0,17,260,1,0,0,0,19,262,1,0,0,0,21,265,1,0,0,0,23,
        268,1,0,0,0,25,271,1,0,0,0,27,274,1,0,0,0,29,277,1,0,0,0,31,279,
        1,0,0,0,33,281,1,0,0,0,35,284,1,0,0,0,37,287,1,0,0,0,39,289,1,0,
        0,0,41,292,1,0,0,0,43,295,1,0,0,0,45,298,1,0,0,0,47,300,1,0,0,0,
        49,303,1,0,0,0,51,306,1,0,0,0,53,309,1,0,0,0,55,312,1,0,0,0,57,315,
        1,0,0,0,59,318,1,0,0,0,61,321,1,0,0,0,63,324,1,0,0,0,65,327,1,0,
        0,0,67,330,1,0,0,0,69,333,1,0,0,0,71,336,1,0,0,0,73,339,1,0,0,0,
        75,342,1,0,0,0,77,345,1,0,0,0,79,347,1,0,0,0,81,350,1,0,0,0,83,353,
        1,0,0,0,85,356,1,0,0,0,87,359,1,0,0,0,89,362,1,0,0,0,91,365,1,0,
        0,0,93,368,1,0,0,0,95,371,1,0,0,0,97,374,1,0,0,0,99,377,1,0,0,0,
        101,380,1,0,0,0,103,383,1,0,0,0,105,386,1,0,0,0,107,388,1,0,0,0,
        109,391,1,0,0,0,111,394,1,0,0,0,113,397,1,0,0,0,115,400,1,0,0,0,
        117,403,1,0,0,0,119,406,1,0,0,0,121,409,1,0,0,0,123,412,1,0,0,0,
        125,415,1,0,0,0,127,418,1,0,0,0,129,421,1,0,0,0,131,424,1,0,0,0,
        133,427,1,0,0,0,135,430,1,0,0,0,137,433,1,0,0,0,139,436,1,0,0,0,
        141,439,1,0,0,0,143,442,1,0,0,0,145,445,1,0,0,0,147,448,1,0,0,0,
        149,450,1,0,0,0,151,453,1,0,0,0,153,456,1,0,0,0,155,459,1,0,0,0,
        157,462,1,0,0,0,159,465,1,0,0,0,161,468,1,0,0,0,163,471,1,0,0,0,
        165,474,1,0,0,0,167,477,1,0,0,0,169,480,1,0,0,0,171,483,1,0,0,0,
        173,486,1,0,0,0,175,489,1,0,0,0,177,492,1,0,0,0,179,495,1,0,0,0,
        181,498,1,0,0,0,183,501,1,0,0,0,185,503,1,0,0,0,187,506,1,0,0,0,
        189,509,1,0,0,0,191,512,1,0,0,0,193,515,1,0,0,0,195,518,1,0,0,0,
        197,521,1,0,0,0,199,524,1,0,0,0,201,527,1,0,0,0,203,530,1,0,0,0,
        205,533,1,0,0,0,207,536,1,0,0,0,209,539,1,0,0,0,211,542,1,0,0,0,
        213,545,1,0,0,0,215,548,1,0,0,0,217,551,1,0,0,0,219,554,1,0,0,0,
        221,557,1,0,0,0,223,560,1,0,0,0,225,563,1,0,0,0,227,566,1,0,0,0,
        229,569,1,0,0,0,231,572,1,0,0,0,233,575,1,0,0,0,235,578,1,0,0,0,
        237,581,1,0,0,0,239,583,1,0,0,0,241,242,5,72,0,0,242,2,1,0,0,0,243,
        244,5,72,0,0,244,245,5,101,0,0,245,4,1,0,0,0,246,247,5,76,0,0,247,
        248,5,105,0,0,248,6,1,0,0,0,249,250,5,66,0,0,250,251,5,101,0,0,251,
        8,1,0,0,0,252,253,5,66,0,0,253,10,1,0,0,0,254,255,5,67,0,0,255,12,
        1,0,0,0,256,257,5,78,0,0,257,14,1,0,0,0,258,259,5,79,0,0,259,16,
        1,0,0,0,260,261,5,70,0,0,261,18,1,0,0,0,262,263,5,78,0,0,263,264,
        5,101,0,0,264,20,1,0,0,0,265,266,5,78,0,0,266,267,5,97,0,0,267,22,
        1,0,0,0,268,269,5,77,0,0,269,270,5,103,0,0,270,24,1,0,0,0,271,272,
        5,65,0,0,272,273,5,108,0,0,273,26,1,0,0,0,274,275,5,83,0,0,275,276,
        5,105,0,0,276,28,1,0,0,0,277,278,5,80,0,0,278,30,1,0,0,0,279,280,
        5,83,0,0,280,32,1,0,0,0,281,282,5,67,0,0,282,283,5,108,0,0,283,34,
        1,0,0,0,284,285,5,65,0,0,285,286,5,114,0,0,286,36,1,0,0,0,287,288,
        5,75,0,0,288,38,1,0,0,0,289,290,5,67,0,0,290,291,5,97,0,0,291,40,
        1,0,0,0,292,293,5,83,0,0,293,294,5,99,0,0,294,42,1,0,0,0,295,296,
        5,84,0,0,296,297,5,105,0,0,297,44,1,0,0,0,298,299,5,86,0,0,299,46,
        1,0,0,0,300,301,5,67,0,0,301,302,5,114,0,0,302,48,1,0,0,0,303,304,
        5,77,0,0,304,305,5,110,0,0,305,50,1,0,0,0,306,307,5,70,0,0,307,308,
        5,101,0,0,308,52,1,0,0,0,309,310,5,67,0,0,310,311,5,111,0,0,311,
        54,1,0,0,0,312,313,5,78,0,0,313,314,5,105,0,0,314,56,1,0,0,0,315,
        316,5,67,0,0,316,317,5,117,0,0,317,58,1,0,0,0,318,319,5,90,0,0,319,
        320,5,110,0,0,320,60,1,0,0,0,321,322,5,71,0,0,322,323,5,97,0,0,323,
        62,1,0,0,0,324,325,5,71,0,0,325,326,5,101,0,0,326,64,1,0,0,0,327,
        328,5,65,0,0,328,329,5,115,0,0,329,66,1,0,0,0,330,331,5,83,0,0,331,
        332,5,101,0,0,332,68,1,0,0,0,333,334,5,66,0,0,334,335,5,114,0,0,
        335,70,1,0,0,0,336,337,5,75,0,0,337,338,5,114,0,0,338,72,1,0,0,0,
        339,340,5,82,0,0,340,341,5,98,0,0,341,74,1,0,0,0,342,343,5,83,0,
        0,343,344,5,114,0,0,344,76,1,0,0,0,345,346,5,89,0,0,346,78,1,0,0,
        0,347,348,5,90,0,0,348,349,5,114,0,0,349,80,1,0,0,0,350,351,5,78,
        0,0,351,352,5,98,0,0,352,82,1,0,0,0,353,354,5,77,0,0,354,355,5,111,
        0,0,355,84,1,0,0,0,356,357,5,84,0,0,357,358,5,99,0,0,358,86,1,0,
        0,0,359,360,5,82,0,0,360,361,5,117,0,0,361,88,1,0,0,0,362,363,5,
        82,0,0,363,364,5,104,0,0,364,90,1,0,0,0,365,366,5,80,0,0,366,367,
        5,100,0,0,367,92,1,0,0,0,368,369,5,65,0,0,369,370,5,103,0,0,370,
        94,1,0,0,0,371,372,5,67,0,0,372,373,5,100,0,0,373,96,1,0,0,0,374,
        375,5,73,0,0,375,376,5,110,0,0,376,98,1,0,0,0,377,378,5,83,0,0,378,
        379,5,110,0,0,379,100,1,0,0,0,380,381,5,83,0,0,381,382,5,98,0,0,
        382,102,1,0,0,0,383,384,5,84,0,0,384,385,5,101,0,0,385,104,1,0,0,
        0,386,387,5,73,0,0,387,106,1,0,0,0,388,389,5,88,0,0,389,390,5,101,
        0,0,390,108,1,0,0,0,391,392,5,67,0,0,392,393,5,115,0,0,393,110,1,
        0,0,0,394,395,5,66,0,0,395,396,5,97,0,0,396,112,1,0,0,0,397,398,
        5,76,0,0,398,399,5,97,0,0,399,114,1,0,0,0,400,401,5,67,0,0,401,402,
        5,101,0,0,402,116,1,0,0,0,403,404,5,80,0,0,404,405,5,114,0,0,405,
        118,1,0,0,0,406,407,5,78,0,0,407,408,5,100,0,0,408,120,1,0,0,0,409,
        410,5,80,0,0,410,411,5,109,0,0,411,122,1,0,0,0,412,413,5,83,0,0,
        413,414,5,109,0,0,414,124,1,0,0,0,415,416,5,69,0,0,416,417,5,117,
        0,0,417,126,1,0,0,0,418,419,5,71,0,0,419,420,5,100,0,0,420,128,1,
        0,0,0,421,422,5,84,0,0,422,423,5,98,0,0,423,130,1,0,0,0,424,425,
        5,68,0,0,425,426,5,121,0,0,426,132,1,0,0,0,427,428,5,72,0,0,428,
        429,5,111,0,0,429,134,1,0,0,0,430,431,5,69,0,0,431,432,5,114,0,0,
        432,136,1,0,0,0,433,434,5,84,0,0,434,435,5,109,0,0,435,138,1,0,0,
        0,436,437,5,89,0,0,437,438,5,98,0,0,438,140,1,0,0,0,439,440,5,76,
        0,0,440,441,5,117,0,0,441,142,1,0,0,0,442,443,5,72,0,0,443,444,5,
        102,0,0,444,144,1,0,0,0,445,446,5,84,0,0,446,447,5,97,0,0,447,146,
        1,0,0,0,448,449,5,87,0,0,449,148,1,0,0,0,450,451,5,82,0,0,451,452,
        5,101,0,0,452,150,1,0,0,0,453,454,5,79,0,0,454,455,5,115,0,0,455,
        152,1,0,0,0,456,457,5,73,0,0,457,458,5,114,0,0,458,154,1,0,0,0,459,
        460,5,80,0,0,460,461,5,116,0,0,461,156,1,0,0,0,462,463,5,65,0,0,
        463,464,5,117,0,0,464,158,1,0,0,0,465,466,5,72,0,0,466,467,5,103,
        0,0,467,160,1,0,0,0,468,469,5,84,0,0,469,470,5,108,0,0,470,162,1,
        0,0,0,471,472,5,80,0,0,472,473,5,98,0,0,473,164,1,0,0,0,474,475,
        5,66,0,0,475,476,5,105,0,0,476,166,1,0,0,0,477,478,5,80,0,0,478,
        479,5,111,0,0,479,168,1,0,0,0,480,481,5,65,0,0,481,482,5,116,0,0,
        482,170,1,0,0,0,483,484,5,82,0,0,484,485,5,110,0,0,485,172,1,0,0,
        0,486,487,5,70,0,0,487,488,5,114,0,0,488,174,1,0,0,0,489,490,5,82,
        0,0,490,491,5,97,0,0,491,176,1,0,0,0,492,493,5,65,0,0,493,494,5,
        99,0,0,494,178,1,0,0,0,495,496,5,84,0,0,496,497,5,104,0,0,497,180,
        1,0,0,0,498,499,5,80,0,0,499,500,5,97,0,0,500,182,1,0,0,0,501,502,
        5,85,0,0,502,184,1,0,0,0,503,504,5,78,0,0,504,505,5,112,0,0,505,
        186,1,0,0,0,506,507,5,80,0,0,507,508,5,117,0,0,508,188,1,0,0,0,509,
        510,5,65,0,0,510,511,5,109,0,0,511,190,1,0,0,0,512,513,5,67,0,0,
        513,514,5,109,0,0,514,192,1,0,0,0,515,516,5,66,0,0,516,517,5,107,
        0,0,517,194,1,0,0,0,518,519,5,67,0,0,519,520,5,102,0,0,520,196,1,
        0,0,0,521,522,5,69,0,0,522,523,5,115,0,0,523,198,1,0,0,0,524,525,
        5,70,0,0,525,526,5,109,0,0,526,200,1,0,0,0,527,528,5,77,0,0,528,
        529,5,100,0,0,529,202,1,0,0,0,530,531,5,78,0,0,531,532,5,111,0,0,
        532,204,1,0,0,0,533,534,5,76,0,0,534,535,5,114,0,0,535,206,1,0,0,
        0,536,537,5,82,0,0,537,538,5,102,0,0,538,208,1,0,0,0,539,540,5,68,
        0,0,540,541,5,98,0,0,541,210,1,0,0,0,542,543,5,83,0,0,543,544,5,
        103,0,0,544,212,1,0,0,0,545,546,5,66,0,0,546,547,5,104,0,0,547,214,
        1,0,0,0,548,549,5,72,0,0,549,550,5,115,0,0,550,216,1,0,0,0,551,552,
        5,77,0,0,552,553,5,116,0,0,553,218,1,0,0,0,554,555,5,68,0,0,555,
        556,5,115,0,0,556,220,1,0,0,0,557,558,5,82,0,0,558,559,5,103,0,0,
        559,222,1,0,0,0,560,561,5,67,0,0,561,562,5,110,0,0,562,224,1,0,0,
        0,563,564,5,78,0,0,564,565,5,104,0,0,565,226,1,0,0,0,566,567,5,70,
        0,0,567,568,5,108,0,0,568,228,1,0,0,0,569,570,5,77,0,0,570,571,5,
        99,0,0,571,230,1,0,0,0,572,573,5,76,0,0,573,574,5,118,0,0,574,232,
        1,0,0,0,575,576,5,84,0,0,576,577,5,115,0,0,577,234,1,0,0,0,578,579,
        5,79,0,0,579,580,5,103,0,0,580,236,1,0,0,0,581,582,7,0,0,0,582,238,
        1,0,0,0,583,585,7,1,0,0,584,586,7,2,0,0,585,584,1,0,0,0,586,587,
        1,0,0,0,587,585,1,0,0,0,587,588,1,0,0,0,588,240,1,0,0,0,2,0,587,
        0
    ]

class tucanLexer(Lexer):

    atn = ATNDeserializer().deserialize(serializedATN())

    decisionsToDFA = [ DFA(ds, i) for i, ds in enumerate(atn.decisionToState) ]

    T__0 = 1
    T__1 = 2
    T__2 = 3
    T__3 = 4
    T__4 = 5
    T__5 = 6
    T__6 = 7
    T__7 = 8
    T__8 = 9
    T__9 = 10
    T__10 = 11
    T__11 = 12
    T__12 = 13
    T__13 = 14
    T__14 = 15
    T__15 = 16
    T__16 = 17
    T__17 = 18
    T__18 = 19
    T__19 = 20
    T__20 = 21
    T__21 = 22
    T__22 = 23
    T__23 = 24
    T__24 = 25
    T__25 = 26
    T__26 = 27
    T__27 = 28
    T__28 = 29
    T__29 = 30
    T__30 = 31
    T__31 = 32
    T__32 = 33
    T__33 = 34
    T__34 = 35
    T__35 = 36
    T__36 = 37
    T__37 = 38
    T__38 = 39
    T__39 = 40
    T__40 = 41
    T__41 = 42
    T__42 = 43
    T__43 = 44
    T__44 = 45
    T__45 = 46
    T__46 = 47
    T__47 = 48
    T__48 = 49
    T__49 = 50
    T__50 = 51
    T__51 = 52
    T__52 = 53
    T__53 = 54
    T__54 = 55
    T__55 = 56
    T__56 = 57
    T__57 = 58
    T__58 = 59
    T__59 = 60
    T__60 = 61
    T__61 = 62
    T__62 = 63
    T__63 = 64
    T__64 = 65
    T__65 = 66
    T__66 = 67
    T__67 = 68
    T__68 = 69
    T__69 = 70
    T__70 = 71
    T__71 = 72
    T__72 = 73
    T__73 = 74
    T__74 = 75
    T__75 = 76
    T__76 = 77
    T__77 = 78
    T__78 = 79
    T__79 = 80
    T__80 = 81
    T__81 = 82
    T__82 = 83
    T__83 = 84
    T__84 = 85
    T__85 = 86
    T__86 = 87
    T__87 = 88
    T__88 = 89
    T__89 = 90
    T__90 = 91
    T__91 = 92
    T__92 = 93
    T__93 = 94
    T__94 = 95
    T__95 = 96
    T__96 = 97
    T__97 = 98
    T__98 = 99
    T__99 = 100
    T__100 = 101
    T__101 = 102
    T__102 = 103
    T__103 = 104
    T__104 = 105
    T__105 = 106
    T__106 = 107
    T__107 = 108
    T__108 = 109
    T__109 = 110
    T__110 = 111
    T__111 = 112
    T__112 = 113
    T__113 = 114
    T__114 = 115
    T__115 = 116
    T__116 = 117
    T__117 = 118
    TWO_TO_NINE = 119
    GREATER_THAN_NINE = 120

    channelNames = [ u"DEFAULT_TOKEN_CHANNEL", u"HIDDEN" ]

    modeNames = [ "DEFAULT_MODE" ]

    literalNames = [ "<INVALID>",
            "'H'", "'He'", "'Li'", "'Be'", "'B'", "'C'", "'N'", "'O'", "'F'", 
            "'Ne'", "'Na'", "'Mg'", "'Al'", "'Si'", "'P'", "'S'", "'Cl'", 
            "'Ar'", "'K'", "'Ca'", "'Sc'", "'Ti'", "'V'", "'Cr'", "'Mn'", 
            "'Fe'", "'Co'", "'Ni'", "'Cu'", "'Zn'", "'Ga'", "'Ge'", "'As'", 
            "'Se'", "'Br'", "'Kr'", "'Rb'", "'Sr'", "'Y'", "'Zr'", "'Nb'", 
            "'Mo'", "'Tc'", "'Ru'", "'Rh'", "'Pd'", "'Ag'", "'Cd'", "'In'", 
            "'Sn'", "'Sb'", "'Te'", "'I'", "'Xe'", "'Cs'", "'Ba'", "'La'", 
            "'Ce'", "'Pr'", "'Nd'", "'Pm'", "'Sm'", "'Eu'", "'Gd'", "'Tb'", 
            "'Dy'", "'Ho'", "'Er'", "'Tm'", "'Yb'", "'Lu'", "'Hf'", "'Ta'", 
            "'W'", "'Re'", "'Os'", "'Ir'", "'Pt'", "'Au'", "'Hg'", "'Tl'", 
            "'Pb'", "'Bi'", "'Po'", "'At'", "'Rn'", "'Fr'", "'Ra'", "'Ac'", 
            "'Th'", "'Pa'", "'U'", "'Np'", "'Pu'", "'Am'", "'Cm'", "'Bk'", 
            "'Cf'", "'Es'", "'Fm'", "'Md'", "'No'", "'Lr'", "'Rf'", "'Db'", 
            "'Sg'", "'Bh'", "'Hs'", "'Mt'", "'Ds'", "'Rg'", "'Cn'", "'Nh'", 
            "'Fl'", "'Mc'", "'Lv'", "'Ts'", "'Og'" ]

    symbolicNames = [ "<INVALID>",
            "TWO_TO_NINE", "GREATER_THAN_NINE" ]

    ruleNames = [ "T__0", "T__1", "T__2", "T__3", "T__4", "T__5", "T__6", 
                  "T__7", "T__8", "T__9", "T__10", "T__11", "T__12", "T__13", 
                  "T__14", "T__15", "T__16", "T__17", "T__18", "T__19", 
                  "T__20", "T__21", "T__22", "T__23", "T__24", "T__25", 
                  "T__26", "T__27", "T__28", "T__29", "T__30", "T__31", 
                  "T__32", "T__33", "T__34", "T__35", "T__36", "T__37", 
                  "T__38", "T__39", "T__40", "T__41", "T__42", "T__43", 
                  "T__44", "T__45", "T__46", "T__47", "T__48", "T__49", 
                  "T__50", "T__51", "T__52", "T__53", "T__54", "T__55", 
                  "T__56", "T__57", "T__58", "T__59", "T__60", "T__61", 
                  "T__62", "T__63", "T__64", "T__65", "T__66", "T__67", 
                  "T__68", "T__69", "T__70", "T__71", "T__72", "T__73", 
                  "T__74", "T__75", "T__76", "T__77", "T__78", "T__79", 
                  "T__80", "T__81", "T__82", "T__83", "T__84", "T__85", 
                  "T__86", "T__87", "T__88", "T__89", "T__90", "T__91", 
                  "T__92", "T__93", "T__94", "T__95", "T__96", "T__97", 
                  "T__98", "T__99", "T__100", "T__101", "T__102", "T__103", 
                  "T__104", "T__105", "T__106", "T__107", "T__108", "T__109", 
                  "T__110", "T__111", "T__112", "T__113", "T__114", "T__115", 
                  "T__116", "T__117", "TWO_TO_NINE", "GREATER_THAN_NINE" ]

    grammarFileName = "tucan.g4"

    def __init__(self, input=None, output:TextIO = sys.stdout):
        super().__init__(input, output)
        self.checkVersion("4.11.1")
        self._interp = LexerATNSimulator(self, self.atn, self.decisionsToDFA, PredictionContextCache())
        self._actions = None
        self._predicates = None


