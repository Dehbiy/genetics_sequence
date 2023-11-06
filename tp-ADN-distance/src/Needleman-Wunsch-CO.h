/**
 * \file Needleman-Wunsch-CO.h
 * \brief recursive cache oblivious implementation of Needleman-Wunsch global alignment algorithm that computes the distance
 * between two genetic sequences
 * \date 28/10/2023
 * \author Yakoub Dehbi (Ensimag, Grenoble-INP - University Grenoble-Alpes) yakoub.dehbi@grenoble-inp.org
 * \author Youssef Elaasri (Ensimag, Grenoble-INP - University Grenoble-Alpes) youssef.elaasri@grenoble-inp.org
 */


#include <stdlib.h>


long EditDistance_NW_Cache_Oblivious(char *A, size_t lengthA, char *B, size_t lengthB);