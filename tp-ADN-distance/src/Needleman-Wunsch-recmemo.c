/**
 * \file Needleman-Wunsch-recmemo.c
 * \brief recursive implementation with memoization of Needleman-Wunsch global alignment algorithm that computes the distance between two genetic sequences 
 * \version 0.1
 * \date 03/10/2022 
 * \author Jean-Louis Roch (Ensimag, Grenoble-INP - University Grenoble-Alpes) jean-louis.roch@grenoble-inp.fr
 *
 * Documentation: see Needleman-Wunsch-recmemo.h
 * Costs of basic base opertaions (SUBSTITUTION_COST, SUBSTITUTION_UNKNOWN_COST, INSERTION_COST) are
 * defined in Needleman-Wunsch-recmemo.h
 */


#include "Needleman-Wunsch-recmemo.h"
#include <stdio.h>  
#include <stdlib.h> 
#include <string.h> /* for strchr */
// #include <ctype.h> /* for toupper */

#include "characters_to_base.h" /* mapping from char to base */

/*****************************************************************************/
   
/* Context of the memoization : passed to all recursive calls */
/** \def NOT_YET_COMPUTED
 * \brief default value for memoization of minimal distance (defined as an impossible value for a distance, -1).
 */
#define NOT_YET_COMPUTED -1L 
#define K 4
#define ZERO 4


/** \struct NW_MemoContext
 * \brief data for memoization of recursive Needleman-Wunsch algorithm 
*/
struct NW_MemoContext 
{
    char *X ; /*!< the longest genetic sequences */
    char *Y ; /*!< the shortest genetic sequences */
    size_t M; /*!< length of X */
    size_t N; /*!< length of Y,  N <= M */
    long **memo; /*!< memoization table to store memo[0..M][0..N] (including stopping conditions phi(M,j) and phi(i,N) */
} ;

/*
 *  static long EditDistance_NW_RecMemo(struct NW_MemoContext *c, size_t i, size_t j) 
 * \brief  EditDistance_NW_RecMemo :  Private (static)  recursive function with memoization \
 * direct implementation of Needleman-Wursch extended to manage FASTA sequences (cf TP description)
 * \param c : data passed for recursive calls that includes the memoization array 
 * \param i : starting position of the left sequence :  c->X[ i .. c->M ] 
 * \param j : starting position of the right sequence :  c->Y[ j .. c->N ] 
 */ 
static long EditDistance_NW_RecMemo(struct NW_MemoContext *c, size_t i, size_t j) 
/* compute and returns phi(i,j) using data in c -allocated and initialized by EditDistance_NW_Rec */
{
   if (c->memo[i][j] == NOT_YET_COMPUTED)
   {  
      long res ;
      char Xi = c->X[i] ;
      char Yj = c->Y[j] ;
      if (i == c->M) /* Reach end of X */
      {  if (j == c->N) res = 0;  /* Reach end of Y too */
         else res = (isBase(Yj) ? INSERTION_COST : 0) + EditDistance_NW_RecMemo(c, i, j+1) ;
      }
      else if (j == c->N) /* Reach end of Y but not end of X */
      {  res = (isBase(Xi) ? INSERTION_COST : 0) + EditDistance_NW_RecMemo(c, i+1, j) ;
      }
      else if (! isBase(Xi))  /* skip ccharacter in Xi that is not a base */
      {  ManageBaseError( Xi ) ;
         res = EditDistance_NW_RecMemo(c, i+1, j) ;
      }
      else if (! isBase(Yj))  /* skip ccharacter in Yj that is not a base */
      {  ManageBaseError( Yj ) ;
         res = EditDistance_NW_RecMemo(c, i, j+1) ;
      }
      else  
      {  /* Note that stopping conditions (i==M) and (j==N) are already stored in c->memo (cf EditDistance_NW_Rec) */ 
         long min = /* initialization  with cas 1*/
                   ( isUnknownBase(Xi) ?  SUBSTITUTION_UNKNOWN_COST 
                          : ( isSameBase(Xi, Yj) ? 0 : SUBSTITUTION_COST ) 
                   )
                   + EditDistance_NW_RecMemo(c, i+1, j+1) ; 
         { long cas2 = INSERTION_COST + EditDistance_NW_RecMemo(c, i+1, j) ;      
           if (cas2 < min) min = cas2 ;
         }
         { long cas3 = INSERTION_COST + EditDistance_NW_RecMemo(c, i, j+1) ;      
           if (cas3 < min) min = cas3 ; 
         }
         res = min ;
      }
       c->memo[i][j] = res ;
   }
   return c->memo[i][j] ;
}

/* EditDistance_NW_Rec :  is the main function to call, cf .h for specification 
 * It allocates and initailizes data (NW_MemoContext) for memoization and call the 
 * recursivefunction EditDistance_NW_RecMemo 
 * See .h file for documentation
 */
long EditDistance_NW_Rec(char* A, size_t lengthA, char* B, size_t lengthB)
{
   _init_base_match() ;
   struct NW_MemoContext ctx;
   if (lengthA >= lengthB) /* X is the longest sequence, Y the shortest */
   {  ctx.X = A ;
      ctx.M = lengthA ;
      ctx.Y = B ;
      ctx.N = lengthB ;
   }
   else
   {  ctx.X = B ;
      ctx.M = lengthB ;
      ctx.Y = A ;
      ctx.N = lengthA ;
   }
   size_t M = ctx.M ;
   size_t N = ctx.N ;
   {  /* Allocation and initialization of ctx.memo to NOT_YET_COMPUTED*/
      /* Note: memo is of size (N+1)*(M+1) but is stored as (M+1) distinct arrays each with (N+1) continuous elements 
       * It would have been possible to allocate only one big array memezone of (M+1)*(N+1) elements 
       * and then memo as an array of (M+1) pointers, the memo[i] being the address of memzone[i*(N+1)].
       */ 
      ctx.memo = (long **) malloc ( (M+1) * sizeof(long *)) ;
      if (ctx.memo == NULL) { perror("EditDistance_NW_Rec: malloc of ctx_memo" ); exit(EXIT_FAILURE); }
      for (int i=0; i <= M; ++i) 
      {  ctx.memo[i] = (long*) malloc( (N+1) * sizeof(long));
         if (ctx.memo[i] == NULL) { perror("EditDistance_NW_Rec: malloc of ctx_memo[i]" ); exit(EXIT_FAILURE); }
         for (int j=0; j<=N; ++j) ctx.memo[i][j] = NOT_YET_COMPUTED ;
      }
   }    
   
   /* Compute phi(0,0) = ctx.memo[0][0] by calling the recursive function EditDistance_NW_RecMemo */
   long res = EditDistance_NW_RecMemo( &ctx, 0, 0 ) ;
    
   { /* Deallocation of ctx.memo */
      for (int i=0; i <= M; ++i) free( ctx.memo[i] ) ;
      free( ctx.memo ) ;
   }
   return res ;
}

long EditDistance_NW_It(char* A, size_t lengthA, char* B, size_t lengthB){
   _init_base_match() ;
   if (lengthB < lengthA){
      char * C = A;
      A = B; B = C;
      size_t size = lengthA;
      lengthA = lengthB; lengthB = size;
   }
   long phi[lengthA+1];
   int j = lengthB;
   int i = lengthA;

   phi[i] = 0;

   i--;

   for (i; i > -1; i--){
      phi[i] = 2 * (isBase(A[i])) + phi[i+1];
   }
   long prec = phi[0];
   j--;

   for(j; j > -1; j--){
      phi[0] = prec;
      prec = 2 * (isBase(B[j])) + phi[lengthA];
      for(i = lengthA - 1; i > -1; i--){
         if(isBase(B[j]) == 0){
            phi[i + 1] = prec;
            prec = phi[i];
         }

         else if( isBase(A[i]) == 0){
            phi[i + 1] = prec;            
         }

         else{
            int sigma =  isUnknownBase(A[i]) ?  SUBSTITUTION_UNKNOWN_COST 
                          : ( isSameBase(A[i], B[j]) ? 0 : SUBSTITUTION_COST ) 
                   ;
            
            int phi1 = sigma + phi[i+1];
            int phi2 = 2 + prec;  
            int phi3 = 2 + phi[i];  
            phi[i + 1] = prec;

            prec = (phi1 < phi2) ? ((phi1 < phi3) ? phi1 : phi3) : ((phi2 < phi3) ? phi2 : phi3);
         }
      } 
   }
   return prec;
}
long EditDistance_NW_Cache_Aware(char *A, size_t lengthA, char *B, size_t lengthB){
   _init_base_match();

   // Ensure that sequenceA (A) is the longer sequence.
   if (lengthB < lengthA){
      char * temp = A;
      A = B;
      B = temp;
      size_t tempLength = lengthA;
      lengthA = lengthB;
      lengthB = tempLength;
   }
   
   // Calculate the number of iterations based on K or the length of sequenceB.
   int realK = K < lengthB ? K : lengthB;

   // Initialize arrays for phi and ksi.
   long phi[lengthA + 1];
   long ksi[realK + 1];

   int i = lengthA;
   int j = realK - 1;
   phi[i] = 0;
   int precKsi = 0;

   i--;

   // Calculate phi values for sequenceA.
   for (i; i > -1; i--){
      phi[i] = 2 * (isBase(A[i])) + phi[i + 1];
   }

   int ksiK = 0;
   long precPhi = phi[(int) lengthA - realK > 0 ? lengthA - realK : 0];

   // Iterate through sequenceB with a step size of realK.
   for(int BLength = lengthB; BLength > 0; BLength -= realK){
      ksi[realK] = ksiK;

      // Calculate ksi values for sequenceB.
      for (j = realK - 1; j > -1; j--){
         ksi[j] = 2 * (isBase(B[j])) + ksi[j + 1];
      }
      ksiK = ksi[realK - (int) lengthB > 0 ? realK - lengthB : 0];

      // Iterate through sequenceA with a step size of realK.
      for(int ALength = lengthA; ALength > 0; ALength -= realK){
         int ksiKPlus = phi[ALength - realK > 0 ? ALength - realK : 0];

         for(j = BLength; (j > BLength - realK && j > 0); j--){

            if(j != BLength) phi[ALength - realK > 0 ? ALength - realK : 0] = precPhi;

            // Handle cases where characters are not bases.
            if(isBase(B[j-1]) == 0 ){
               phi[ALength] = ksi[realK - BLength + j - 1];
               precPhi = phi[ALength-1];
            }
            else if( isBase(A[ALength-1]) == 0){
               phi[ALength] = ksi[realK - BLength + j - 1];
               precPhi = ksi[realK - BLength + j - 1];            
            }
            else{
               // Calculate costs for substitution and insertion/deletion.
               int sigma =  isUnknownBase(A[ALength-1]) ?  SUBSTITUTION_UNKNOWN_COST 
                          : ( isSameBase(A[ALength-1], B[j-1]) ? 0 : SUBSTITUTION_COST ) 
                  ;
            
               int phi1 = sigma + precKsi;
               int phi2 = 2 + ksi[realK - BLength + j - 1];  
               int phi3 = 2 + phi[ALength-1]; 
               precPhi = (phi1 < phi2) ? ((phi1 < phi3) ? phi1 : phi3) : ((phi2 < phi3) ? phi2 : phi3);
               phi[ALength] = ksi[realK - BLength + j - 1];
            }

            if (j == BLength) phi[ALength] = ksi[0];
            
            // Iterate through sequenceA with a step size of realK.
            for(i = ALength-1; (i > ALength - realK && i >0); i--){
               if(isBase(B[j-1]) == 0){
                  phi[i] = precPhi;
                  precPhi = phi[i-1];
               }
               else if( isBase(A[i-1]) == 0){
                  phi[i] = precPhi;            
               }
               else{
                  // Calculate costs for substitution and insertion/deletion.
                  int sigma =  isUnknownBase(A[i-1]) ?  SUBSTITUTION_UNKNOWN_COST 
                              : ( isSameBase(A[i-1], B[j-1]) ? 0 : SUBSTITUTION_COST ) 
                        ;
                  int phi1 = sigma + phi[i];
                  int phi2 = 2 + precPhi;  
                  int phi3 = 2 + phi[i-1];  
                  phi[i] = precPhi;
                  precPhi = (phi1 < phi2) ? ((phi1 < phi3) ? phi1 : phi3) : ((phi2 < phi3) ? phi2 : phi3);
               }
            }
            precKsi = ksi[realK - BLength + j - 1];
            ksi[realK - BLength + j - 1] = precPhi;
         }
         ksi[realK] = ksiKPlus;
         phi[ALength - realK > 0 ? ALength - realK : 0] = precPhi;
         precKsi = precPhi;
      }
      precKsi = ksiK;
   }

   return precPhi;
}



void EditDistance_NW_Cache_Oblivious_it(char *A, size_t lengthA, 
                                       char *B, size_t lengthB,
                                       int beginI, int endI,
                                       int beginJ, int endJ, 
                                       long* phi, long* ksi){
   
   for(int j =  endJ; j >= beginJ; j--){
      for(int i = endI; i >= beginI; i--){

      }
   }

}








void EditDistance_NW_Cache_Oblivious_rec(char *A, size_t lengthA,
                                       char *B, size_t lengthB,
                                       int beginI, int endI,
                                       int beginJ, int endJ,
                                       long* phi, long* ksi){
   if ((lengthA < K) && (lengthB < K)){
      EditDistance_NW_Cache_Oblivious_it(A, lengthA, B, lengthB, beginI, endI, beginJ, endJ, phi, ksi);
      return;
   }

   EditDistance_NW_Cache_Oblivious_rec(A, lengthA,
                                       B, lengthB,
                                       beginI + lengthA/2, endI, 
                                       beginJ + lengthB/2, endJ,
                                       phi, ksi);

   EditDistance_NW_Cache_Oblivious_rec(A, lengthA,
                                       B, lengthB,
                                       beginI, beginI + lengthA/2 - 1,
                                       beginJ + lengthB/2, endJ,
                                       phi, ksi);
   
   EditDistance_NW_Cache_Oblivious_rec(A, lengthA,
                                       B, lengthB,
                                       beginI + lengthA/2, endI,
                                       beginJ, beginJ + lengthB/2 - 1,
                                       phi, ksi);
   
   EditDistance_NW_Cache_Oblivious_rec(A, lengthA,
                                       B,  lengthB, 
                                       beginI, beginI + lengthA/2 - 1,
                                       beginJ, beginJ + lengthB/2 - 1,
                                       phi, ksi);

}




long EditDistance_NW_Cache_Oblivious(char *A, size_t lengthA, char *B, size_t lengthB){

   _init_base_match() ;
   
   if(lengthA < lengthB){
      long phi[lengthA];
      long ksi[lengthB];

      EditDistance_NW_Cache_Oblivious_rec(A, lengthA,
                                          B, lengthB,
                                          ZERO, lengthA,
                                          ZERO, lengthB,
                                          phi, ksi);

      return phi[0];

   }

   else{
      long phi[lengthB];
      long ksi[lengthA];

      EditDistance_NW_Cache_Oblivious_rec(B, lengthB,
                                          A, lengthA,
                                          ZERO, lengthB,
                                          ZERO, lengthA,
                                          phi, ksi);
      
      return phi[0];

   }

}