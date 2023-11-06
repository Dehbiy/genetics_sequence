#include "Needleman-Wunsch-recmemo.h"
#include <stdio.h>  
#include "Needleman-Wunsch-CO.h" 
#include <string.h> /* for strchr */
// #include <ctype.h> /* for toupper */

#include "characters_to_base.h" /* mapping from char to base */

/*****************************************************************************/
   
/* Context of the memoization : passed to all recursive calls */
/** \def NOT_YET_COMPUTED
 * \brief default value for memoization of minimal distance (defined as an impossible value for a distance, -1).
 */
#define NOT_YET_COMPUTED -1L 
#define K 120
#define ZERO 0


long EditDistance_NW_Cache_Oblivious_it(char *A, size_t lengthA, 
                                       char *B, size_t lengthB,
                                       int beginI, int endI,
                                       int beginJ, int endJ, 
                                       long* phi, long* ksi,
                                       long corner){

   //precPhi is the variable containing the value of the eastern element of the current one.
   //After each iteration it should be stocked in the in the array phi.
   long precPhi = 42;


    //The block defined between the lines beginJ and endJ both included
   for(int j =  endJ; j >= beginJ; j--){

        //Initializing precPhi is tricky
        //It should at first contain the corresponding element in ksi (the one one the same line), since it is the eastern element that actually belongs to the previous block
        //precPhi isn't relevant when the calculated element is on the eastern border.
      precPhi = ksi[j];


      //The block defined between the columns beginJ and endJ both included
      for(int i = endI; i >= beginI; i--){
         if(i == lengthA){
            // First we treat the case of the top right element (always being zero)
            if(j == lengthB){
                //update precPhi and ksi[lengthB]
               precPhi = ZERO;
               ksi[j] = ZERO;
            }

            //In case the treated element is on the eastern border, it becomes clear that precPhi is not used to calculate the next precPhi
            else{
                //phi[i] is the northen element, its new value is  actually ksi[j] that is being calculated here
               ksi[j] = 2 * isBase(B[j]) + phi[i];
               //updating phi; the case that should be updated is beginI. Notice that the value put in phi[beginI] here is not precPhi but ksi[j + 1].
               //In fact, the treated element doesn't have an eastern element to it.
               phi[j==endJ? 0 : beginI] = ksi[j + 1];
               precPhi = ksi[j];
            }

         }
            // Northern border
         else if(j == lengthB){
            //this one is very ituitive.
            //update phi
            phi[i + 1] = precPhi;
            //calculate the next precPhi 
            precPhi = 2 * isBase(A[i]) + precPhi;
         }

         else if(isBase(A[i]) == 0){
            //if Xi is not base the the current element is actually equal to its eastern brother
            phi[i + 1] = precPhi;
         }

         else if(isBase(B[j]) == 0){
            //if Xi is not base the the current element is actually equal to its northern brother
            phi[i + 1] = precPhi;
            precPhi = phi[i];
         }



         else{
            // sigma checks if Xi is equal to Yj
            int sigma =  isUnknownBase(A[i]) ?  SUBSTITUTION_UNKNOWN_COST 
                              : ( isSameBase(A[i], B[j]) ? 0 : SUBSTITUTION_COST ) 
                        ;

            //Tricky part! in order to calculate phi1 the program should access the north eastern element
            //In the general case it's phi[i+1]. Since it hasn't been updated, is still has the correct value
            //However, for the first element is the bottom left block, the north eastern eastern element cannot be found in either arrays phi and ksi.
            //that's way it's being transported via the argument corner          
            int phi1 = sigma + ((corner == -1) ? phi[i+1]: corner);
            //Once the corner in used, it gets set to the default value -1
            corner = -1;
            //Precphi is the eastern value
            int phi2 = 2 + precPhi;
            //phi[i] is the northern value
            int phi3 = 2 + phi[i];  
            //lastly we update phi
            phi[i + 1] = precPhi;
            //calculate precPhi as the smallest phi
            precPhi = (phi1 < phi2) ? ((phi1 < phi3) ? phi1 : phi3) : ((phi2 < phi3) ? phi2 : phi3);
         }
      }
        //once we are done with each line of the block, ksi[j] is update to the last element calculated.
      ksi[j] = precPhi;
        //lastly we update phi to hop to the next line
      phi[beginI] = precPhi;
   }

   return precPhi;
}








long EditDistance_NW_Cache_Oblivious_rec(char *A, size_t lengthA,
                                       char *B, size_t lengthB,
                                       int beginI, int endI,
                                       int beginJ, int endJ,
                                       long* phi, long* ksi, long corner){
    //if the block is small enough the compute it
   if (((endI - beginI  < K) && (endJ - beginJ  < K)) || (endJ == beginJ)  || (endI == beginI)){
      return EditDistance_NW_Cache_Oblivious_it(A, lengthA, B, lengthB, beginI, endI, beginJ, endJ, phi, ksi, corner);
      
   }
    //else devide it to four sections
   corner =  EditDistance_NW_Cache_Oblivious_rec(A, lengthA,
                                       B, lengthB,
                                       (beginI + endI)/2 + 1, endI, 
                                       (beginJ + endJ)/2 + 1, endJ,
                                       phi, ksi, corner);

   EditDistance_NW_Cache_Oblivious_rec(A, lengthA,
                                       B, lengthB,
                                       beginI, ((beginI + endI)/2),
                                       (beginJ + endJ)/2 + 1, endJ,
                                       phi, ksi, -1);
   
   EditDistance_NW_Cache_Oblivious_rec(A, lengthA,
                                       B, lengthB,
                                       (beginI + endI)/2 + 1, endI,
                                       beginJ, ((beginJ + endJ)/2),
                                       phi, ksi, -1);
   //we pass the result of the first one as the corner of the last one
   return EditDistance_NW_Cache_Oblivious_rec(A, lengthA,
                                       B,  lengthB, 
                                       beginI, ((beginI + endI)/2),
                                       beginJ, ((beginJ + endJ)/2),
                                       phi, ksi, corner);

}




long EditDistance_NW_Cache_Oblivious(char *A, size_t lengthA, char *B, size_t lengthB){

   _init_base_match() ;
   //A is the smalled sequence and B is the biggest
   if(lengthA < lengthB){
      long* phi = malloc(sizeof(long) * (lengthA+1));
      long* ksi = malloc(sizeof(long) * (lengthB+1));

      long result = EditDistance_NW_Cache_Oblivious_rec(A, lengthA,
                                          B, lengthB,
                                          ZERO, lengthA,
                                          ZERO, lengthB,
                                          phi, ksi, -1);
      free(phi); free(ksi);
      return result;


   }

   else{
      long* phi = malloc(sizeof(long) * (lengthB+1));
      long* ksi = malloc(sizeof(long) * (lengthA+1));

      long result = EditDistance_NW_Cache_Oblivious_rec(B, lengthB,
                                          A, lengthA,
                                          ZERO, lengthB,
                                          ZERO, lengthA,
                                          phi, ksi, -1);
      free(phi); free(ksi);
      return result;
      
   }

}