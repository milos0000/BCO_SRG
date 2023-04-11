/*******************************************************************
 *                                                                 *
 *      FILE:      normal.c                                        *
 *                                                                 *
 *      FUNCTION:  normalization function to convert randomly      *
 *                 generated number into an integer within         *
 *                 given range                                     *
 *                                                                 *
 *      ARGUMENT:  r    - random number (int)                      *
 *                 a    - left border of the interval (int)        *
 *                 b    - right border of the interval (int)       *
 *                                                                 *
 *      RETURN:    aux - normalized random value                   *
 *                                                                 *
 *      AUTHOR: Tatjana Davidovic                                  *
 *              Mathematical Institute                             *
 *              Beograd, 2008-9.                                   *
 *******************************************************************/
#include <cstdio>
#include <cstdlib>



int normal(int r, int a, int b)
{

   int aux; 
   
   double c;
      
    c = (double) r/(double) RAND_MAX;
    c *= (double) (b - a + 1);
               
    aux = (int) c;
    aux += a;
                        
    return (aux);
}


