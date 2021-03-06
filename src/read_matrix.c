#include "config.h"

#include <stdio.h>
#include <stdlib.h>
#include "gadem.h"

#ifdef HAVE_MALLOC_H
#include <malloc.h>
#endif

#define  MAX_MATRIX_LENGTH 100

double **alloc_double_double(int ,int );

double **read_initial_pwm(int *matrixLen,char *fileName) {

   register int m,n;
   int numCol,numRow;
   double colSum;
   double **pwm;
   FILE *fp;
   int checkfscanf;

   pwm=alloc_double_double(MAX_MATRIX_LENGTH,4);

   fp=fopen(fileName,"r");
   if (!fp) { error(fileName);}
   /*if (!fp) { perror(fileName); exit(0); }*/

   checkfscanf = fscanf(fp,"%d %d",&numRow,&numCol);
   if (numRow !=4) {
      error("\n\nError: please check documentation for input format.\n");
      /*Rprintf("\n\nError: please check documentation for input format.\n"); exit(0); */
   }
   /*if (numRow !=4) { Rprintf("\n\nError: please check documentation for input format.\n"); exit(0); }*/
   if (numCol<5)  Rprintf("\n\nWarning: motif length is %d SHORT\n",numCol); 

   for (m=0; m<4; m++) {
      for (n=0; n<numCol; n++) {
         checkfscanf = fscanf(fp,"%lf",&pwm[n][m]);
         if (pwm[n][m]<0) { 
            error("\n\nError: elements in PWM must be positive. Please see examples on gapwm website\n"); 
            /*Rprintf("\n\nError: elements in PWM must be positive. Please see examples on gapwm website\n"); */
            /*exit(0); */
         } 
      }
   }
   fclose(fp);

   Rprintf("\nInitial PWM:\n");
   for (m=0; m<4; m++) {
      for (n=0; n<numCol; n++) {
        if (n<numCol-1) Rprintf("%5.3f\t",pwm[n][m]);
        else            Rprintf("%5.3f\n",pwm[n][m]); 
      }
   }

   for (n=0; n<numCol; n++) {
      colSum=0; for (m=0; m<4; m++) colSum  +=pwm[n][m];
      for (m=0; m<4; m++) pwm[n][m] /=colSum;
   }
   *matrixLen=numCol;

   return (pwm); 
}
