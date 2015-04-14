#include <stdio.h>
#include <string.h>

struct {
   float AIST[500];
   float PERC[500];
   float ABUNOR[500];
   int NUCADR[110][2];
   int IZMAX, L110;
} extern abunuc_;

struct {
   int IRES[600];
   int ITYPE[600];
   double BR[600];
   int IMODE[600][5];
   int NLP;
} extern atabdata_;

extern void initis_();
extern void setdky_();

int main(int n_arg, char *arg[])
{
   FILE *f;
   f = fopen("tabnuc.h", "w");
   if (f == 0) {
      printf("HADGEN DATA PARSER FATAL: cannot open tabnuc.h for output\n");
      return 1;
   } else {
      initis_();
      int i;
      for(i=0; i<abunuc_.IZMAX; i++) {
         int AD1, AD2;
         AD1 = abunuc_.NUCADR[i][0];
         AD2 = abunuc_.NUCADR[i][1];
         int j;
         for(j=AD1; j<=AD2; j++) {
            fprintf(f, "   _hdg_add_isotope(%i., %3.2f, %3.4f);\n", 
               i+1, abunuc_.AIST[j-1], abunuc_.PERC[j-1]);
         }
         fprintf(f, "   _hdg_done();\n");
      } 
      fclose(f);
   }

   f = fopen("atab.h", "w");
   if (f == 0) {
      printf("HADGEN DATA PARSER FATAL: cannot open atab.h for output\n");
      return 1;
   } else {
      setdky_();
      int i;
      for(i=0; i<atabdata_.NLP; i++) {
         fprintf(f, "   _hdg_setdky(%6i, %3i, %6i, %6i, %6i, %6i, %6i, %15e);\n",
            atabdata_.IRES[i], atabdata_.ITYPE[i], 
            atabdata_.IMODE[i][0], atabdata_.IMODE[i][1],
            atabdata_.IMODE[i][2], atabdata_.IMODE[i][3],
            atabdata_.IMODE[i][4], atabdata_.BR[i]);
      }
      fclose(f);
   }
   return 0;
}
