#include "iac_proj/common/gen.h"
#include <stdlib.h>

#define ERROR_RETURN(MESS) {\
                ERRerror("RandomNumCore",1,ERR_ORIG,MESS);\
                if (out_vals!=NULL) ARRfree((float *)out_vals);\
                return 0;\
                }


int
RAND_update(OMobj_id RandomNumCore_id, OMevent_mask event_mask, int seq_num)
{
   /***********************/
   /*  Declare variables  */
   /***********************/
   int  num_vals;
   int  rseed;
   double  min_val;
   double  max_val;
   int  out_vals_size = 0;
   float *out_vals = NULL; 

   int i;
   float diff;


   /***********************/
   /*  Get input values   */
   /***********************/
   /* Get num_vals's value */ 
   if (OMget_name_int_val(RandomNumCore_id, OMstr_to_name("num_vals"), &num_vals) != OM_STAT_SUCCESS) {
      ERROR_RETURN("unable to get number of values");
   }

   if (num_vals<1) {
      /* No values to be generated */
      return 0;
   }
   
   /* Get rseed's value */ 
   if (OMget_name_int_val(RandomNumCore_id, OMstr_to_name("rseed"), &rseed) == OM_STAT_SUCCESS) {
      srand( (unsigned)rseed );
   }

   /* Get min_val's value */
   if (OMget_name_real_val(RandomNumCore_id, OMstr_to_name("min_val"), &min_val) != OM_STAT_SUCCESS) {
      ERROR_RETURN("unable to get minimum value");
   }

   /* Get max_val's value */
   if (OMget_name_real_val(RandomNumCore_id, OMstr_to_name("max_val"), &max_val) != OM_STAT_SUCCESS) {
      ERROR_RETURN("unable to get maximum value");
   }

   /* If min is greater than max then swap values */
   if (min_val > max_val) {
      double temp = max_val;
      max_val = min_val;
      min_val = temp;
   }

   /* Pre-calculate difference between min and max */
   diff = max_val - min_val;


   /* Get pointer to out_vals array */
   out_vals = (float *)OMret_name_array_ptr(RandomNumCore_id,
                                            OMstr_to_name("out_vals"), OM_GET_ARRAY_WR,
                                            &out_vals_size, NULL);
   if ((out_vals==NULL) || (out_vals_size!=num_vals)) {
       ERROR_RETURN("unable to access output array");
   }



   /***********************/
   /*  Set output values  */
   /***********************/

   for (i=0; i<num_vals; i++) {
       out_vals[i] = min_val + diff * ((float)rand()/(float)RAND_MAX);
   }


   /* Free allocated arrays */

   if (out_vals != NULL) 
      ARRfree((float *)out_vals);

   return(1);
}

