#include "iac_proj/common/gen.h"


int
TF_update(OMobj_id TrueFalseCore_id, OMevent_mask event_mask, int seq_num)
{
   /***********************/
   /*  Declare variables  */
   /***********************/
   int  state;

   /***********************/
   /*  Get input values   */
   /***********************/
   /* Get state's value */ 
   if (OMget_name_int_val(TrueFalseCore_id, OMstr_to_name("state"), &state) != OM_STAT_SUCCESS) 
      state = 0;


   /***********************/
   /* Function's Body     */
   /***********************/
   /* ERRerror("TF_update",1,ERR_ORIG,"I'm in function: TF_update generated from method: TrueFalseCore.TF_update\n"); */


   /***********************/
   /*  Set output values  */
   /***********************/
   if (state == 0) {
      OMset_name_int_val(TrueFalseCore_id, OMstr_to_name("state_out"), 0);
      OMset_name_int_val(TrueFalseCore_id, OMstr_to_name("state_false"), 1);
   } else {
      OMset_name_int_val(TrueFalseCore_id, OMstr_to_name("state_out"), 1);
      OMset_name_int_val(TrueFalseCore_id, OMstr_to_name("state_true"), 1);
   }


   return(1);
}

