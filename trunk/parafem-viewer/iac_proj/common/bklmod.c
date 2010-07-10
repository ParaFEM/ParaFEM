#include "iac_proj/common/gen.h"


int
BKL_update(OMobj_id BreakLinkCore_id, OMevent_mask event_mask, int seq_num)
{
   /***********************/
   /*  Declare variables  */
   /***********************/
   OMobj_id in_id, make_id, break_id, out_id;

   /***********************/
   /*  Get input values   */
   /***********************/
   in_id = OMfind_subobj(BreakLinkCore_id, OMstr_to_name("in"), OM_OBJ_RD);
   if (OMis_null_obj(in_id)) {
      ERRerror("BKL_update",1,ERR_ORIG, "Error searching for in.");
      return 0;
   }

   make_id = OMfind_subobj(BreakLinkCore_id, OMstr_to_name("make"), OM_OBJ_RD);
   if (OMis_null_obj(make_id)) {
      ERRerror("BKL_update",1,ERR_ORIG, "Error searching for make.");
      return 0;
   }

   break_id = OMfind_subobj(BreakLinkCore_id, OMstr_to_name("break"), OM_OBJ_RD);
   if (OMis_null_obj(break_id)) {
      ERRerror("BKL_update",1,ERR_ORIG, "Error searching for break.");
      return 0;
   }

   out_id = OMfind_subobj(BreakLinkCore_id, OMstr_to_name("out"), OM_OBJ_RW);
   if (OMis_null_obj(out_id)) {
      ERRerror("BKL_update",1,ERR_ORIG, "Error searching for out.");
      return 0;
   }



   /***********************/
   /* Function's Body     */
   /***********************/
   /* ERRerror("BKL_update",1,ERR_ORIG,"I'm in function: BKL_update generated from method: BreakLinkCore.BKL_update\n"); */


   /***********************/
   /*  Set output values  */
   /***********************/
   if (OMchanged(make_id, seq_num)) { /* make triggered */
      /*Connect out to in*/
      OMset_obj_ref(out_id, in_id, 0);
   }
   else if (OMchanged(break_id, seq_num)) { /* break triggered */
      /*Disconnect out from in*/
      OMset_obj_ref(out_id, OMnull_obj, 0);
   }

   return(1);
}

