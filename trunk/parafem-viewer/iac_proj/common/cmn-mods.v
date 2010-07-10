// INTERNATIONAL AVS CENTRE - WARRANTY DISCLAIMER
// Please read the file DISCLAIMER for conditions associated with this file.
// avs@iavsc.org, www.iavsc.org

flibrary CommonMods  < build_dir="iac_proj/common",
                       out_src_file="gen.c",
                       out_hdr_file="gen.h" >
{

   module TrueFalseCore <src_file="tfmod.c"> {
      int+IPort2 trigger = 0;
      int+IPort2 state = 0;

      omethod+req TF_update(
         .state+read+notify+req,
         .trigger+notify,
         .state_out+nonotify+write,
         .state_true+nonotify+write,
         .state_false+nonotify+write
      ) = "TF_update";

      boolean+OPort2 state_out;
      boolean+OPort2 state_true;
      boolean+OPort2 state_false;
   };


   module BreakLinkCore <src_file="bklmod.c"> {
      ilink in;
      int+IPort2 make;
      int+IPort2 break;

      omethod+req BKL_update(
         .in+read+notify+req,
         .make+read+notify,
         .break+read+notify,
         .out+nonotify+write
      ) = "BKL_update";

      olink out;
   };



   // Common Parameter Block for the RandomNumCore module
   // Enables easy connection of separate components
   group+OPort RandomNumParams {
      int+Port2   num_vals;   // Number of random number that should be generated
      float+Port2 min_val;    // Minimum value of random numbers
      float+Port2 max_val;    // Maximim value of random numbers
      int rseed;              // Seed value for random number generator
   };


   // Low-level module.
   // Acts as a wrapper for the C code that generates random numbers
   module RandomNumCore <src_file="randmod.c"> {

      // Reference to external parameter block
      RandomNumParams+IPort2 &RandomNumParams;

      // Links into external parameter block
      int num_vals => RandomNumParams.num_vals;
      float min_val => RandomNumParams.min_val;
      float max_val => RandomNumParams.max_val;
      int rseed => RandomNumParams.rseed;

      // Defines name of C function and specifies how it
      // interacts with the module parameters
      omethod+req+notify_inst RAND_update(
         .num_vals+read+notify+req,
         .min_val+read+notify+req,
         .max_val+read+notify+req,
         .rseed+read+notify,
         .out_vals+write
      ) = "RAND_update";

      // Output array containing generated random numbers
      float+OPort2 out_vals[num_vals];
   };

};

