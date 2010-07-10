// INTERNATIONAL AVS CENTRE - WARRANTY DISCLAIMER
// Please read the file DISCLAIMER for conditions associated with this file.
// avs@iavsc.org, www.iavsc.org

flibrary LogicWidgets <compile_subs=0> {


   // Evaluates an input condition and outputs the evaluated condition.
   // If the condition evaluates to true then the condition true output
   // is triggered else the condition false output is triggered.
   macro IfTrue {
      ilink trigger;
      ilink in;

      IAC_PROJ.Common.CommonMods.TrueFalseCore TrueFalseCore {
         trigger => <-.trigger;
         state => in;
      };

      olink cond => TrueFalseCore.state_out;
      olink cond_true => TrueFalseCore.state_true;
      olink cond_false => TrueFalseCore.state_false;
   };

   // Evaluates an input condition and outputs the evaluated condition.
   // If the condition evaluates to false then the condition true output
   // is triggered else the condition false output is triggered.
   IfTrue IfFalse {
      TrueFalseCore {
         state => (!in);
      };
   };


   // If the input values are equal then the condition true output is
   // triggered else the condition false output is triggered.
   macro IfEqualTo {
      ilink in1;
      ilink in2;

      IAC_PROJ.Common.CommonMods.TrueFalseCore TrueFalseCore {
         state => (in1 == in2);
      };

      olink cond => TrueFalseCore.state_out;
      olink cond_true => TrueFalseCore.state_true;
      olink cond_false => TrueFalseCore.state_false;
   };

   // If the input values are not equal then the condition true output is
   // triggered else the condition false output is triggered.
   IfEqualTo IfNotEqualTo {
      TrueFalseCore {
         state => (in1 != in2);
      };
   };

   // If the first input value is less than the second then the condition
   // true output is triggered else the condition false output is triggered.
   IfEqualTo IfLessThan {
      TrueFalseCore {
         state => (in1 < in2);
      };
   };

   // If the first input value is greater than the second then the condition
   // true output is triggered else the condition false output is triggered.
   IfEqualTo IfGreaterThan {
      TrueFalseCore {
         state => (in1 > in2);
      };
   };

   // If the first input value is less than or equal to the second then the condition
   // true output is triggered else the condition false output is triggered.
   IfEqualTo IfLessThanEqual {
      TrueFalseCore {
         state => (in1 <= in2);
      };
   };

   // If the first input value is greater than or equal to the second then the condition
   // true output is triggered else the condition false output is triggered.
   IfEqualTo IfGreaterThanEqual {
      TrueFalseCore {
         state => (in1 >= in2);
      };
   };

   // If the first input value logical ANDed with the second value evaluates to true then
   // the condition true output is triggered else the condition false output is triggered.
   IfEqualTo IfLogicalAND {
      TrueFalseCore {
         state => (in1 && in2);
      };
   };

   // If the first input value logical ORed with the second value evaluates to true then
   // the condition true output is triggered else the condition false output is triggered.
   IfEqualTo IfLogicalOR {
      TrueFalseCore {
         state => (in1 || in2);
      };
   };

   // If the first input value bitwise ANDed with the second value evaluates to true then
   // the condition true output is triggered else the condition false output is triggered.
   IfEqualTo IfBitwiseAND {
      TrueFalseCore {
         state => (in1 & in2);
      };
   };

   // If the first input value bitwise ORed with the second value evaluates to true then
   // the condition true output is triggered else the condition false output is triggered.
   IfEqualTo IfBitwiseOR {
      TrueFalseCore {
         state => (in1 | in2);
      };
   };



   // Toggle between 0 and 1 each time a trigger input is recieved
   group Toggle {
      ilink trigger;

      IAC_PROJ.Common.CommonMods.TrueFalseCore TrueFalseCore {
         trigger => <-.trigger;
         state => (!<-.state);
         state_out => <-.state;
      };
      int state = 0;

      olink out => state;
      olink not_out => (!state);
      olink value_true => TrueFalseCore.state_true;
      olink value_false => TrueFalseCore.state_false;
   };


   // Sets an output value to true or false depending upon on and off trigger inputs
   macro Switch {
      ilink on;
      ilink off;

      GMOD.copy_on_change set_on {
         input = 1;
         trigger+IPort2 => <-.on;
         output => <-.out;
      };
      GMOD.copy_on_change set_off {
         input = 0;
         trigger+IPort2 => <-.off;
         output => <-.out;
      };

      int+OPort2 out = 0;
   };



   // Executes a v command string if the input condition evaluates to true
   group parse_v_on_true {
      ilink in;
      IAC_PROJ.Common.CommonMods.TrueFalseCore TrueFalseCore {
         state => <-.in;
      };
      GMOD.parse_v parse_v {
         v_commands+IPort3;
         trigger => <-.TrueFalseCore.state_true;
         active+IPort3;
         on_inst+IPort3;
         relative+IPort3;
      };
   };

   // Executes a v command string if the input condition evaluates to false
   parse_v_on_true parse_v_on_false {
      parse_v {
         trigger => <-.TrueFalseCore.state_false;
      };
   };

   // Executes a v command string if the input values are equal
   macro parse_v_on_equal {
      ilink in1;
      ilink in2;

      ilink v_commands;
      ilink active;
      ilink on_inst;
      ilink relative;

      parse_v_on_true parse_v_on_true {
         TrueFalseCore {
            state => (in1 == in2);
         };
         parse_v {
            v_commands => <-.<-.v_commands;
            active => <-.<-.active;
            on_inst => <-.<-.on_inst;
            relative => <-.<-.relative;
         };
      };

   };

   // Executes a v command string if the input values are not equal
   parse_v_on_equal parse_v_on_not_equal {
      parse_v_on_true {
         TrueFalseCore {
            state => (in1 != in2);
         };
      };
   };

   // Executes a v command string if the first input value is less than the second
   parse_v_on_equal parse_v_on_less_than {
      parse_v_on_true {
         TrueFalseCore {
            state => (in1 < in2);
         };
      };
   };

   // Executes a v command string if the first input value is greater than the second
   parse_v_on_equal parse_v_on_greater_than {
      parse_v_on_true {
         TrueFalseCore {
            state => (in1 > in2);
         };
      };
   };


};

