// INTERNATIONAL AVS CENTRE - WARRANTY DISCLAIMER
// Please read the file DISCLAIMER for conditions associated with this file.
// avs@iavsc.org, www.iavsc.org

flibrary NetworkWidgets <compile_subs=0> {

   // Takes a trigger input, and increments a output value
   group Increment {
      ilink in;
      ilink reset;

      int start_value = -1;

      GMOD.copy_on_change copy_on_change {
         input => <-.value;
         trigger<NEportLevels={2,0}> => <-.in;
         output = <-.start_value;
      };

      GMOD.parse_v parse_v {
         v_commands => "copy_on_change.output = <-.start_value;";
         trigger => <-.reset;
         relative => (<-.<-);
      };

      int value<NEportLevels={1,2}> => (copy_on_change.output + 1);
   };


   // Takes a trigger input, and decrements a output value
   Increment Decrement {
      start_value = 1;
      value => (copy_on_change.output - 1);
   };


   // Uses an Increment widget and a Decrement widget to create a
   // widget that can step a value up and down in integer sequence.
   macro Stepper {
      ilink up;
      ilink down;
      ilink reset;
      int start<NEportLevels={2,1}> = 0;

      Decrement decrement {
         in => <-.down;
         reset => <-.reset;
      };
      Increment increment {
         in => <-.up;
         reset => <-.reset;
      };

      int incVal<NEportLevels=1> =>.increment.value;
      int decVal<NEportLevels=1> =>.decrement.value;
    
      int value<NEportLevels={1,2}> => ((start+incVal)+decVal);
   };



   // Detects which of 2 integers is triggered last
   // and passes the value of that integer through to the output
   macro last_of_2 {
      int out<NEportLevels={1,2}>;

      int in1 <NEportLevels={2,1}>;
      int in2 <NEportLevels={2,1}>;

      GMOD.copy_on_change set1 {
            trigger => <-.in1;
            input =>   <-.in1;
            output<NEportLevels={2,2}> => <-.out;
      };
      GMOD.copy_on_change set2 {
            trigger => <-.in2;
            input =>   <-.in2;
            output<NEportLevels={2,2}> => <-.out;
      };
   };

   // Detects which of 3 integers is triggered last
   // and passes the value of that integer through to the output
   last_of_2 last_of_3 {
      int in3 <NEportLevels={2,1}>;

      GMOD.copy_on_change set3 {
            trigger => <-.in3;
            input =>   <-.in3;
            output<NEportLevels={2,2}> => <-.out;
      };
   };


   // Detects which of an array of integers is triggered last
   // and passes the index number and value of that integer through to the output
   // NOTE: Does not currently work correctly.
/*
   macro last_of_array {
      int in<NEportLevels={2,1}>[];
      int size<NEportLevels={1,1}> => array_size(in);

      group trigger_array[size] {
         GMOD.copy_on_change copy {
            trigger<NEportLevels={3,1}> => <-.<-.in[<-.index];
            input => (<-.index);
            output<NEportLevels={2,2}> => <-.<-.selected;
         };
         int index<NEportLevels=1> => index_of(trigger_array);
      };

      int selected<NEportLevels={1,2}>;
      int out<NEportLevels={1,2}> => in[selected];
   };
*/



   // Modifies the behaviour of the standard copy_on_change module.
   // Normally copy_on_change updates the output when the input OR trigger is updated.
   // This module changes this so that the module only updates the output when the trigger is updated.
   // Updating the input does not copy the value to the output.
   module copy_on_trigger {
      prim+read+IPort2 trigger; // IAC: the only way to get this module to trigger
      prim+Iparam+nonotify input;   // IAC: nonotify added to prevent input triggering copy
      prim+Oparam output;
      int+Iparam on_inst = 1;
      omethod+notify_val+notify_inst copy_on_change = "GMODcopy_on_change";
   };



   // Implements a link that can be broke and unbroken 
   // when triggers are recieved on two inputs.
   // Useful for disconnecting areas of a network when they
   // are not needed.
   macro BreakableLink {
      ilink in<NEportLevels={2,1}>;
      ilink connect<NEportLevels={2,1}>;
      ilink disconnect<NEportLevels={2,1}>;

      IAC_PROJ.Common.CommonMods.BreakLinkCore BreakLinkCore {
         in => <-.in;
         make => <-.connect;
         break => <-.disconnect;
      };

      olink out<NEportLevels={1,2}> => BreakLinkCore.out;
   };


   // Implements a link that is toggled between broken and unbroken 
   // states when a trigger is recieved on the input.
   // Useful for disconnecting areas of a network when they
   // are not needed.
   macro ToggleableLink {
      ilink in;
      ilink trigger;

      IAC_PROJ.Common.CommonMods.TrueFalseCore TrueFalseCore {
         trigger => <-.trigger;
         state => (!<-.state);
         state_out => <-.state;
      };
      int state = 0;

      IAC_PROJ.Common.CommonMods.BreakLinkCore BreakLinkCore {
         in => <-.in;
         make => <-.TrueFalseCore.state_true;
         break => <-.TrueFalseCore.state_false;
      };

      olink out => .BreakLinkCore.out;
   };



   // Implements a link that is first made and then unbroken 
   // when a trigger is recieved on the input.
   // Useful for briefly connecting areas of a network so that
   // data can be updated.
   macro TemporaryLink {
      ilink in<NEportLevels={2,1}>;
      ilink trigger<NEportLevels={2,1}>;

      GMOD.parse_v connect_v {
         v_commands = "out => .in; out=>;";
         trigger => <-.trigger;
         relative => <-;
         sync = 1;
         mode = 1;
      };

      olink out<NEportLevels={1,2}>;
   };



   // Takes two primitive inputs, sorts them and outputs them
   // as a higher number and a lower number.
   group SortTwoValues {
      prim value1<NEportLevels={2,1}>;
      prim value2<NEportLevels={2,1}>;
      prim array[2] => {.value1,.value2};
      prim low<NEportLevels={1,2}> => min_array(.array);
      prim high<NEportLevels={1,2}> => max_array(.array);
   };


   // Sets an individual element of an array,
   // based on an index number and an update value
   macro SetArrayElement {
      link index<NEportLevels={2,1}>;
      link value<NEportLevels={2,1}>;
      link data_array<NEportLevels={2,2}>;
 
      prim ref;
 
      GMOD.parse_v parse_v {
         v_commands => "ref => data_array[index]; ref = value;";
         trigger => <-.index;
         relative => (<-.<-);
      };
   };


   // Utility macro that kills off parent of macro 
   // when triggered
   macro kill_parent {
      int trigger <NEportLevels={2,1}> = 0;
      string name_of_parent => name_of(<-.<-);

      GMOD.parse_v parse_v {
         v_commands => ("-" + <-.name_of_parent) + ";";
         trigger => <-.trigger;
         on_inst = 0;
         relative => (<-.<-.<-.<-);
      };
   };


   // Utility macro that kills off parent of macro 
   // when triggered
   kill_parent kill_grandparent {
      name_of_parent => name_of(<-.<-.<-);

      parse_v {
         relative => (<-.<-.<-.<-.<-);
      };
   };


   // Utility macro that kills off parent of macro 
   // when triggered
   kill_parent kill_greatgrandparent {
      name_of_parent => name_of(<-.<-.<-.<-);

      parse_v {
         relative => (<-.<-.<-.<-.<-.<-);
      };
   };

};

