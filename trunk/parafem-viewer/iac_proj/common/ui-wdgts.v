// INTERNATIONAL AVS CENTRE - WARRANTY DISCLAIMER
// Please read the file DISCLAIMER for conditions associated with this file.
// avs@iavsc.org, www.iavsc.org

flibrary UIWidgets <compile_subs=0> {

   group+GMOD.hconnect SmartTypeinParams <NEportLevels={0,1}> {
      int total_width <NEportLevels={2,2},export=1> = 250;
      int label_width <NEportLevels={2,2},export=1> = 70;
      int offset      <NEportLevels={2,2},export=1> = 10;
      accept = "smartTypein";  // accept all offer = "UIcomp" children
      accept_name = ".";       // connect to parent_shell itself
   };

   macro SmartTextTypein {
      ilink  UIparent;
      ilink  slabel;
      ilink  stext;

      SmartTypeinParams &params<NEportLevels={2,1}>;

      UIpanel panel {
         parent => UIparent;
         width => <-.params.total_width;
         height => <-.text.height;
      };

      UIlabel text_label {
         parent => panel;
         label => slabel;
         alignment = 2;
         x = 0;
         y => ((<-.text.height - .height) / 2);
         width => <-.params.label_width;
      };

      UItext text {
         parent => panel;
         &text => stext;
         x => <-.params.label_width + <-.params.offset;
         y = 0;
         width => <-.panel.width - .x;
         height => ((UIdata.UIfonts[0].lineHeight + 0) * 1.5);
      };

      GMOD.hconnect UIhconnect {
         direction = 0;              // connection goes from child to parent
         offer = "UIcomp";           // connect to "UIComp" class
         offer_name = "<-.UIparent"; // connect UIparent subobject of SmartTypein
      };

      GMOD.hconnect Smart_hconnect {
         direction = 0;              // connection goes from child to parent
         offer = "smartTypein";      // connect to "smartTypein" class
         offer_name = "<-.params";   // connect UIparent subobject of SmartTypein
         skip_levels = 1;
      };
   };


   macro SmartFieldTypein {
      ilink  UIparent;
      ilink  flabel;
      ilink  fval;
      ilink  fmin;
      ilink  fmax;

      SmartTypeinParams &params<NEportLevels={2,1}>;

      UIpanel panel {
         parent => UIparent;
         width => <-.params.total_width;
         height => <-.field.height;
      };

      UIlabel field_label {
         parent => panel;
         label => flabel;
         alignment = 2;
         x = 0;
         y => ((<-.field.height - .height) / 2);
         width => <-.params.label_width;
      };

      UIfield field {
         parent => panel;
         min => fmin;
         max => fmax;
         value => fval;
         x => <-.params.label_width + <-.params.offset;
         y = 0;
         width => <-.panel.width - .x;
      };

      GMOD.hconnect UIhconnect {
         direction = 0;              // connection goes from child to parent
         offer = "UIcomp";           // connect to "UIComp" class
         offer_name = "<-.UIparent"; // connect UIparent subobject of SmartTypein
      };

      GMOD.hconnect Smart_hconnect {
         direction = 0;              // connection goes from child to parent
         offer = "smartTypein";      // connect to "smartTypein" class
         offer_name = "<-.params";   // connect UIparent subobject of SmartTypein
         skip_levels = 1;
      };
   };


   group VideoLoopParams <NEportLevels={0,1}> {
      int reset         <NEportLevels={2,2},export=1> = 0;
      int reset_back    <NEportLevels={2,2},export=1> = 0;
      int run           <NEportLevels={2,2},export=1> = 0;
      int run_back      <NEportLevels={2,2},export=1> = 0;
      int step          <NEportLevels={2,2},export=1> = 0;
      int step_back     <NEportLevels={2,2},export=1> = 0;
      int cycle         <NEportLevels={2,2},export=1> = 0;
      double start_val  <NEportLevels={2,2},export=1>;
      double end_val    <NEportLevels={2,2},export=1>;
      double incr       <NEportLevels={2,2},export=1> = 1;
      double count      <NEportLevels={2,2},export=1>;
   };

   GMOD.loop VideoLoop {
      VideoLoopParams &VideoLoopParams<NEportLevels={2,1}>;
      reset       <NEportLevels={0,2}> => VideoLoopParams.reset;
      reset_back  <NEportLevels={0,2}> => VideoLoopParams.reset_back;
      run         <NEportLevels={0,0}> => VideoLoopParams.run;
      run_back    <NEportLevels={0,0}> => VideoLoopParams.run_back;
      step        <NEportLevels={0,0}> => VideoLoopParams.step;
      step_back   <NEportLevels={0,0}> => VideoLoopParams.step_back;
      cycle       <NEportLevels={0,0}> => VideoLoopParams.cycle;
      done;
      start_val   <NEportLevels={0,0}> => VideoLoopParams.start_val;
      end_val     <NEportLevels={0,0}> => VideoLoopParams.end_val;
      incr        <NEportLevels={0,0}> => VideoLoopParams.incr;
      count       <NEportLevels={0,2}> => VideoLoopParams.count;
   };


   macro VideoControl {

      ilink  UIparent<NEx=99.,NEy=66.>;
      int  button_size<NEx=99.,NEy=209.> = 28;

      VideoLoopParams &VideoLoopParams<NEportLevels={2,1},NEx=407.,NEy=264.>;

      UIpanel panel<NEx=605.,NEy=44.> {
         parent => UIparent;
         width => parent.width;
         height => button_size + 8;
      };

      UIbutton ResetUIbutton<NEx=187.,NEy=121.> {
         parent => <-.panel;
         x = 0;
         y = 4;
         height => <-.button_size;
         width => <-.button_size;
         label => "";
         labelPixmap {
            filename<NEdisplayMode="open"> = "iac_proj/common/lend.xbm";
         };
         do<NEportLevels={2,2}> => <-.VideoLoopParams.reset;
      };

      UIbutton PlayBackUIbutton<NEx=341.,NEy=121.> {
         parent => <-.panel;
         x => ((<-.ResetUIbutton.x + <-.button_size) + 4);
         y = 4;
         height => <-.button_size;
         width => <-.button_size;
         label => "";
         labelPixmap {
            filename<NEdisplayMode="open"> = "iac_proj/common/lplay.xbm";
         };
         do<NEportLevels={2,2}> => <-.VideoLoopParams.run_back;
      };

      UIbutton StepBackUIbutton<NEx=495.,NEy=121.> {
         parent => <-.panel;
         x => ((<-.PlayBackUIbutton.x + <-.button_size) + 4);
         y = 4;
         height => <-.button_size;
         width => <-.button_size;
         label => "";
         labelPixmap {
            filename<NEdisplayMode="open"> = "iac_proj/common/larrow.xbm";
         };
         do<NEportLevels={2,2}> => <-.VideoLoopParams.step_back;
      };

      UIbutton StopUIbutton<NEx=649.,NEy=121.> {
         parent => <-.panel;
         x => ((<-.StepBackUIbutton.x + <-.button_size) + 4);
         y = 4;
         height => <-.button_size;
         width => <-.button_size;
         label => "";
         labelPixmap {
            filename<NEdisplayMode="open"> = "iac_proj/common/stop.xbm";
         };
         do<NEportLevels={2,2}> => (!<-.VideoLoopParams.run);
      };

      UIbutton StepUIbutton<NEx=341.,NEy=165.> {
         parent => <-.panel;
         x => ((<-.StopUIbutton.x + <-.button_size) + 4);
         y = 4;
         height => <-.button_size;
         width => <-.button_size;
         label => "";
         labelPixmap {
            filename<NEdisplayMode="open"> = "iac_proj/common/rarrow.xbm";
         };
         do<NEportLevels={2,2}> => <-.VideoLoopParams.step;
      };

      UIbutton PlayUIbutton<NEx=495.,NEy=165.> {
         parent => <-.panel;
         x => ((<-.StepUIbutton.x + <-.button_size) + 4);
         y = 4;
         height => <-.button_size;
         width => <-.button_size;
         label => "";
         labelPixmap {
            filename<NEdisplayMode="open"> = "iac_proj/common/rplay.xbm";
         };
         do<NEportLevels={2,2}> => <-.VideoLoopParams.run;
      };

      UIbutton ResetBackUIbutton<NEx=649.,NEy=165.> {
         parent => <-.panel;
         x => ((<-.PlayUIbutton.x + <-.button_size) + 4);
         y = 4;
         height => <-.button_size;
         width => <-.button_size;
         label => "";
         labelPixmap {
            filename<NEdisplayMode="open"> = "iac_proj/common/rend.xbm";
         };
         do<NEportLevels={2,2}> => <-.VideoLoopParams.reset_back;
      };

      GMOD.hconnect UIhconnect {
         direction = 0;              // connection goes from child to parent
         offer = "UIcomp";           // connect to "UIComp" class
         offer_name = "<-.UIparent"; // connect UIparent subobject of SmartTypein
      };
   };


   UIlabel titleLabel {
      width => parent.width;
      alignment = "left";
      label => "##Your title here##";

      GMOD.hconnect UIhconnect {
         direction = 0;              // connection goes from child to parent
         offer = "UIcomp";           // connect to "UIComp" class
         offer_name = "<-.parent";   // connect parent subobject of SmartTypein
      };
   };

   titleLabel mainTitleLabel {
      color {
         backgroundColor = "blue";
         foregroundColor = "white";
      };
   };


   macro IAC_StandardUI {
      UImod_panel+GMOD.hconnect panel {
         message = "Select ##Your Mod name## control panel.";
         title => name_of(<-.<-.<-);
         parent<NEportLevels={4,0}>;
         accept = "UIcomp";  // accept all offer = "UIcomp" children
         accept_name = ".";  // connect to parent_shell itself
      };

      mainTitleLabel mainTitleLabel {
         parent => <-.panel;
      };

      SmartTypeinParams SmartTypeinParams;
   };

};

