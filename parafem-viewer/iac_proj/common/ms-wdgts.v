// INTERNATIONAL AVS CENTRE - WARRANTY DISCLAIMER
// Please read the file DISCLAIMER for conditions associated with this file.
// avs@iavsc.org, www.iavsc.org

flibrary MiscWidgets <compile_subs=0> {

   macro random_num {

      IAC_PROJ.Common.CommonMods.RandomNumParams RandomNumParams {
         num_vals = 10;
         min_val = 0;
         max_val = 1;
      };

      IAC_PROJ.Common.CommonMods.RandomNumCore RandomNumCore {
         RandomNumParams => <-.RandomNumParams;
      };


      IAC_PROJ.Common.UIWidgets.IAC_StandardUI IAC_StandardUI {
         panel {
            message = "Select Random Number control panel.";
         };

         mainTitleLabel {
            label = "Random Number Generator";
         };

         IAC_PROJ.Common.CommonMods.RandomNumParams &RandomNumParams<NEportLevels={2,1}> => <-.RandomNumParams;

         UIslider num_vals_slider {
            parent => <-.panel;
            value => <-.RandomNumParams.num_vals;
            title => "Number of Random Values";
            min = 1.;
            max = 50.;
            mode = "integer";
            decimalPoints = 0;
         };
         UIslider min_val_slider {
            parent => <-.panel;
            value => <-.RandomNumParams.min_val;
            title => "Minimum Random Value";
            min = -100.;
            max = 100.;
            decimalPoints = 2;
         };
         UIslider max_val_slider {
            parent => <-.panel;
            value<NEportLevels={2,2}> => <-.RandomNumParams.max_val;
            min = -100.;
            max = 100.;
            decimalPoints = 2;
            title => "Maximum Random Value";
         };
      };

      olink out_vals => .RandomNumCore.out_vals;
   };

};

