// INTERNATIONAL AVS CENTRE - WARRANTY DISCLAIMER
// Please read the file DISCLAIMER for conditions associated with this file.
// avs@iavsc.org, www.iavsc.org

flibrary IOWidgets <compile_subs=0> {

  // IndexedFileParams
  //   Parameter Group
  //     Holds the parameters necessary to generate an indexed filename

  group IndexedFileParams <NEportLevels={0,1}> {
      string dir<NEportLevels={2,2},export=1>;
      string filename_stub<NEportLevels={2,2},export=1>;
      string extension<NEportLevels={2,2},export=1>;
      int index_width<NEportLevels={2,2},export=1>;
  };


   // CreateIndexedFilenameUI
   //   UI Macro
   //   UImod_panel containing 3 SmartTextTypein macros, 1 SmartFieldTypein macro
   //   and associated strings. Used to setup and alter the IndexedFileParams group.

  IAC_PROJ.Common.UIWidgets.IAC_StandardUI CreateIndexedFilenameUI<NEx=55.,NEy=77.> {

      IAC_PROJ.Common.IOWidgets.IndexedFileParams &IndexedFileParams<NEx=198.,NEy=11.,NEportLevels={2,1}>;

      panel<NEx=33.,NEy=99.> {
         message = "Select Indexed Filename control panel.";
         title => name_of(<-.<-.<-) + ": Filename Strings";
      };
      mainTitleLabel<NEx=187.,NEy=220.> {
         label = "Indexed Filename Control Panel";
      };
      SmartTypeinParams<NEx=231.,NEy=374.>;

      IAC_PROJ.Common.UIWidgets.SmartTextTypein DirectorySmartTextTypein<NEx=462.,NEy=88.> {
         UIparent => <-.panel;
         slabel => <-.dirString;
         stext => <-.IndexedFileParams.dir;
      };
      IAC_PROJ.Common.UIWidgets.SmartTextTypein StubSmartTextTypein<NEx=462.,NEy=132.> {
         UIparent => <-.panel;
         slabel => <-.stubString;
         stext => <-.IndexedFileParams.filename_stub;
      };
      IAC_PROJ.Common.UIWidgets.SmartTextTypein ExtSmartTextTypein<NEx=462.,NEy=176.> {
         UIparent => <-.panel;
         slabel => <-.extString;
         stext => <-.IndexedFileParams.extension;
      };

      IAC_PROJ.Common.UIWidgets.SmartFieldTypein SmartFieldTypein<NEx=462.,NEy=220.> {
         UIparent => <-.panel;
         flabel => <-.widthString;
         fval => <-.IndexedFileParams.index_width;
         fmin = 0;
         fmax = 20;
      };

      string dirString<NEportLevels=1,NEx=671.,NEy=88.> = "Directory";
      string stubString<NEportLevels=1,NEx=671.,NEy=132.> = "Filename Stub";
      string extString<NEportLevels=1,NEx=671.,NEy=176.> = "File Extension";
      string widthString<NEportLevels=1,NEx=671.,NEy=220.> = "Width of Index";
  };


  // CreateIndexedFilename
  //   Functional Group
  //     Generates the required filename by using a string formatting and string conatenation.

  group CreateIndexedFilename {
      IndexedFileParams &FileParams <NEportLevels={2,0}>;

      int index<NEportLevels={2,0},export=1>;

      string format => "%0" + FileParams.index_width + "d";

      string filename<NEportLevels={1,2}> => FileParams.dir + FileParams.filename_stub + str_format(.format,.index) + "." + FileParams.extension;
  };


  // createIndexedFilename
  //   User Macro
  //   Combines the dynamic filename generation module with a User Interface
  //   and a set of parameters

  macro createIndexedFilename<NEx=187.,NEy=88.> {
      ilink index<NEx=187.,NEy=33.>;

      IAC_PROJ.Common.IOWidgets.IndexedFileParams IndexedFileParams<NEx=396.,NEy=77.> {
         index_width = 0;
      };
      IAC_PROJ.Common.IOWidgets.CreateIndexedFilename CreateIndexedFilename<NEx=187.,NEy=165.> {
         FileParams => <-.IndexedFileParams;
         index => <-.index;
      };
      IAC_PROJ.Common.IOWidgets.CreateIndexedFilenameUI CreateIndexedFilenameUI<NEx=374.,NEy=165.> {
         IndexedFileParams => <-.IndexedFileParams;
      };

      olink filename<NEx=99.,NEy=286.> => .CreateIndexedFilename.filename;
  };


  // AnimFilenameUI
  //   UI Macro
  //   UImod_panel containing 3 SmartFieldTypein macros, 1 toggle control,
  //   2 buttons and associated strings.
  //   Used to setup and alter the AnimLoopParams group.

   IAC_PROJ.Common.UIWidgets.IAC_StandardUI AnimFilenameUI {
      panel<NEx=22.,NEy=132.> {
         message = "Select AnimFilename control panel.";
         title => name_of(<-.<-.<-) + ": Animation";
      };

      mainTitleLabel {
         label = "Animate Filename Control Panel";
      };
      IAC_PROJ.Common.UIWidgets.titleLabel titleLabel<NEx=132.,NEy=275.> {
         parent => <-.panel;
         label = "Loop Parameters";
      };

      SmartTypeinParams<NEx=165.,NEy=55.>;

      IAC_PROJ.Common.UIWidgets.VideoLoopParams &VideoLoopParams<NEx=539.,NEy=385.,NEportLevels={2,1}>;

      IAC_PROJ.Common.UIWidgets.SmartFieldTypein StartSmartFieldTypein<NEx=473.,NEy=132.> {
         UIparent => <-.panel;
         flabel => <-.startString;
         fval => <-.VideoLoopParams.start_val;
      };
      IAC_PROJ.Common.UIWidgets.SmartFieldTypein EndSmartFieldTypein<NEx=473.,NEy=176.> {
         UIparent => <-.panel;
         flabel => <-.endString;
         fval => <-.VideoLoopParams.end_val;
      };
      IAC_PROJ.Common.UIWidgets.SmartFieldTypein IncrSmartFieldTypein<NEx=473.,NEy=220.> {
         UIparent => <-.panel;
         flabel => <-.incrString;
         fval => <-.VideoLoopParams.incr;
      };


      string startString<NEportLevels=1,NEx=649.,NEy=132.> = "Start Value";
      string endString<NEportLevels={1,1},NEx=649.,NEy=176.> = "End Value";
      string incrString<NEportLevels=1,NEx=649.,NEy=220.> = "Increment Size";


      UItoggle CycleUItoggle<NEx=308.,NEy=374.> {
         parent => <-.panel;
         label => "Cycle Continously";
         width => parent.width;
         set<NEportLevels={2,2}> => <-.VideoLoopParams.cycle;
      };

      UIslider CountUIslider<NEx=616.,NEy=341.> {
         parent => <-.panel;
         width => parent.width;

         value<NEportLevels={2,2}> => <-.VideoLoopParams.count;
         title => "Frame No:";
         min<NEportLevels={2,0}> => <-.VideoLoopParams.start_val;
         max<NEportLevels={2,0}> => <-.VideoLoopParams.end_val;
         decimalPoints = 0;
         increment = 10.;
      };

      IAC_PROJ.Common.UIWidgets.VideoControl VideoControl<NEx=836.,NEy=341.> {
         UIparent => <-.panel;
         VideoLoopParams => <-.VideoLoopParams;
      };

   };


  // AnimFilename
  //   Functional Macro
  //   Combines the altered loop module, VideoLoop, with the dynamic filename
  //   generation module, CreateIndexedFilename.

  macro AnimFilename {

    IAC_PROJ.Common.IOWidgets.IndexedFileParams &IndexedFileParams<NEx=22.,NEy=77.,NEportLevels={2,1}>;
    IAC_PROJ.Common.UIWidgets.VideoLoopParams &VideoLoopParams<NEx=176.,NEy=33.,NEportLevels={2,1}>;

    IAC_PROJ.Common.UIWidgets.VideoLoop VideoLoop<NEx=209.,NEy=88.> {
        VideoLoopParams => <-.VideoLoopParams;
    };

    IAC_PROJ.Common.IOWidgets.CreateIndexedFilename CreateIndexedFilename<NEx=99.,NEy=165.> {
        FileParams => <-.IndexedFileParams;
        index => <-.VideoLoop.count;
    };

    olink filename<NEx=99.,NEy=231.> => .CreateIndexedFilename.filename;
  };


  // animFilename
  //   User Macro
  //   uses both parameter groups, the AnimFilenameGen module, and
  //   the modified GMOD.loop module

  macro animFilename {
    IAC_PROJ.Common.UIWidgets.VideoLoopParams VideoLoopParams <NEx=22.,NEy=22.> {
      start_val = 0;
      end_val = 10;
      incr = 1;
    };

    IAC_PROJ.Common.IOWidgets.IndexedFileParams IndexedFileParams<NEx=22.,NEy=66.> {
      dir = "build3d/";
      filename_stub = "cube";
      extension = "tif";
      index_width = 0;
    };

    IAC_PROJ.Common.IOWidgets.AnimFilenameUI AnimFilenameUI <NEx=286.,NEy=99.> {
      VideoLoopParams => <-.VideoLoopParams;
      UImod_panel {
        title => name_of(<-.<-.<-);
      };
    };

    IAC_PROJ.Common.IOWidgets.CreateIndexedFilenameUI CreateIndexedFilenameUI<NEx=22.,NEy=143.> {
      IndexedFileParams => <-.IndexedFileParams;
    };

    IAC_PROJ.Common.IOWidgets.AnimFilename AnimFilename<NEx=198.,NEy=143.> {
      VideoLoopParams => <-.VideoLoopParams;
      IndexedFileParams => <-.IndexedFileParams;
    };

    link out <NEportLevels={1,2},NEx=198.,NEy=231.> => AnimFilename.filename;
  };

};

