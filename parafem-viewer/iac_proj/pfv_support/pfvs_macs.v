flibrary PFVS_Macs {

  // GROUP: [Application Use]: Parameter Block

  group FEARXP_Param {
    // Model control

    enum input {
      choices = {"cells", "faces"};
    } = "faces";
    enum mesh {
      choices = {"model", "offset"};
    } = "model";
    enum data {
      choices = {"none", "cell", "scalar", "component"};
    } = "none";

    // Input/Model Data Component Control

    enum vec_comp { // node data
      choices = {"x", "y", "z"};
    } = "x";
    int cell_data_comp = 0;
    int+nres vector;

    // Mapped Data Controls
    
    int node_data_comp = 0;
    int glyph_data_comp = 0;
    int iso_data_comp = 0;

    // MinMax Control
    
    int lock_minmax = 0;
    double+nres min;
    double+nres max;

    // Time Step, Sub Step and Animation Controls

    int time = 0;
    int sub_time = 0;
    int time_step = 0;
    int time_sub_step = 0;
    int time_minmax[2];
    int sub_time_minmax[2];

    // Deformation Offset Control
    
    int offset_comp = 0;
    float offset_scale = 1.0;

    // Slice and Cut Control
    int slice = 0;
    int cut = 0;
    int cut_above = 0;
    enum slice_axis {
      choices = {"x", "y", "z", "Free"};
    } = "x";
    float slice_plane = 0;

    // Glyphing Control
    
    int glyph = 0;
    int glyph_geom = 0;
    int glyph_normalize = 0;
    float glyph_scale = 1.0;

    // Additional Style Control
    
    int no_lighting => (.slice && !.cut);

    // AABB Control
    int bounds = 0;
    int bounds_deformed = 0;
    int bounds_slicecut = 1;
  };

  // GROUP: Use as ARRAY: For per component user minmax storage
  
  group DataMinMax {
    double min;
    double max;
    int veclen;
    double min_vec[];
    double max_vec[];
  };

  // GROUP: Use as single master, multiple refs: Holds data and user minmax values for all components
  
  group Data_Info {
    Mesh+Node_Data &in;

    DataMinMax node_data[.in.nnode_data] {
      min+nres => <-.in.node_data[index_of(<-.node_data)].min;
      max+nres => <-.in.node_data[index_of(<-.node_data)].max;
      veclen+nres => <-.in_node_data[index_of(<-.node_data)].veclen;
      min_vec+nres => <-.in.node_data[index_of(<-.node_data)].min_vec;
      max_vec+nres => <-.in.node_data[index_of(<-.node_data)].max_vec;
    };
    DataMinMax cell_data[.in.cell_set[0].ncell_data] {
      min => <-.in.cell_set[0].cell_data[index_of(<-.cell_data)].min;
      max => <-.in.cell_set[0].cell_data[index_of(<-.cell_data)].max;
      // vector cell data not currently supported
    };

    DataMinMax user_node_data[.in.nnode_data] {
      min = 0.0;
      max = 1.0;
      veclen+nres => <-.in.node_data[index_of(<-.user_node_data)].veclen;
      min_vec => init_array(veclen,0.0,0.0);
      max_vec => init_array(veclen,1.0,1.0);
    };
    DataMinMax user_cell_data[.in.cell_set[0].ncell_data] {
      min = 0.0;
      max = 1.0;
    };

    int lock_node_data_minmax[.in.nnode_data] => init_array(.in.nnode_data,0,0);
    int lock_cell_data_minmax[.in.cell_set[0].ncell_data] => init_array(.in.cell_set[0].ncell_data,0,0);

    group lock_vec_node_data_minmax[.in.nnode_data] {
      int+nres veclen => <-.in.node_data[index_of(<-.lock_vec_node_data_minmax)].veclen;
      int lock[veclen] => init_array(veclen,0,0);
    };
    
    string node_data_labels[.in.nnode_data] => .in.node_data.labels;
    string cell_data_labels[.in.cell_set[0].ncell_data] => .in.cell_set[0].cell_data.labels;

    string node_data_units[.in.nnode_data] => .in.node_data.units;
    string cell_data_units[.in.cell_set[0].ncell_data] => .in.cell_set[0].cell_data.units;
  };
  
  // MACRO: [Application Use]: Map Data
  // Take output of viz macros and set what data is used for data mapping
  // Selects between "none", "element", "scalar", "component"

  macro Map_Data {
    ilink param;
    ilink in_info;
    ilink in;

    int show_node_data => (.param.data > 1);
    int show_cell_data => (.param.data == 1);
    int process_all_comps => (<-.param.glyph == 1);
    int node_data_comp => switch(.process_all_comps+1,0,.param.node_data_comp);

    string+nres labels => .in_info.labels;

    Node_Data no_node_data<NEportLevels={0,1}> {
      nnode_data = 0;
      nnodes+nres => <-.in.nnodes;
    };

    FLD_MAP.node_vector node_data {
      in_data+nres => <-.in.node_data[<-.node_data_comp].values;
      out {
	node_data {
	  labels+nres => <-.<-.<-.labels;
	  min+nres => <-.<-.<-.in_info.node_data[param.node_data_comp].min;
	  max+nres => <-.<-.<-.in_info.node_data[param.node_data_comp].max;
	};
      };    
    };
    
    group magnitude {
      Node_Data &in<NEportLevels={2,0}> => <-.node_data.out;
      Node_Data out<NEportLevels={0,2}> {
	nnodes+nres => in.nnodes;
	nnode_data = 1;
	node_data {
	  values => cache(magnitude(<-.<-.in.node_data[0].values));
	  veclen = 1;
	  labels => ("Magnitude" + switch((is_valid(<-.<-.<-.labels) + 1),"",(" of " + <-.<-.<-.labels)));
	  min+nres => <-.<-.<-.in_info.node_data[param.node_data_comp].min;
	  max+nres => <-.<-.<-.in_info.node_data[param.node_data_comp].max;
	};
      };
    };

    ExtractScalar.DVextr_scalar vec_comp {
      in => <-.node_data.out;
      component+nres => <-.param.vec_comp;
      out {
	node_data  {
	  min+nres => <-.<-.<-.in_info.node_data[param.node_data_comp].min_vec[param.vec_comp];
	  max+nres => <-.<-.<-.in_info.node_data[param.node_data_comp].max_vec[param.vec_comp];
	};
      };
    };

    DVM.DVswitch select_scalar {
      in => {<-.node_data.out, <-.magnitude.out};
      index+nres => <-.param.vector;
    };

    DVM.DVswitch select_data {
      // LML: Added 2nd no_node to get cell data; NEED colour/nocolour
      in => {<-.no_node_data, <-.no_node_data, <-.select_scalar.out, <-.vec_comp.out};
      index+nres => <-.param.data;
    };

    FLD_MAP.combine_mesh_data combine_mesh_data {
      in_mesh => <-.in;
      in_nd => <-.select_data.out;
    };

    olink out => .combine_mesh_data.out;
  };
  
  // Datamap import and export
  
  GDM.DataObject DataObjectExportDatamap {
    Datamap<NEportLevels={0,2}>;
  };
  
  GDM.DataObject DataObjectImportDatamap {
    DefaultLinear &Datamap<NEportLevels={2,0}>;
  };

  // Simple AABB object, used for improved slice/cut of AABB
  
  macro QuadBox {
    group+IPort2 &in_field {
      int nspace = 3;
      group coordinates {
	prim min_vec[3];
	prim max_vec[3];
      };
      group+opt xform;
    };
    
    Quad set1 {
      ncells = 6;
      node_connect_list = {0,1,2,3,4,5,6,7,0,1,5,4,3,2,6,7,0,3,7,4,1,5,6,2};
    };
    
    Mesh quadbox {
      int nnodes = 8;
      int nspace = 3;
      coordinates {
	float values[nvals][veclen] => {
	  {
	    <-.<-.in_field.coordinates.min_vec[0],
	    <-.<-.in_field.coordinates.min_vec[1],
	    <-.<-.in_field.coordinates.min_vec[2]
	  },
	  {
	    <-.<-.in_field.coordinates.max_vec[0],
	    <-.<-.in_field.coordinates.min_vec[1],
	    <-.<-.in_field.coordinates.min_vec[2]
	  },
	  {
	    <-.<-.in_field.coordinates.max_vec[0],
	    <-.<-.in_field.coordinates.max_vec[1],
	    <-.<-.in_field.coordinates.min_vec[2]
	  },
	  {
	    <-.<-.in_field.coordinates.min_vec[0],
	    <-.<-.in_field.coordinates.max_vec[1],
	    <-.<-.in_field.coordinates.min_vec[2]
	  },
	  {
	    <-.<-.in_field.coordinates.min_vec[0],
	    <-.<-.in_field.coordinates.min_vec[1],
	    <-.<-.in_field.coordinates.max_vec[2]
	  },
	  {
	    <-.<-.in_field.coordinates.max_vec[0],
	    <-.<-.in_field.coordinates.min_vec[1],
	    <-.<-.in_field.coordinates.max_vec[2]
	  },
	  {
	    <-.<-.in_field.coordinates.max_vec[0],
	    <-.<-.in_field.coordinates.max_vec[1],
	    <-.<-.in_field.coordinates.max_vec[2]
	  },
	  {
	    <-.<-.in_field.coordinates.min_vec[0],
	    <-.<-.in_field.coordinates.max_vec[1],
	    <-.<-.in_field.coordinates.max_vec[2]
	  }
	};
      };
      int ncell_sets = 1;
      cell_set[ncell_sets] => {set1};
    };
    
    DataObjectLite obj {
      in => <-.quadbox;
      Obj.name => name_of(<-.<-.<-);
      Modes {
	// Render the box as a wireframe
	mode = {0, GD_LINES, GD_NO_SURF, 0, 0};
      };
    };
    
    olink out_fld<export_all=2> => quadbox;
    olink out_obj => obj.Obj;
  };

  // Implementation of "missing" set_normals macro
  
  group DV_Param_set_normals<export_all=2> {
    int+Port2 component<animate=1>;
  };

  macro set_normals {
    ilink in_field<export_all=1>;

    DV_Param_set_normals SetNormalsParam<export_all=2> {
      component = 0;
    };

    DVextract_comp DVextract_comp {
      in => <-.in_field;
      component => <-.SetNormalsParam.component;
    };
    
    group extract_normals_data_array {
      Node_Data+nres+IPort2 &in => <-.DVextract_comp.out;
      prim+nres+OPort2 data[] => in.node_data[0].values;
    };
    
    FLD_MAP.node_normals normals_data {
      in_data => <-.extract_normals_data_array.data;
      out.node_data.labels = "normals";
    };
    
    DVcomb_mesh_and_data DVcomb_mesh_and_data {
      in_mesh => <-.in_field;
      in_nd => <-.normals_data.out;
    };

    GMOD.instancer instancer {
      Value => UIpanel.visible;
      active = 2; // don't de-instance when visible = 0
      Group => SetNormalsUI;
    };

    UImod_panel UIpanel {
      parent<NEportLevels={3,0}>;
      title => name_of(<-.<-);
      message = "Select set_alpha panel.";
    };

    macro SetNormalsUI<instanced=0> {
      ilink in_fld => in_field;
      DV_Param_set_normals+IPort2 &param => SetNormalsParam;
      
      DVnode_data_labels DVnode_data_labels {
	in => in_fld;
	int+nres  ncomp => in.nnode_data;
      };
      
      ilink UIpanel  => <-.UIpanel;
      
      UIradioBoxLabel normals_radioBoxLabel {
	parent => <-.UIpanel;
	labels+IPort2 => DVnode_data_labels.labels;
	&selectedItem+IPort2 => param.component;
	visible => <-.UIpanel.visible;
	title = "normals component";
	y =  0;
	width	=> <-.UIpanel.width;
      };
    };

    DataObjectNoTexture obj {
      in => DVcomb_mesh_and_data.out;
      Obj {
	name => name_of(<-.<-.<-);
      };
    };

    olink out_fld<export_all=2> => DVcomb_mesh_and_data.out;
    olink out_obj => obj.obj;
  };
};
