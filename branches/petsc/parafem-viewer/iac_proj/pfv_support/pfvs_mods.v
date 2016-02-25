flibrary PFVS_Mods < build_dir="iac_proj/pfv_support",
                     out_hdr_file="gen.hxx",
                     out_src_file="gen.cxx",
                     cxx_name="" > {

  group Comp_Mapping {
    int selected;
    int mapped;
  };
  
  module component_mapper < src_file="comp_map.cxx" > {
    Comp_Mapping &monitor_comps[];
    int reqd_comps[];
    cxxmethod+req update (
      monitor_comps+read+notify+req,
      reqd_comps+write
    );
  };

  // Additional: Modified modules for cut plane; cuts AABB too
  
  MODS.cut_plane cut_plane_bounds {
    ilink in_bounds<export_all=1>;

    int has_cell_data_bounds => (DVcell_data_labels_bounds.ncomp > 0);

    // Additional Cut for bounds field using same plane
    Cut Cut_Bounds {
      param => CutParam;
      in_fld => in_bounds;
      in_pln => <-.plane;
      DVcut {
	cell_data => switch(has_cell_data_bounds, param.cell_data);
      };
    };

    MODS.external_edges external_edges {
      in_field => <-.Cut_Bounds.out_fld;
    };
    
    // duplicate for bounds input
    DVnode_data_labels DVnode_data_labels_bounds {
      in => in_bounds;
      int+nres  ncomp => in.nnode_data;
    };
    
    DVcell_data_labels DVcell_data_labels_bounds {
      in => in_bounds;
    };

    // output original uncut bounds as dotted line
    DataObject uncut_bounds_obj {
      in => in_bounds;
      Modes {
	// Render the box as a wireframe
	mode = {0, GD_LINES, GD_NO_SURF, 0, 0};
      };
      Props {
	line_style = "Dotted";
	inherit = 0;
      };
      Obj {
	name => name_of(<-.<-.<-);
      };
    };

    // output cut bounds as soldid
    DataObject cut_bounds_obj {
      in => <-.external_edges.out_fld;
      Modes {
	mode = {0,2,1,0,0};
      };
      Obj {
	name => name_of(<-.<-.<-);
	Modes.normals = 1;
      };
    };

    // make plane_obj invisible by default now
    plane_obj {
      Obj {
	visible = 0;
      };
    };
    
    olink out_bounds_fld<export_all=2> => .external_edges.out_fld;
    olink out_bounds_obj_uncut => uncut_bounds_obj.obj;
    olink out_bounds_obj_cut => cut_bounds_obj.obj;
  };

  // Additional: Modified module for double cut using thickness; cuts AABB too

  MODS.cut_plane mirror_cut_plane_bounds {
    ilink in_bounds<export_all=1>;

    // Modify params to cut below
    CutParam {
      above = 0;
    };

    // Additonal parameter block mirroring original parameters
    DV_Param_cut+Port CutParamMirror<export_all=2> {
      component = {0};
      dist => (0.0 - <-.CutParam.dist);
      above = 1;
      cell_data = {0};
    };

    int has_cell_data_bounds => (DVcell_data_labels_bounds.ncomp > 0);

    // Additional Cut for bounds field using same plane
    Cut Cut_Bounds {
      param => CutParam;
      in_fld => in_bounds;
      in_pln => <-.plane;
      DVcut {
	cell_data => switch(has_cell_data_bounds, param.cell_data);
      };
    };

    // Additional SECOND Cut for field using MIRROR plane
    Cut Cut_Mirror {
      param => CutParamMirror;
      in_fld => Cut.out_fld;
      in_pln => <-.plane;
      DVcut {
	cell_data => switch(has_cell_data, param.cell_data);
      };
    };

    // Additional SECOND Cut for bounds field using MIRROR plane
    Cut Cut_Bounds_Mirror {
      param => CutParamMirror;
      in_fld => Cut_Bounds.out_fld;
      in_pln => <-.plane;
      DVcut {
	cell_data => switch(has_cell_data_bounds, param.cell_data);
      };
    };

    MODS.external_edges external_edges {
      in_field => <-.Cut_Bounds_Mirror.out_fld;
    };
    
    // duplicate for bounds input
    DVnode_data_labels DVnode_data_labels_bounds {
      in => in_bounds;
      int+nres  ncomp => in.nnode_data;
    };
    
    DVcell_data_labels DVcell_data_labels_bounds {
      in => in_bounds;
    };

    // output original uncut bounds as dotted line
    DataObject uncut_bounds_obj {
      in => in_bounds;
      Modes {
	// Render the box as a wireframe
	mode = {0, GD_LINES, GD_NO_SURF, 0, 0};
      };
      Props {
	line_style = "Dotted";
	inherit = 0;
      };
      Obj {
	name => str_format( "%s-%s", name_of(<-.<-.<-), "UncutBounds" );
      };
    };

    // output cut bounds as soldid
    DataObject cut_bounds_obj {
      in => <-.external_edges.out_fld;
      Modes {
	mode = {0,2,1,0,0};
      };
      Obj {
	name => str_format( "%s-%s", name_of(<-.<-.<-), "CutBounds" );
	Modes.normals = 1;
      };
    };

    // make plane_obj invisible by default now
    plane_obj {
      Obj {
	visible = 0;
	name => str_format( "%s-%s", name_of(<-.<-.<-), "CutPlaneManipulator" );
      };
    };

    // make out obj show SECOND cut object
    cut_obj {
      in => Cut_Mirror.out_fld;
      Obj {
	name => str_format( "%s-%s", name_of(<-.<-.<-), "CutObject" );
      };
    };

    out_fld => Cut_Mirror.out_fld;
    olink out_bounds_fld<export_all=2> => .external_edges.out_fld;
    olink out_bounds_obj_uncut => uncut_bounds_obj.obj;
    olink out_bounds_obj_cut => cut_bounds_obj.obj;
  };

//  module DmapGetAlpha < src_file="dmap_get_alpha.cxx",
//                        libdeps="DMAP" > {
//
//    DatamapTempl &dmap<NEportLevels={2,0}>;
//    float &values<NEportLevels={2,0}>[];
//    float alphas<NEportLevels={0,2}>[];
//			  
//    cxxmethod+notify_inst+req update {
//      dmap+read+notify+req,
//      values+read+notify+req,
//      alpha+write+nosave
//    };
//  };
};
