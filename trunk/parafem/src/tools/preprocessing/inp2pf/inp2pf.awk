# @(#) inp2pf.awk - Basic conversion of Abaqus Input Deck .inp to ParaFEM input files .d .bnd .lds .dat
# @(#) Usage:awk -f inp2pf.awk <filename.inp>
# Author: Louise M. Lever (louise.lever@manchester.ac.uk)
# Version: 1.1.0
# Date: 2012-04-26

# CHANGES:
# v1.1.0
#   LML: Renamed as inp2pf for consistency
#   LML: Added support for DC3D8 elements; requires further support for thermal problems
#   LML: Redirected all log messages to /dev/stderr
# v1.0.7
#   LML: Added support for C3D20R elements; includes field counting for element lines
#   LML: Added element id renumbering - requires command line argument "-renumber" passing to inp2pf script
#        do_elements_renumber() added
#        do_elements() or do_elements_renumber() called based on MODE set by command line argument
#   LML: Changed do_node output to use printf to force empty values to 0.0
# v1.0.6
#   LML: Added support for nset export to single .nset file
#   LML: Fixed missing Bound count from .dat export
# v1.0.5:
#   LML: Added multiple material support; increments auto-material ID on each repeated *ELEMENT in input deck
#   LML: Fixed "1" removed from each element line in  output; uses mat_id instead; defaults to 1
# v1.0.4:
#   LML: Removed repeated *ELEMENTS keyword in .d for multiple *Elements in input deck
# v1.0.3:
#   LML: Moved partition_mode parameter to new inserted 3rd line of .dat file.
# v1.0.2:
#   LML: Appended partition_mode to 4th line of .dat file; defaults to ParaFEM partitioning mode
# v1.0.1:
#   LML: Fixed building of node set arrays (was missing last index)

# TODO:
#   Use field counting for all element types
#   Add support for more element types
#   Add check to *Nodes processing to prevent empty *NODES fields in output

# FUNCTIONS:
#
# BEGIN
# END
# 
# start_nodes()
# do_nodes()
# start_elements()
# do_elements()
# start_nset()
# do_nset()
# do_nset_gen()
# start_bc()
# do_bc()
# output_bnd()
# start_step()
# output_load()
# output_fixed()
# start_material()
# start_elastic()
# do_elastic()
# start_cload()
# do_cload()
#
# (MAIN BODY)

# BEGIN
# 
# Initialize AWK. Set RS to process DOS and UNIX newlines.
# Set FS to use comma as well as whitespace
# Prepare the various filenames for output.
# Report some settings to stdout
# Set some default ParaFEM .dat values

BEGIN {
  RS = "\r\n|\n";
  FS = "[,]+[ \t]*";
  OFS = " ";

  # default parsing mode (looking for keyword)
  mode = "KEYWORD";

  # process input filename
  abq_file = ARGV[1];
  file_ext_pos = match(abq_file,/.[iI][nN][pP]$/) - 1;
  base_filename = substr(abq_file,0,file_ext_pos);

  # generate output filenames
  d_file = base_filename ".d";
  bnd_file = base_filename ".bnd";
  lds_file = base_filename ".lds";
  fix_file = base_filename ".fix";
  nset_file = base_filename ".nset";
  dat_file = base_filename ".dat";

  # report filenames
  print "Abaqus Input Deck:", abq_file > "/dev/stderr";
  print "Output model file:", d_file > "/dev/stderr";
  print "Output bnd file:", bnd_file > "/dev/stderr";
  print "Output lds file:", lds_file > "/dev/stderr";
  print "Output fix file:", fix_file > "/dev/stderr";
  print "Output nset file:", nset_file > "/dev/stderr";
  print "Output dat file:", dat_file > "/dev/stderr";

  # Open model file <name>.d and generate header
  print "*THREE_DIMENSIONAL" > d_file;

  # Auto-increment replacement element id; ++ on each element line; no reset on additional *ELEMENTS block
  elem_id = 1;
  # Auto-increment replacement node id; ++ on each node line; no reset on additional *NODES block
  node_id = 1;
  # Auto-increment default material id; ++ on additional *ELEMENTS block
  mat_id = 1;

  # old node to new id mapping LUT
  # node_id_lut

  # some hard-coded ParaFEM values
  tol = 1.0e-06;
  limit = 2000;

  # Partition Mode Flag: default ParaFEM = 1; (Post-METIS = 2)
  partition_mode = 1;

  # Initialized counts
  node_count = 0;
  bound_count = 0;
  loaded_count = 0;
  fixed_count = 0;

  # State flags
  started_elements_output = 0;
}

# END
# 
# Input file now completely processed. Report collected values for .dat file.
# Write .dat file

END {
    # report dat file values
    print "Type:", dat_elem_type > "/dev/stderr";
    print "Order:", dat_order > "/dev/stderr";
    print "partition_mode: " partition_mode > "/dev/stderr";
    print "element_count: " element_count, "node_count: " node_count, "bound_count: " bound_count, "nip: " nip, "nod: " nod, "loaded_count: " loaded_count, "fixed_freedoms: " fixed_count > "/dev/stderr";
    print "youngs_mod: " youngs_mod, "poisson_ratio: " poisson_ratio, "tol: " tol, "limit: " limit > "/dev/stderr"; 
    
    # output dat file based on fixed and calculated values
    print "Generating .dat file" > "/dev/stderr";
    print dat_elem_type > dat_file;
    print dat_order > dat_file;
    print partition_mode > dat_file;
    print element_count, node_count, bound_count, nip, nod, loaded_count, fixed_count > dat_file;
    printf "%e %e %e %d\n", youngs_mod, poisson_ratio, tol, limit > dat_file;
    
    # export NSET files
    print node_set_count, "NSETS defined:" > "/dev/stderr";
    if( node_set_count == 6 ) {
	print "Assuming nset problem type of 'kubc' as there are SIX nsets defined" > "/dev/stderr";
	print "'kubc'" > nset_file;
    } else {
	print "Assuming nset problem type of 'none'" > "/dev/stderr";
	print "'none'" > nset_file;
    }
    print node_set_count > nset_file;
    # export each nset in turn
    for( ns=0; ns<node_set_count; ns++ ) {
	tmp_nset_file = "tmp." ns "." base_filename ".nset";
	nsnn = 0; # reset node count for this set
	for( name in node_sets ) {
	    if( node_sets[name] == ns ) {
		print "Exporting NSET: ", ns, name > "/dev/stderr";
		ns_id = ns;
		ns_name = name;
	    }
	}
	# export node indices for this nset
	for( nd in node_list ) {
	    if( match(node_list[nd],":" ns ":" ) ) {
		print nd | "sort -n >> " tmp_nset_file;
		nsnn++;
	    }
	}
	close( "sort -n >> " tmp_nset_file );
	system("echo *NSET " nsnn " " ns_name " >> " nset_file);
	system("cat " tmp_nset_file " >> " nset_file);
	system("rm -f " tmp_nset_file);
    }
    print "Complete." > "/dev/stderr";
}

# Functions
#
# Called below from MAIN BODY

# FUNCTION: start_nodes()
#
# Triggered by *Nodes keyword
# Add *NODES keyword to output .d file
# Enter NODE mode to trigger do_nodes() until new keyword found

function start_nodes() {
  print "Processing Nodes" > "/dev/stderr";
  print "*NODES" > d_file;
  if( !renumber ) {
      mode = "NODE";
  } else {
      mode = "NODE-RENUMBER";
  }
}

function do_nodes() {
  gsub(/^ */,"");
  printf "%d %g %g %g\n", $1, $2, $3, $4 > d_file;
  node_count++;
}

function do_nodes_renumber() {
  gsub(/^ */,"");
  printf "%d %g %g %g\n", node_id, $2, $3, $4 > d_file;
  node_id_lut[$1] = node_id;
  node_count++;
  node_id++;
}

function start_elements() {
  print "Processing Elements" > "/dev/stderr";
    
  # Only want to output *ELEMENTS once in D FILE
  if( started_elements_output == 0 ) {
    print "*ELEMENTS" > d_file;
    # mat_id = 1; set in START
    started_elements_output = 1;
  } else {
    # else Just carry on adding element entries
    # but increment the auto-material ID
    mat_id++;	
  }
  split( $2,elem_type,"=" );
  print "Processing Element Type", elem_type[2] > "/dev/stderr";
  if( elem_type[2] == "C3D4" ) {
    dat_elem_type = "'tetrahedron'";
    nip = 1; nod = 4;
  } else if( elem_type[2] == "C3D10" ) {
    dat_elem_type = "'tetrahedron'";
    nip = 8; nod = 10;
  } else if( elem_type[2] == "C3D10R" ) {
    dat_elem_type = "'tetrahedron'";
    nip = 1; nod = 10;
  } else if( elem_type[2] == "C3D8") {
    dat_elem_type = "'hexahedron'";
    nip = 8; nod = 8;
  } else if( elem_type[2] == "C3D8R") {
    dat_elem_type = "'hexahedron'";
    nip = 1; nod = 8;
  } else if( elem_type[2] == "C3D8I") {
    dat_elem_type = "'hexahedron'";
    nip = 1; nod = 8;
  } else if( elem_type[2] == "C3D20") {
    dat_elem_type = "'hexahedron'";
    nip = 27; nod = 20;
  } else if( elem_type[2] == "C3D20R") {
    dat_elem_type = "'hexahedron'";
    nip = 8; nod = 20;
  } else if( elem_type[2] == "DC3D8") {
    dat_elem_type = "'hexahedron'";
    nip = 8; nod = 8;
  }
  dat_order = 2; # abaqus ordering

  if( !renumber ) {
      mode = "ELEMENT";
  } else {
      mode = "ELEMENT-RENUMBER";
  }
}

# ---------------------------------------------------------------------------------------------------------------------
# do_elements()
# do_elements_renumber()
# ---------------------------------------------------------------------------------------------------------------------
# For the given Abaqus element type, output 1) the given ID or new ID, 2) the ParaFEM element type code (e.g., 3 20 1),
# 3) the given node indices and the current material ID.
#
# Notes: The ParaFEM element type codes comprise "NDIM NNODES 1" for standard types, and "NDIM NNODES 2" for those
# ending with R. #LML: expand on this
#
# LML: TODO: Needs additional element type support
# ---------------------------------------------------------------------------------------------------------------------

function do_elements() {
  if( elem_type[2] == "C3D4" ) {
    gsub(/^ */,"");
    print " ", $1, "3 4 1", $2, $3, $4, $5, mat_id > d_file;
    element_count++;
  } else if( elem_type[2] == "C3D8I" ) {
    gsub(/^ */,"");
    print " ", $1, "3 8 1", $2, $3, $4, $5, $6, $7, $8, $9, mat_id > d_file;
    element_count++;
  } else if( elem_type[2] == "C3D8" ) {
    gsub(/^ */,"");
    print " ", $1, "3 8 1", $2, $3, $4, $5, $6, $7, $8, $9, mat_id > d_file;
    element_count++;
  } else if( elem_type[2] == "DC3D8" ) {
    gsub(/^ */,"");
    print " ", $1, "3 8 1", $2, $3, $4, $5, $6, $7, $8, $9, mat_id > d_file;
    element_count++;
  } else if( elem_type[2] == "C3D20R" ) {
    gsub(/^ */,"");
    node_remain = 20;
    field_remain = NF - 1; # subtract one for first line (the element id)
    field = 2
    printf "%d %s ", $1, "3 20 2" > d_file;
    while( node_remain ) {
	if( field_remain ) {
	    if( $field ) {
		printf "%d ", $field > d_file;
		node_remain--;
	    }
	    field_remain--;
	    field++;
	} else {
	    getline;
	    gsub(/^ */,"");
	    field_remain = NF;
	    field = 1;
	}
    }
    printf "%d\n", mat_id > d_file;
    element_count++;
  } else {
    print "Element Type", elem_type[2], "Not Supported" > "/dev/stderr";
  }
}

function do_elements_renumber() {
  if( elem_type[2] == "C3D4" ) {
    gsub(/^ */,"");
    print " ", elem_id, "3 4 1", node_id_lut[$2], node_id_lut[$3], node_id_lut[$4], node_id_lut[$5], mat_id > d_file;
    element_count++;
    elem_id++;
  } else if( elem_type[2] == "C3D8I" ) {
    gsub(/^ */,"");
    print " ", elem_id, "3 8 1", node_id_lut[$2], node_id_lut[$3], node_id_lut[$4], node_id_lut[$5], node_id_lut[$6], node_id_lut[$7], node_id_lut[$8], node_id_lut[$9], mat_id > d_file;
    element_count++;
    elem_id++;
  } else if( elem_type[2] == "C3D8" ) {
    gsub(/^ */,"");
    print " ", elem_id, "3 8 1", node_id_lut[$2], node_id_lut[$3], node_id_lut[$4], node_id_lut[$5], node_id_lut[$6], node_id_lut[$7], node_id_lut[$8], node_id_lut[$9], mat_id > d_file;
    element_count++;
    elem_id++;
  } else if( elem_type[2] == "C3D20R" ) {
    gsub(/^ */,"");
    node_remain = 20;
    field_remain = NF - 1; # subtract one for first line (the element id)
    field = 2
    printf "%d %s ", $1, "3 20 2" > d_file;
    while( node_remain ) {
	if( field_remain ) {
	    if( $field ) {
		printf "%d ", node_id_lut[$field] > d_file;
		node_remain--;
	    }
	    field_remain--;
	    field++;
	} else {
	    getline;
	    gsub(/^ */,"");
	    field_remain = NF;
	    field = 1;
	}
    }
    printf "%d\n", mat_id > d_file;
    element_count++;
  } else {
    print "Element Type", elem_type[2], "Not Supported" > "/dev/stderr";
  }
}

function start_nset() {
  split($2,nset_name,"=");
  node_set_id = node_set_count++;
  node_sets[nset_name[2]] = node_set_id; 
  print "Processing Nset", nset_name[2], node_sets[nset_name[2]], node_set_id > "/dev/stderr";
  mode = "NSET";
  for( of=3; of<=NF; of++ ) {
    if( $of == "generate" ) mode="NSET-GEN";
  }
  return;
}

function do_nset() {
  gsub(/^ */,"");
  for(i=1;i<=NF;i++) {
    if( $i ) {
      if($i in node_list) {
        node_list[$i] = node_list[$i] node_set_id ":"; 
      } else {
        node_list[$i] = ":" node_set_id ":";
      }
    }
  }
}

function do_nset_gen() {
  gsub(/^ */,"");
  print "Generating NSET " node_set_id "[" $1 "," $2 "];" $3 > "/dev/stderr";
  for( ind=$1; ind<=$2; ind+=$3 ) {
    if( ind ) {
      if(ind in node_list) {
        node_list[ind] = node_list[ind] node_set_id ":"; 
      } else {
        node_list[ind] = ":" node_set_id ":";
      }
    }
  }
}

function start_bc() {
  if( step_count == 0 ) {
    print "Processing BC for Step 'Initial'" > "/dev/stderr";
  } else {
    print "Processing BC for Step " step_count > "/dev/stderr";
  }
  mode = "BOUNDARY";
}

function do_bc() {
  bc_nset_id = node_sets[$1];
  if( step_count == 0 ) { # initial step BC; hold XYZ?
    # check for ENCASTRE and PINNED
    if(($2 == "ENCASTRE") || ($2 == "PINNED")) {
      mask = "XYZ";
    } else {
      # determine if we require X, Y and Z holding
      mask = "";
      if($2 == 1) mask = "X";
      if(($2 <= 2) && ($3 >= 2)) mask = mask "Y";
      if($3 == 3) mask = mask "Z";
    }
    print "Restraining " mask " in " $1 ":" bc_nset_id " node set" > "/dev/stderr";
    for( nd in node_list ) { # append mask to bound_list
      if( match(node_list[nd],":" bc_nset_id ":" ) ) {
        bound_list[nd] = bound_list[nd] mask;
      }
    }
  } else { # Fixed Displacement -steps
    if($2 == 1) {
      print "Fixing X in " $1 " node set" > "/dev/stderr";
      for( nd in node_list ) {
        if( match(node_list[nd],":" bc_nset_id ":") ) {
          fix_i_list[nd] = $4;
        }
      }
    }
    if(($2 <= 2) && ($3 >= 2)) {
      print "Loading Y in " $1 " node set" > "/dev/stderr";
      for( nd in node_list ) {
        if( match(node_list[nd],":" bc_nset_id ":") ) {
          fix_j_list[nd] = $4;
        }
      }
    }
    if($3 == 3) {
      print "Fixing Z in " $1 " node set" > "/dev/stderr";
      for( nd in node_list ) {
        if( match(node_list[nd],":" bc_nset_id ":") ) {
          fix_k_list[nd] = $4;
        }
      }
    }

    return;
  }
}

function output_bnd() {
  for( nd in bound_list ) {
    print " ", nd, match(bound_list[nd],"X")?0:1, match(bound_list[nd],"Y")?0:1,match(bound_list[nd],"Z")?0:1 | "sort -t'\t' -n > " bnd_file;
    ++bound_count;
  }
  return;
}

function start_step() {
  if( step_count == 0 ) { output_bnd(); } # Step 'Initial' complete
  split($2,step_name,"=");
  step[step_name[2]] = ++step_count;
}

# FUNCTION: output_load()
#
# Outputs each set of loads per time-step, preceeded by the time-stamp.
# First step after the initial step, starts the output and requires "cat >" to clear existing content.
# Subsequent steps append to the file.
# NOTE: The time-stamp is disabled for the first-step as p121 DOES NOT SUPPORT multiple time-steps
#       Further time-steps should not be encountered but will append WITH the time-stamp

function output_load() {
  if( step_count == 1 ) {
    # print step_count | "cat > " lds_file;
    printf "" | "cat > " lds_file;
    close( "cat > " lds_file );
  } else {
    print step_count | "cat >> " lds_file;
    close( "cat >> " lds_file );
  }
  for( nd in node_list ) {
    if( (nd in load_i_list) || (nd in load_j_list) || (nd in load_k_list) ) {
      printf " %d %e %e %e\n", nd, (nd in load_i_list)?load_i_list[nd]:0.0, (nd in load_j_list)?load_j_list[nd]:0.0, (nd in load_k_list)?load_k_list[nd]:0.0 | "sort -n >> " lds_file;
      ++loaded_count;
    }
  }
  close( "sort -n >> " lds_file );
}

function output_fixed() {
  if( step_count == 1 ) {
    # print step_count | "cat > " fix_file;
    printf "" | "cat > " fix_file;
    close( "cat > " fix_file );
  } else {
    print step_count | "cat >> " fix_file;
    close( "cat >> " fix_file );
  }
  for(nd in node_list ) {
    if (nd in fix_i_list) {
      printf "%d 1 %e\n", nd, fix_i_list[nd] | "sort -n >> " fix_file;
      ++fixed_count;
    }
    if(nd in fix_j_list) {
      printf "%d 2 %e\n", nd, fix_j_list[nd] | "sort -n >> " fix_file;
      ++fixed_count;
    }
    if(nd in fix_k_list) {
      printf "%d 3 %e\n", nd, fix_k_list[nd] | "sort -n >> " fix_file;
      ++fixed_count;
    }
  }
  close( "sort -n >> " fix_file );
}
function start_material() {
  split($2,mat_name,"=");
  material[mat_name[2]] = ++mat_count;
}

function start_elastic() {
  print "Processing ELASTIC properties for material " mat_name[2] > "/dev/stderr";
  mode = "ELASTIC";
}

function do_elastic() {
  gsub(/^ */,"");
  youngs_mod = $1;
  poisson_ratio = $2;

  # need to store for multiple materials
}

function start_cload() {
  if( step_count == 0 ) {
    print "Processing CLOAD for Step 'Initial'" > "/dev/stderr";
  } else {
    print "Processing CLOAD for Step " step_count > "/dev/stderr";
  }
  mode = "CLOAD";
}

function do_cload() {
  cload_nset_id = node_sets[$1];
  if( step_count == 0 ) { # initial step CLOAD; hold XYZ?
    # determine if we require X, Y and Z holding
    mask = "";
    if($2 == 1) mask = "X";
    if($2 == 2) mask = "Y";
    if($2 == 3) mask = "Z";
    print "Concetrated Force " mask " in " $1 " node set" > "/dev/stderr";
    for( nd in node_list ) { # append mask to bound_list
      if( match(node_list[nd],":" cload_nset_id ":") ) {
        bound_list[nd] = bound_list[nd] mask;
      }
    }
  } else { # load-steps; load IJK?
    print "STEP " step_count;
    if($2 == 1) {
      print "Concentrated Force X in " $1 " node set" > "/dev/stderr";
      for( nd in node_list ) {
        if( match(node_list[nd],":" cload_nset_id ":") ) {
          load_i_list[nd] = -$3;
        }
      }
    }
    if($2 == 2) {
      print "Concentrated Force Y in " $1 " node set" > "/dev/stderr";
      for( nd in node_list ) {
        if( match(node_list[nd],":" cload_nset_id ":") ) {
          load_j_list[nd] = -$3;
        }
      }
    }
    if($2 == 3) {
      print "Concentrated Force Z in " $1 " node set" > "/dev/stderr";
      for( nd in node_list ) {
        if( match(node_list[nd],":" cload_nset_id ":") ) {
          load_k_list[nd] = $3;
        }
      }
    }

    return;
  }
}

# --------------------------------------------------------------------------------
# MAIN RULES
# --------------------------------------------------------------------------------

# ignore comments
/^\*\*/ { next; }

# match these keywords and enter appropriate mode
/^*(Node Output|NODE OUTPUT)/ { next; }
/^*(Element Output|ELEMENT OUTPUT)/ { next; }
/^*(Node|NODE)/ { start_nodes(); next; }
/^*(Element|ELEMENT)/ { start_elements(); next; }
/^*(Nset|NSET)/ { start_nset(); next; }
/^*(Boundary|BOUNDARY)/ { start_bc(); next; }
/^*(Step|STEP)/ { start_step(); next; }
/^*(End Step|END STEP)/ { output_load(); output_fixed(); next; }
/^*(Material|MATERIAL)/ { start_material(); next; }
/^*(Elastic|ELASTIC)/ { start_elastic(); next; }
/^*(Cload|CLOAD)/ { start_cload(); next; }

# bounce unhandled keywords
mode == "KEYWORD" && /^*[A-Z]/ { printf "** Line Ignored (%d fields): %s\n", NF, $0 > "/dev/stderr"; }

# revert to KEYWORD mode if unhandled keyword found
mode != "KEYWORD" && /^*[A-Z]/ { mode = "KEYWORD"; }

# process non-keyword lines based on mode
mode == "NODE" { do_nodes(); }
mode == "NODE-RENUMBER" { do_nodes_renumber(); }
mode == "ELEMENT" { do_elements(); }
mode == "ELEMENT-RENUMBER" { do_elements_renumber(); }
mode == "NSET" { do_nset(); }
mode == "NSET-GEN" { do_nset_gen(); }
mode == "BOUNDARY" { do_bc(); }
mode == "ELASTIC" { do_elastic(); }
mode == "CLOAD" { do_cload(); }

# END OF MAIN
