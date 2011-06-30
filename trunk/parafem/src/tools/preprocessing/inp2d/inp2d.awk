# @(#) inp2d.awk - Basic conversion of Abaqus Input Deck .inp to ParaFEM input files .d .bnd .lds .dat
# @(#) Usage:awk -f inp2d.awk <filename.inp>
# Author: Louise M. Lever (louise.lever@manchester.ac.uk)
# Version: 1.0.4

# CHANGES:
# v1.0.4:
#   LML: Removed repeated *ELEMENTS keyword in .d for multiple *Elements in input deck
# v1.0.3:
#   LML: Moved partition_mode parameter to new inserted 3rd line of .dat file.
# v1.0.2:
#   LML: Appended partition_mode to 4th line of .dat file; defaults to ParaFEM partitioning mode
# v1.0.1:
#   LML: Fixed building of node set arrays (was missing last index)

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
  dat_file = base_filename ".dat";

  # report filenames
  print "Abaqus Input Deck:", abq_file;
  print "Output model file:", d_file;
  print "Output bnd file:", bnd_file;
  print "Output lds file:", lds_file;
  print "Output fix file:", fix_file;
  print "Output dat file:", dat_file;

  # Open model file <name>.d and generate header
  print "*THREE_DIMENSIONAL" > d_file;

  # some hard-coded ParaFEM values
  tol = 1.0e-06;
  limit = 2000;

  # Partition Mode Flag: default ParaFEM = 1; (Post-METIS = 2)
  partition_mode = 1;

  # Initialized counts
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
  print "Type:", dat_elem_type;
  print "Order:", dat_order;
  print "partition_mode: " partition_mode;
  print "element_count: " element_count, "node_count: " node_count, "bound_count: " bound_count, "nip: " nip, "nod: " nod, "loaded_count: " loaded_count, "fixed_freedoms: " fixed_count;
  print "youngs_mod: " youngs_mod, "poisson_ratio: " poisson_ratio, "tol: " tol, "limit: " limit; 

  # output dat file based on fixed and calculated values
  print "Generating .dat file";
  print dat_elem_type > dat_file;
  print dat_order > dat_file;
  print partition_mode > dat_file;
  print element_count, node_count, bound_count, nip, nod, loaded_count, fixed_count > dat_file;
  printf "%e %e %e %d\n", youngs_mod, poisson_ratio, tol, limit > dat_file;
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
  print "Processing Nodes";
  print "*NODES" > d_file;
  mode = "NODE";
}

function do_nodes() {
  gsub(/^ */,"");
  print " ", $1, $2, $3, $4 > d_file;
  node_count++;
}

function start_elements() {
  print "Processing Elements";
  if( started_elements_output == 0 ) {
    print "*ELEMENTS" > d_file;
    started_elements_output = 1;
  } # else Just carry on adding element entries
  split( $2,elem_type,"=" );
  print "Element Type", elem_type[2];
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
  }
  dat_order = 2; # abaqus ordering

  mode = "ELEMENT";
}

function do_elements() {
  if( elem_type[2] == "C3D4" ) {
    gsub(/^ */,"");
    print " ", $1, "3 4 1", $2, $3, $4, $5, "1" > d_file;
    element_count++;
  } else if( elem_type[2] == "C3D8I" ) {
    gsub(/^ */,"");
    print " ", $1, "3 8 1", $2, $3, $4, $5, $6, $7, $8, $9, "1" > d_file;
    element_count++;
  } else if( elem_type[2] == "C3D8" ) {
    gsub(/^ */,"");
    print " ", $1, "3 8 1", $2, $3, $4, $5, $6, $7, $8, $9, "1" > d_file;
    element_count++;
  } else {
    print "Element Type", elem_type[2], "Not Supported";
  }
}

function start_nset() {
  split($2,nset_name,"=");
  node_set_id = node_set_count++;
  node_sets[nset_name[2]] = node_set_id; 
  print "Processing Nset", nset_name[2], node_sets[nset_name[2]], node_set_id;
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
  print "Generating NSET " node_set_id "[" $1 "," $2 "];" $3;
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
    print "Processing BC for Step 'Initial'";
  } else {
    print "Processing BC for Step " step_count;
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
    print "Restraining " mask " in " $1 ":" bc_nset_id " node set";
    for( nd in node_list ) { # append mask to bound_list
      print ">> " nd "/" bc_nset_id " " node_list[nd], match(node_list[nd],":" bc_nset_id ":")?" MATCH":" -";
      if( match(node_list[nd],":" bc_nset_id ":" ) ) {
        bound_list[nd] = bound_list[nd] mask;
      }
    }
  } else { # Fixed Displacement -steps
    if($2 == 1) {
      print "Fixing X in " $1 " node set";
      for( nd in node_list ) {
        if( match(node_list[nd],":" bc_nset_id ":") ) {
          fix_i_list[nd] = $4;
        }
      }
    }
    if(($2 <= 2) && ($3 >= 2)) {
      print "Loading Y in " $1 " node set";
      for( nd in node_list ) {
        if( match(node_list[nd],":" bc_nset_id ":") ) {
          fix_j_list[nd] = $4;
        }
      }
    }
    if($3 == 3) {
      print "Fixing Z in " $1 " node set";
      for( nd in node_list ) {
        if( match(node_list[nd],":" bc_nset_id ":") ) {
          fix_k_list[nd] = $4;
          print "--" fix_k_list[nd];
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
  print "Processing ELASTIC properties for material " mat_name[2] ;
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
    print "Processing CLOAD for Step 'Initial'";
  } else {
    print "Processing CLOAD for Step " step_count;
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
    print "Concetrated Force " mask " in " $1 " node set";
    for( nd in node_list ) { # append mask to bound_list
      if( match(node_list[nd],":" cload_nset_id ":") ) {
        bound_list[nd] = bound_list[nd] mask;
      }
    }
  } else { # load-steps; load IJK?
    print "STEP " step_count;
    if($2 == 1) {
      print "Concentrated Force X in " $1 " node set";
      for( nd in node_list ) {
        if( match(node_list[nd],":" cload_nset_id ":") ) {
          load_i_list[nd] = -$3;
        }
      }
    }
    if($2 == 2) {
      print "Concentrated Force Y in " $1 " node set";
      for( nd in node_list ) {
        if( match(node_list[nd],":" cload_nset_id ":") ) {
          load_j_list[nd] = -$3;
        }
      }
    }
    if($2 == 3) {
      print "Concentrated Force Z in " $1 " node set";
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
mode == "KEYWORD" && /^*[A-Z]/ { print NF, $0; }

# revert to KEYWORD mode if unhandled keyword found
mode != "KEYWORD" && /^*[A-Z]/ { mode = "KEYWORD"; }

# process non-keyword lines based on mode
mode == "NODE" { do_nodes(); }
mode == "ELEMENT" { do_elements(); }
mode == "NSET" { do_nset(); }
mode == "NSET-GEN" { do_nset_gen(); }
mode == "BOUNDARY" { do_bc(); }
mode == "ELASTIC" { do_elastic(); }
mode == "CLOAD" { do_cload(); }

# END OF MAIN
