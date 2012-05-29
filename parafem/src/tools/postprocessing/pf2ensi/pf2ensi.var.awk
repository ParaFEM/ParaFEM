# @(#) pf2ensi.var.awk - Basic conversion of ParaFEM data file to EnSight Variable file
# @(#) Usage:awk -f pf2ensi.var.awk <filename>.{bnd,dis,fix,lds,nset,pri,rea,str,vms,ttr,flx}
# Author: Louise M. Lever (louise.lever@manchester.ac.uk)
# Version: 1.0.2
# Date: 2012-05-29

# CHANGES:
# v1.0.2
#   LML: Added support for .mat and element-scalar output
#   LML: Uses temporary file of element type and element id lists. Replaces id with loaded values
#        from .mat file into multiscalar ensi output files.
# v1.0.1
#   LML: Added support for .ttr and .flx (temperature-scalar, and flux-vector) files

BEGIN {
    RS = "\r\n|\n";
    OFS = " ";
    
    # default parsing mode - mode determied by file-type and/or header
    var_type = "NONE"; # NONE, NDBND, DISPL, NDFIX, NDLDS, NDSET, (PRIMAX, PRISC, PRIMI), NDREA, STRESS, NDVMS, NDTTR, NDFLX, ELMAT
    mode = "NONE"; # NONE, NODE, ELEMENT
    time_mode = "NONE"; # NONE, STEPS
    dtype = "NONE"; # NONE, SCALAR, MULTISCALAR, BND-SCALAR, VECTOR, TENSOR
    dset = "NONE"; # NONE, PARTIAL, FULL
    dopt = "NONE"; # NONE, BND, FIX, SET

    # --------------------------------------------------------------------------------
    # process input filename
    # --------------------------------------------------------------------------------

    vec_file = ARGV[1];
    # determine type and get extension position
    if( (file_ext_pos = match(vec_file,/.[dD][iI][sS]$/) - 1) != -1 ) {
	var_type = "DISPL";
	mode = "NODE";
	time_mode = "STEPS";
	dtype = "VECTOR";
	dset = "FULL";
    } else if( (file_ext_pos = match(vec_file,/.[pP][rR][iI]$/) - 1) != -1 ) {
	var_type = "PRIMX";
	var2_type = "PRISC";
	var3_type = "PRIMI";
	mode = "NODE";
	time_mode = "STEPS";
	dtype = "MULTISCALAR";
	dset = "FULL";
    } else if( (file_ext_pos = match(vec_file,/.[rR][eE][aA]$/) - 1) != -1 ) {
	var_type = "NDREA";
	mode = "NODE";
	time_mode = "STEPS";
	dtype = "VECTOR";
	dset = "FULL";
    } else if( (file_ext_pos = match(vec_file,/.[sS][tT][rR]$/) - 1) != -1 ) {
	var_type = "STRESS";
	mode = "NODE";
	time_mode = "STEPS";
	dtype = "TENSOR";
	dset = "FULL";
    } else if( (file_ext_pos = match(vec_file,/.[bB][nN][dD]$/) - 1) != -1 ) {
	var_type = "NDBND";
	mode = "NODE";
	dtype = "SCALAR";
	dset = "PARTIAL";
	dopt = "BND";
    } else if( (file_ext_pos = match(vec_file,/.[lL][dD][sS]$/) - 1) != -1 ) {
	var_type = "NDLDS";
	mode = "NODE";
	dtype = "VECTOR";
	dset = "PARTIAL";
    } else if( (file_ext_pos = match(vec_file,/.[vV][mM][sS]$/) - 1) != -1 ) {
	var_type = "NDVMS";
	mode = "NODE";
	time_mode = "STEPS";
	dtype = "SCALAR";
	dset = "FULL";
    } else if( (file_ext_pos = match(vec_file,/.[tT][tT][rR]$/) - 1) != -1 ) {
	var_type = "NDTTR";
	mode = "NODE";
	time_mode = "STEPS";
	dtype = "SCALAR";
	dset = "FULL";
    } else if( (file_ext_pos = match(vec_file,/.[fF][iI][xX]$/) - 1) != -1 ) {
	var_type = "NDFIX";
	mode = "NODE";
	dtype = "VECTOR";
	dset = "PARTIAL";
	dopt = "FIX";
    } else if( (file_ext_pos = match(vec_file,/.[fF][lL][xX]$/) - 1) != -1 ) {
	var_type = "NDFLX";
	mode = "NODE";
	time_mode = "STEPS";
	dtype = "VECTOR";
	dset = "FULL";
    } else if( (file_ext_pos = match(vec_file,/.[mM][aA][tT]$/) - 1) != -1 ) {
	var_type = "ELMAT";
	mode = "ELEMENT";
	time_mode = "NONE";
	dtype = "MULTISCALAR";
	dset = "FULL";
	dopt = "MAT";
    } else if( (file_ext_pos = match(vec_file,/.[nN][sS][eE][tT]$/) - 1) != -1 ) {
	var_type = "NDSET";
	mode = "NODE";
	dtype = "SCALAR";
	dset = "PARTIAL";
	dopt = "SET";
    }
    base_filename = substr(vec_file,0,file_ext_pos);
    
    # --------------------------------------------------------------------------------
    # generate output filenames
    # --------------------------------------------------------------------------------

    if( time_mode == "NONE" ) {
	if( var_type == "NDSET" ) {
	    var_base_filename = base_filename ".ensi." var_type "_";
	} else if( var_type == "ELMAT" ) {
	    var_base_filename = base_filename ".ensi." var_type "_";
	} else {
	    var_base_filename = base_filename ".ensi." var_type;
	}
    } else {
	var_base_filename = base_filename ".ensi." var_type "-";
	if( dtype == "MULTISCALAR" ) {
	    var2_base_filename = base_filename ".ensi." var2_type "-";
	    var3_base_filename = base_filename ".ensi." var3_type "-";
	}
    }

    # --------------------------------------------------------------------------------
    # generate temporary filenames
    # --------------------------------------------------------------------------------

    # partialid list
    if( dset == "PARTIAL" ) {
	tmp_pid_file = "tmp.pid." base_filename;
    }

    # scalar-1; i vector component; i tensor component
    tmp_icomp_file = "tmp.icomp." base_filename;
    # scalar-2, scalar-3 components; j,k vector components; tensor j,k components
    if( dtype == "VECTOR" || dtype == "TENSOR" || dtype == "MULTISCALAR" ) {
	tmp_jcomp_file = "tmp.jcomp." base_filename;
	tmp_kcomp_file = "tmp.kcomp." base_filename;
    }
    # tensor s,t,u components
    if( dtype == "TENSOR" ) {
	tmp_scomp_file = "tmp.scomp." base_filename;
	tmp_tcomp_file = "tmp.tcomp." base_filename;
	tmp_ucomp_file = "tmp.ucomp." base_filename;
    }

    # element template file - copy for element scalar data files and replace index with value
    elem_tmpl_file = "tmp." base_filename ".elem_tmpl"

    # --------------------------------------------------------------------------------
    # initialize static (and non-keyword based) datasets
    # --------------------------------------------------------------------------------

    if( var_type == "NDBND" ) {
    	start_bnd();
    } else if( var_type == "NDLDS" ) {
    	start_lds();
    } else if( var_type == "NDFIX" ) {
    	start_fix();
	cur_pid = -1;
    } else if( var_type == "NDSET" ) {
	init_nset();
	nset_count = 0;
	nset_expected = 0;
    }

    # initialize other vars
    time_step = 0;
    node_count = 0;
    elem_count = 0;
}

END {
    # end specific types 
    if( var_type == "NDFIX" ) {
	# for NDFIX, output last node
	do_fix_vector_partial_cur();
    } else if( var_type == "ELMAT" ) {
	end_mat_element_multiscalar()
    } else {
	# end generic types
	# complete last time step and close files
	if( dset == "PARTIAL" ) {
	    if( dtype == "SCALAR" ) {
		end_partial_scalar();
	    } else if( dtype == "VECTOR" ) {
		end_partial_vector();
	    }
	} else { # dset == FULL
	    if( dtype == "SCALAR" ) {
		end_scalar();
	    } else if( dtype == "VECTOR" ) {
		end_vector();
	    } else if( dtype == "TENSOR" ) {
		end_tensor();
	    } else if( dtype == "MULTISCALAR" ) {
		end_multi_scalar();
	    }
	}
    }

    print "Completed." > "/dev/stderr";

    # return via stdout to calling script:
    #   NSET names if NDSET file OR
    #   number of time steps for all other file types

    if( var_type == "NDSET" ) {
	print all_nset_names;
    } if( var_type == "ELMAT" ) {
	print all_prop_labels;
    } else {
	print time_step;
    }

    # All Done; EXIT
}

# --------------------------------------------------------------------------------
# FUNCTIONS
# --------------------------------------------------------------------------------

# --------------------------------------------------------------------------------
# Variable type specific functions:
# --------------------------------------------------------------------------------

# ---- BND ----------------------------------------------------------------------------

function start_bnd() {
    print "Processing BND (partial scalar) data" > "/dev/stderr";
    # no time steps

    # Open variable file for this partial scalar data; single time step
    var_file = var_base_filename;
    print "Generating variable file: " var_file > "/dev/stderr";

    print "Alya Ensight Gold --- Scalar per-partial-node variable file" > var_file;
    print "part" > var_file;
    print "1" > var_file;
    print "coordinates partial" > var_file;
}

function do_bnd_scalar_partial_data() {
    gsub(/^ */,"");
    print $1 > tmp_pid_file;
    print int(!$4*4 + !$3*2 + !$2) > tmp_icomp_file;
    node_count++;
}

# ---- DIS ----------------------------------------------------------------------------

function start_dis() {
    print "Processing DISPLACEMENT vector data" > "/dev/stderr";
    if( time_step > 0 ) {
	end_vector();
    }
    start_vector();
    getline; # LML: skip time step in file - assume increments by 1 each step
}

# ---- FIX ----------------------------------------------------------------------------

function start_fix() {
    print "Processing FIX (partial vector) data" > "/dev/stderr";
    # no time steps

    # Open variable file for this partial scalar data; single time step
    var_file = var_base_filename;
    print "Generating variable file: " var_file > "/dev/stderr";

    print "Alya Ensight Gold --- Vector per-partial-node variable file" > var_file;
    print "part" > var_file;
    print "1" > var_file;
    print "coordinates partial" > var_file;
}

function do_fix_vector_partial_cur() {
    if( cur_pid > 0 ) {
	print cur_pid > tmp_pid_file;
	print cur_i > tmp_icomp_file;
	print cur_j > tmp_jcomp_file;
	print cur_k > tmp_kcomp_file;
    }
    cur_i = 0.0;
    cur_j = 0.0;
    cur_k = 0.0;
}

function do_fix_vector_partial_data() {
    gsub(/^ */,"");
    if( $1 > cur_pid ) {
	do_fix_vector_partial_cur();
	node_count++;
	cur_pid = $1;
    }
    if( $2 == 1 ) {
	cur_i = $3;
    } else if( $2 == 2 ) {
	cur_j = $3;
    } else {
	cur_k = $3;
    }
}

# ---- LDS ----------------------------------------------------------------------------

function start_lds() {
    print "Processing LDS (partial vector) data" > "/dev/stderr";
    # no time steps

    # Open variable file for this partial scalar data; single time step
    var_file = var_base_filename;
    print "Generating variable file: " var_file > "/dev/stderr";

    print "Alya Ensight Gold --- Vector per-partial-node variable file" > var_file;
    print "part" > var_file;
    print "1" > var_file;
    print "coordinates partial" > var_file;
}

# ---- NSET ---------------------------------------------------------------------------

function init_nset() {
    getline;
    nset_problem = $1;
    getline;
    nset_expected = $1;
    print "Processing " nset_expected " NSETs for problem type: " nset_problem > "/dev/stderr";
}

function start_nset() {
    print "Processing NSET (partial scalar) data" > "/dev/stderr";
    # no time steps

    # get other values on *NSET line
    nset_ncount = $2;
    nset_name = $3;

    # end previous NSET if not first
    if( nset_count > 0 ) {
	end_partial_scalar();
	all_nset_names = all_nset_names " " nset_name;
    } else {
	all_nset_names = nset_name;
    }

    nset_count++;
    node_count = 0;

    # Open variable file for this partial scalar data; single time step
    var_file = var_base_filename nset_count "_" nset_name;
    print "Generating variable file: " var_file > "/dev/stderr";

    print "Alya Ensight Gold --- Scalar per-partial-node variable file" > var_file;
    print "part" > var_file;
    print "1" > var_file;
    print "coordinates partial" > var_file;
}

function do_nset_scalar_partial_data() {
    gsub(/^ */,"");
    print $1 > tmp_pid_file;
    print 1 > tmp_icomp_file;
    node_count++;
}

# ---- PRI ----------------------------------------------------------------------------

function start_pri() {
    print "Processing PRINCIPAL STRESS multiple scalar data" > "/dev/stderr";
    if( time_step > 0 ) {
	end_multi_scalar();
    }
    start_multi_scalar();
    getline; # LML: skip time step in file - assume increments by 1 each step
}

# ---- REA ----------------------------------------------------------------------------

function start_rea() {
    print "Processing NODAL REACTION vector data" > "/dev/stderr";
    if( time_step > 0 ) {
	end_vector();
    }
    start_vector();
    getline; # LML: skip time step in file - assume increments by 1 each step
}

# ---- STR ----------------------------------------------------------------------------

function start_str() {
    print "Processing STRESS tensor data" > "/dev/stderr";
    if( time_step > 0 ) {
	end_tensor();
    }
    start_tensor();
    getline; # LML: skip time step in file - assume increments by 1 each step
}

# ---- VMS ----------------------------------------------------------------------------

function start_vms() {
    print "Processing MISES STRESS scalar data (has dummy vector data)" > "/dev/stderr";
    if( time_step > 0 ) {
	end_scalar();
    }
    start_scalar();
    getline; # LML: skip time step in file - assume increments by 1 each step
}

# ---- TTR ----------------------------------------------------------------------------

function start_ttr() {
    print "Processing TEMPERATURE scalar data" > "/dev/stderr";
    if( time_step > 0 ) {
	end_scalar();
    }
    start_scalar();
    getline; # LML: skip time step in file - assume increments by 1 each step
}

# ---- FLX ----------------------------------------------------------------------------

function start_flx() {
    print "Processing FLUX vector data" > "/dev/stderr";
    if( time_step > 0 ) {
	end_vector();
    }
    start_vector();
    getline; # LML: skip time step in file - assume increments by 1 each step
}

# ---- MAT ----------------------------------------------------------------------------

function start_mat() {
    print "Processing MAT multiscalar (element) data" > "/dev/stderr";
    # no time steps

    # get other values on *MATERIAL line
    num_mats = $2;
    num_props = $3;

    # get labels
    getline;
    for( prop=1; prop<=num_props; prop++ ) {
	mat_prop_label[prop] = $prop;
	all_prop_labels = all_prop_labels " " $prop;
    }

    elem_count = 0;

    # Open variable file for this element scalar data; no time step
    for( prop=1; prop<=num_props; prop++ ) {
	var_file[prop] = var_base_filename prop "_" mat_prop_label[prop];
	print "Generating variable file: " var_file[prop] > "/dev/stderr";
	
	print "Alya Ensight Gold --- Scalar per-element variable file" > var_file[prop];
	print "part" > var_file[prop];
	print "1" > var_file[prop];
    }
}

function do_mat_element_multiscalar_data() {
    gsub(/^ */,"");
    pfld = 2;
    for( prop=1; prop<=num_props; prop++ ) {
	elem_id = $1;
	elem_prop[elem_id, prop] = $pfld;
	pfld++;
    }
    elem_count++;
}

function end_mat_element_multiscalar() {
    while( (getline < elem_tmpl_file) > 0 ) {
	if( ($1 == "tetra4") || ($1 == "tetra10") || ($1 == "hexa8") || ($1 == "hexa20") ) {
	    for( prop=1; prop<=num_props; prop++ ) {
		print $1 > var_file[prop];
	    }
	} else {
	    for( prop=1; prop<=num_props; prop++ ) {
		print elem_prop[$1, prop] > var_file[prop];
	    }
	}
    }
    for( prop=1; prop<=num_props; prop++ ) {
	close(var_file[prop]);
    }
}

# --------------------------------------------------------------------------------
# Data type specific functions:
# --------------------------------------------------------------------------------

# ---- Scalar --------------------------------------------------------------------

function start_scalar() {
    # incremement to next time step
    time_step++;

    # Open variable file for this step <name>.ensi.case and generate header
    var_file = sprintf( "%s%.6d", var_base_filename, time_step );
    print "Generating variable file: " var_file > "/dev/stderr";

    print "Alya Ensight Gold --- Scalar per-node variable file" > var_file
    print "part" > var_file;
    print "1" > var_file;
    print "coordinates" > var_file;
}

function do_scalar_data() {
    gsub(/^ */,"");
    print $2 > tmp_icomp_file;
    node_count++;
}

# See also:
# Variable specific function (above): do_bnd_scalar_partial_data()

function end_scalar() {
    # close the var_file and temporary files
    close(var_file);
    close(tmp_icomp_file);

    # append the scalar component to the var_file
    system("cat " tmp_icomp_file " >> " var_file);

    # clean up and remove the tmp files
    system("rm -f " tmp_icomp_file);
}

function end_partial_scalar() {
    # add the partial node count before node id list and data
    print node_count > var_file;
    # close the var_file and temporary files
    close(var_file);
    close(tmp_pid_file);
    close(tmp_icomp_file);

    # append each pid and scalar var to the var_file
    system("cat " tmp_pid_file " >> " var_file);
    system("cat " tmp_icomp_file " >> " var_file);

    # clean up and remove the tmp files
    system("rm -f " tmp_pid_file " " tmp_icomp_file);
}

# ---- Multi-Scalar ----------------------------------------------------------------------------

function start_multi_scalar() {
    # incremement to next time step
    time_step++;

    # Open variable file for this step <name>.ensi.case and generate header
    var_file = sprintf( "%s%.6d", var_base_filename, time_step );
    var2_file = sprintf( "%s%.6d", var2_base_filename, time_step );
    var3_file = sprintf( "%s%.6d", var3_base_filename, time_step );
    print "Generating variable file: " var_file > "/dev/stderr";
    print "Generating variable file: " var2_file > "/dev/stderr";
    print "Generating variable file: " var3_file > "/dev/stderr";

    print "Alya Ensight Gold --- Scalar per-node variable file" > var_file
    print "part" > var_file;
    print "1" > var_file;
    print "coordinates" > var_file;

    print "Alya Ensight Gold --- Scalar per-node variable file" > var2_file
    print "part" > var2_file;
    print "1" > var2_file;
    print "coordinates" > var2_file;

    print "Alya Ensight Gold --- Scalar per-node variable file" > var3_file
    print "part" > var3_file;
    print "1" > var3_file;
    print "coordinates" > var3_file;
}

function do_multi_scalar_data() {
    gsub(/^ */,"");
    print $2 > tmp_icomp_file;
    print $3 > tmp_jcomp_file;
    print $4 > tmp_kcomp_file;
    node_count++;
}

function end_multi_scalar() {
    # close the var_file and temporary files
    close(var_file);
    close(tmp_icomp_file);
    close(tmp_jcomp_file);
    close(tmp_kcomp_file);

    # append each scalar component to the varX_files
    system("cat " tmp_icomp_file " >> " var_file);
    system("cat " tmp_jcomp_file " >> " var2_file);
    system("cat " tmp_kcomp_file " >> " var3_file);

    # clean up and remove the tmp files
    system("rm -f " tmp_icomp_file " " tmp_jcomp_file " " tmp_kcomp_file);
}

# ---- Vector ----------------------------------------------------------------------------

function start_vector() {
    # incremement to next time step
    time_step++;

    # Open variable file for this step <name>.ensi.case and generate header
    var_file = sprintf( "%s%.6d", var_base_filename, time_step );
    print "Generating variable file: " var_file > "/dev/stderr";

    print "Alya Ensight Gold --- Vector per-node variable file" > var_file
    print "part" > var_file;
    print "1" > var_file;
    print "coordinates" > var_file;
}

function do_vector_data() {
    gsub(/^ */,"");
    print $2 > tmp_icomp_file;
    print $3 > tmp_jcomp_file;
    print $4 > tmp_kcomp_file;
    node_count++;
}

function do_vector_partial_data() {
    gsub(/^ */,"");
    print $1 > tmp_pid_file;
    print $2 > tmp_icomp_file;
    print $3 > tmp_jcomp_file;
    print $4 > tmp_kcomp_file;
    node_count++;
}

# See also:
# Variable specific function (above): function do_fix_vector_partial_cur()
# Variable specific function (above): function do_fix_vector_partial_data()

function end_vector() {
    # close the var_file and temporary files
    close(var_file);
    close(tmp_icomp_file);
    close(tmp_jcomp_file);
    close(tmp_kcomp_file);

    # append each vector component to the var_file
    system("cat " tmp_icomp_file " >> " var_file);
    system("cat " tmp_jcomp_file " >> " var_file);
    system("cat " tmp_kcomp_file " >> " var_file);

    # clean up and remove the tmp files
    system("rm -f " tmp_icomp_file " " tmp_jcomp_file " " tmp_kcomp_file);
}

function end_partial_vector() {
    # add the partial node count before node id list and data
    print node_count > var_file;
    # close the var_file and temporary files
    close(var_file);
    close(tmp_pid_file);
    close(tmp_icomp_file);
    close(tmp_jcomp_file);
    close(tmp_kcomp_file);

    # append each node id and vector component to the var_file
    system("cat " tmp_pid_file " >> " var_file);
    system("cat " tmp_icomp_file " >> " var_file);
    system("cat " tmp_jcomp_file " >> " var_file);
    system("cat " tmp_kcomp_file " >> " var_file);

    # clean up and remove the tmp files
    system("rm -f " tmp_pid_file " " tmp_icomp_file " " tmp_jcomp_file " " tmp_kcomp_file);
}

# ---- Tensor ----------------------------------------------------------------------------

function start_tensor() {
    # incremement to next time step
    time_step++;

    # Open variable file for this step <name>.ensi.case and generate header
    var_file = sprintf( "%s%.6d", var_base_filename, time_step );
    print "Generating variable file: " var_file > "/dev/stderr";

    print "Alya Ensight Gold --- Tensor symm per-node variable file" > var_file
    print "part" > var_file;
    print "1" > var_file;
    print "coordinates" > var_file;
}

function do_tensor_data() {
    gsub(/^ */,"");
    print $2 > tmp_icomp_file;
    print $3 > tmp_jcomp_file;
    print $4 > tmp_kcomp_file;
    print $5 > tmp_scomp_file;
    print $6 > tmp_tcomp_file;
    print $7 > tmp_ucomp_file;
    node_count++;
}

function end_tensor() {
    # close the var_file and temporary files
    close(var_file);
    close(tmp_icomp_file);
    close(tmp_jcomp_file);
    close(tmp_kcomp_file);
    close(tmp_scomp_file);
    close(tmp_tcomp_file);
    close(tmp_ucomp_file);

    # append each vector component to the var_file
    system("cat " tmp_icomp_file " >> " var_file);
    system("cat " tmp_jcomp_file " >> " var_file);
    system("cat " tmp_kcomp_file " >> " var_file);
    system("cat " tmp_scomp_file " >> " var_file);
    system("cat " tmp_tcomp_file " >> " var_file);
    system("cat " tmp_ucomp_file " >> " var_file);

    # clean up and remove the tmp files
    system("rm -f " tmp_icomp_file " " tmp_jcomp_file " " tmp_kcomp_file);
    system("rm -f " tmp_scomp_file " " tmp_tcomp_file " " tmp_ucomp_file);
}

function end_partial_tensor() {
 # NOT IMPLEMENTED YET
}

# --------------------------------------------------------------------------------
# Error reporting functions
# --------------------------------------------------------------------------------

function report_unknown_type() {
    print "Unexpected input and/or wrong file type." > "/dev/stderr";
    exit;
}

# --------------------------------------------------------------------------------
# MAIN RULES
# --------------------------------------------------------------------------------

# ignore comments and blank lines
/^#/ { next; }
/^$/ { next; }

# ignore NUM keywords
/^*NUM/ { next; }

# match these keywords and call appropriate start function
/^*DISPLACEMENT/ { start_dis(); next; }
/^*NODAL REACTIONS/ { start_rea(); next; }
/^*PRINCIPAL STRESS/ { start_pri(); next; }
/^*MISES STRESS/ { start_vms(); next; }
/^*STRESS/ { start_str(); next; }
/^*TEMPERATURE/ { start_ttr(); next; }
/^*FLUX/ { start_flx(); next; }
/^*NSET/ { start_nset(); next; }
/^*MATERIAL/ { start_mat(); next; }
# no keywords for BND, FIX and LDS - are manually "started" in BEGIN

# process non-keyword lines based on mode
mode == "NONE" { report_unknown_type(); }
# scalar inputs
mode == "NODE" && dtype == "SCALAR" && dset == "FULL" && dopt == "NONE" { do_scalar_data(); }
mode == "NODE" && dtype == "SCALAR" && dset == "PARTIAL" && dopt == "BND" { do_bnd_scalar_partial_data(); }
mode == "NODE" && dtype == "SCALAR" && dset == "PARTIAL" && dopt == "SET" { do_nset_scalar_partial_data(); }
# multiscalar inputs
mode == "NODE" && dtype == "MULTISCALAR" && dset == "FULL" && dopt == "NONE" { do_multi_scalar_data(); }
mode == "ELEMENT" && dtype == "MULTISCALAR" && dset == "FULL" && dopt == "MAT" { do_mat_element_multiscalar_data(); }
# vector inputs
mode == "NODE" && dtype == "VECTOR" && dset == "FULL" { do_vector_data(); }
mode == "NODE" && dtype == "VECTOR" && dset == "PARTIAL" && dopt == "NONE" { do_vector_partial_data(); }
mode == "NODE" && dtype == "VECTOR" && dset == "PARTIAL" && dopt == "FIX" { do_fix_vector_partial_data(); }
# tensor inputs
mode == "NODE" && dtype == "TENSOR" && dset == "FULL" { do_tensor_data(); }

# END OF MAIN
