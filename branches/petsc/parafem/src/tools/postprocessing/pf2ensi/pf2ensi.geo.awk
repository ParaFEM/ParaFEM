# @(#) pf2ensi.geo.awk - Basic conversion of ParaFEM .d model file to EnSight Model .geo
# @(#) Usage:awk -f [-v id_given=1] pf2ensi.geo.awk <filename.d> ;; as called by
# @9#) Usage:pf2ensi [-id-given] <filename>
# Author: Louise M. Lever (louise.lever@manchester.ac.uk)
# Version: 1.0.2
# Date: 2012-05-29

# CHANGES:
# v1.0.2
#   LML: Added support for .mat and element-scalar output
#   LML: Generates a temporary file with element type and element id lists. Used as
#        input to any Variable processing that outputs element data.
# v1.0.1
#   LML: Added support for .ttr and .flx (temperature-scalar, and flux-vector) files

# --------------------------------------------------------------------------------
# ARGUMENTS:
#   ARGV[1] : ParaFEM Model .d filename
#
# OPTIONAL PASSED VARIABLES:
#   flip_normals : (0) (default) to output as is;
#                  (1) to change element connectivity to flip normals
#   id_given     : (0) (default) to not produce id given lists
#                  (1) to produce id given lists for node and elements
# --------------------------------------------------------------------------------

# --------------------------------------------------------------------------------
# CODE SECTIONS:
#   BEGIN
#   END
#   Functions:
#     start_nodes()
#     do_nodes()
#     start_elements
#     do_elements()
#   MAIN RULES
# --------------------------------------------------------------------------------

# --------------------------------------------------------------------------------
# BEGIN
# --------------------------------------------------------------------------------
# 
# --------------------------------------------------------------------------------

BEGIN {
    RS = "\r\n|\n";
    OFS = " ";
    
    # default parsing mode (looking for keyword)
    mode = "KEYWORD";

    # flag initialization
    id_given = 0; # use node id given and element id given in GEO file

    # process input filename
    d_file = ARGV[1];
    file_ext_pos = match(d_file,/.[dD]$/) - 1;
    base_filename = substr(d_file,0,file_ext_pos);
    
    # generate output filenames
    geo_file = base_filename ".ensi.geo";
    matid_file = base_filename ".ensi.MATID";

    if( gen_elem_tmpl ) {
	elem_tmpl_file = "tmp." base_filename ".elem_tmpl"
    }

    # --------------------------------------------------------------------------------
    # generate temporary filenames
    # --------------------------------------------------------------------------------
    
    # node id given list
    if( id_given ) {
	tmp_idord_file = "tmp.idord." geo_file;
    }
    # xord, yord, zord lists
    tmp_xord_file = "tmp.xord." geo_file;
    tmp_yord_file = "tmp.yord." geo_file;
    tmp_zord_file = "tmp.zord." geo_file;
    # elem id given lists for each elem type
    if( id_given || gen_elem_tmpl ) {
	tmp_idtetra4_file = "tmp.idtetra4." geo_file;
	tmp_idtetra10_file = "tmp.idtetra10." geo_file;
	tmp_idhexa8_file = "tmp.idhexa8." geo_file;
	tmp_idhexa20_file = "tmp.idhexa20." geo_file;
    }
    # connectivity lists
    tmp_tetra4_file = "tmp.tetra4." geo_file;
    tmp_tetra10_file = "tmp.tetra10." geo_file;
    tmp_hexa8_file = "tmp.hexa8." geo_file;
    tmp_hexa20_file = "tmp.hexa20." geo_file;
    # material id lists
    tmp_matidtetra4_file = "tmp.matidtetra4." geo_file;
    tmp_matidtetra10_file = "tmp.matidtetra10." geo_file;
    tmp_matidhexa8_file = "tmp.matidhexa8." geo_file;
    tmp_matidhexa20_file = "tmp.matidhexa20." geo_file;
    # flags
    suppress_elem_errors = 0;
}

# --------------------------------------------------------------------------------
# END
# --------------------------------------------------------------------------------
# 
# --------------------------------------------------------------------------------

END {
    # close all files generated so far, so system commands can be called cleanly
    if( id_given ) {
	close(tmp_idord_file);
    }
    close(tmp_xord_file);
    close(tmp_yord_file);
    close(tmp_zord_file);
    if( id_given || gen_elem_tmpl ) {
	close(tmp_idtetra4_file);
	close(tmp_idhexa8_file);
	close(tmp_idtetra10_file);
	close(tmp_idhexa20_file);
    }
    close(tmp_tetra4_file);
    close(tmp_hexa8_file);
    close(tmp_tetra10_file);
    close(tmp_hexa20_file);
    close(tmp_matidtetra4_file);
    close(tmp_matidhexa8_file);
    close(tmp_matidtetra10_file);
    close(tmp_matidhexa20_file);
    
    # append the node id given count and list to the geo file
    system("echo " node_count " >> " geo_file);
    if( id_given ) {
	system("cat " tmp_idord_file " >> " geo_file);
    }
    # append each xord, yord, zord list to the geo file
    system("cat " tmp_xord_file " >> " geo_file);
    system("cat " tmp_yord_file " >> " geo_file);
    system("cat " tmp_zord_file " >> " geo_file);
    
    if( elem_tetra4_count > 0 ) {
	system("echo tetra4 >> " geo_file);
	system("echo tetra4 >> " matid_file);
	# append the tetra4 label to element template
	if( gen_elem_tmpl ) {
	    system("echo tetra4 >> " elem_tmpl_file);
	}
	# append the tetra4 elem count and list to the geo file
	system("echo " elem_tetra4_count " >> " geo_file);
	if( id_given ) {
	    system("cat " tmp_idtetra4_file " >> " geo_file);
	}
	# append the tetra4 id list to element template
	if( gen_elem_tmpl ) {
	    system("cat " tmp_idtetra4_file " >> " elem_tmpl_file);
	}
	# append the tetra4 connectivity to the geo file
	system("cat " tmp_tetra4_file " >> " geo_file);
	# append the tetra4 material list to the matid file
	system("cat " tmp_matidtetra4_file " >> " matid_file);
    }
    
    if( elem_tetra10_count > 0 ) {
	system("echo tetra10 >> " geo_file);
	system("echo tetra10 >> " matid_file);
	# append the tetra10 label to element template
	if( gen_elem_tmpl ) {
	    system("echo tetra10 >> " elem_tmpl_file);
	}
	# append the tetra10 elem count and list to the geo file
	system("echo " elem_tetra10_count " >> " geo_file);
	if( id_given ) {
	    system("cat " tmp_idtetra10_file " >> " geo_file);
	}
	# append the tetra10 id list to element template
	if( gen_elem_tmpl ) {
	    system("cat " tmp_idtetra10_file " >> " elem_tmpl_file);
	}
	# append the tetra10 connectivity to the geo file
	system("cat " tmp_tetra10_file " >> " geo_file);
	# append the tetra10 material list to the matid file
	system("cat " tmp_matidtetra10_file " >> " matid_file);
    }
    
    if( elem_hexa8_count > 0 ) {
	system("echo hexa8 >> " geo_file);
	system("echo hexa8 >> " matid_file);
	# append the hexa8 label to element template
	if( gen_elem_tmpl ) {
	    system("echo hexa8 >> " elem_tmpl_file);
	}
	# append the hexa8 elem count and list to the geo file
	system("echo " elem_hexa8_count " >> " geo_file);
	if( id_given ) {
	    system("cat " tmp_idhexa8_file " >> " geo_file);
	}
	# append the hexa8 id list to element template
	if( gen_elem_tmpl ) {
	    system("cat " tmp_idhexa8_file " >> " elem_tmpl_file);
	}
	# append the hexa8 connectivity to the geo file
	system("cat " tmp_hexa8_file " >> " geo_file);
	# append the hexa8 material list to the matid file
	system("cat " tmp_matidhexa8_file " >> " matid_file);
    }
    
    if( elem_hexa20_count > 0 ) {
	system("echo hexa20 >> " geo_file);
	system("echo hexa20 >> " matid_file);
	# append the hexa20 label to element template
	if( gen_elem_tmpl ) {
	    system("echo hexa20 >> " elem_tmpl_file);
	}
	# append the hexa20 elem count and list to the geo file
	system("echo " elem_hexa20_count " >> " geo_file);
	if( id_given ) {
	    system("cat " tmp_idhexa20_file " >> " geo_file);
	}
	# append the hexa20 id list to element template
	if( gen_elem_tmpl ) {
	    system("cat " tmp_idhexa20_file " >> " elem_tmpl_file);
	}
	# append the hexa20 connectivity to the geo file
	system("cat " tmp_hexa20_file " >> " geo_file);
	# append the hexa20 material list to the matid file
	system("cat " tmp_matidhexa20_file " >> " matid_file);
    }
    
    # clean up and remove the tmp files
    if( id_given ) {
	system("rm -f " tmp_idord_file);
    }
    if( id_given || gen_elem_tmpl ) {
	system("rm -f " tmp_idtetra4_file " " tmp_idtetra10_file " " tmp_idhexa8_file " " tmp_idhexa20_file);
    }
    system("rm -f " tmp_xord_file " " tmp_yord_file " " tmp_zord_file);
    system("rm -f " tmp_tetra4_file " " tmp_tetra10_file " " tmp_hexa8_file " " tmp_hexa20_file);
    system("rm -f " tmp_matidtetra4_file " " tmp_matidtetra10_file " " tmp_matidhexa8_file " " tmp_matidhexa20_file);
    
    print "Completed." > "/dev/stderr";

    # return data to calling script on stdout
    num_nodes = node_count;
    num_cells = elem_tetra4_count + elem_tetra10_count + elem_hexa8_count + elem_hex20_count;
    num_mats = length(mat_count);
    for( matidd in mat_count ) {
	mat_ids = mat_ids ":" matidd;
    }

    print num_nodes, num_cells, num_mats, mat_ids;

    # All Done; EXIT
}

# Functions

function start_nodes() {
    print "Processing Nodes" > "/dev/stderr";
    mode = "NODE";
}

function do_nodes() {
    gsub(/^ */,"");
    if( id_given ) {
	print $1 > tmp_idord_file;
    }
    print $2 > tmp_xord_file;
    print $3 > tmp_yord_file;
    print $4 > tmp_zord_file;
    node_count++;
}

function start_elements() {
    print "Processing Elements" > "/dev/stderr";
    mode = "ELEMENT";
}

function do_elements() {
    gsub(/^ */,"");
    
    if( $2 != "3" ) {
	if( !suppress_elem_errors ) {
	    print "2D elements not supported" > "/dev/stderr";
	    suppress_elem_errors = 1;
	}
    } else {
	if( !flip_normals ) {
	    if( $3 == "4" ) { # tetra4/tet2
		if( id_given || gen_elem_tmpl ) {
		    print $1 > tmp_idtetra4_file;
		}
		print $5, $6, $7, $8 > tmp_tetra4_file;
		print $9 > tmp_matidtetra4_file;
		mat_count[$9]++;
		elem_tetra4_count++;
	    } else if( $3 == "8" ) { #hexa8/hex
		if( id_given || gen_elem_tmpl ) {
		    print $1 > tmp_idhexa8_file;
		}
		print $5, $6, $7, $8, $9, $10, $11, $12 > tmp_hexa8_file;
		print $13 > tmp_matidhexa8_file;
		mat_count[$13]++;
		elem_hexa8_count++;
	    } else if( $3 == "10" ) { #tetra10/tet
		if( id_given || gen_elem_tmpl ) {
		    print $1 > tmp_idtetra10_file;
		}
		print $5, $6, $7, $8, $9, $10, $11, $12, $13, $14 > tmp_tetra10_file;
		print $15 > tmp_matidtetra10_file;
		mat_count[$15]++;
		elem_tetra10_count++;
	    } else if( $3 == "20" ) { #hexa20/hex2
		if( id_given || gen_elem_tmpl ) {
		    print $1 > tmp_idhexa20_file;
		}
		print $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $15, $16, $17, $18, $19, $20, $21, $22, $23, $24 > tmp_hexa20_file;
		print $25 > tmp_matidhexa20_file;
		mat_count[$25]++;
		elem_hexa20_count++;
	    }
	} else {
	    if( $3 == "4" ) { # FLIP tetra4/tet2 ;; { 0,_2_,_1_,3 }
		if( id_given || gen_elem_tmpl ) {
		    print $1 > tmp_idtetra4_file;
		}
		print $5, $7, $6, $8 > tmp_tetra4_file;
		print $9 > tmp_matidtetra4_file;
		mat_count[$9]++;
		elem_tetra4_count++;
	    }
	}
    }
}

# --------------------------------------------------------------------------------
# MAIN RULES
# --------------------------------------------------------------------------------

# ignore comments and blank lines
/^#/ { next; }
/^$/ { next; }

# match these keywords and enter appropriate mode
/^*NODES/ { start_nodes(); next; }
/^*ELEMENTS/ { start_elements(); next; }

# process non-keyword lines based on mode
mode == "NODE" { do_nodes(); }
mode == "ELEMENT" { do_elements(); }

# END OF MAIN
