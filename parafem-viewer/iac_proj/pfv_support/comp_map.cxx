
#define XP_WIDE_API	/* Use Wide APIs */

#include "iac_proj/pfv_support/gen.hxx"

#include<istream>
#include<iostream>
#include<fstream>
using namespace std;

/* --------------------------------------------------------------------------------
 * CXXMETHOD: component_mapper()
 * --------------------------------------------------------------------------------
 * This module method monitors an array of watched component selections i.e, the
 * selected component settings that various modules may use. Mappers that process
 * these input components may select which components are processed. If further
 * downstream Mappers need to process a different component that has not been
 * processed by the upstream Mappers then it will not be able to execute. Such
 * selection of components is employed to reduce memory and CPU requirements by
 * only processing what is needed rather than ALL components.
 *
 * An example is a FIELD containing 3 components (Principle Stress, Displacements
 * and Temperature). One module may access the Temperature component for slicing
 * and cutting while another module uses Displacements to render glyphs. If the
 * input to the glyph module is changed to use the slice output, then the Slice
 * Mapper must have processed BOTH required components for the glyph operation to
 * succeed.
 *
 * In Mappers, only the mapped components are output. This modifies the component
 * index in the array of components that a downstream Mappers sees. Given that a
 * user wishes to select from the ORIGINAL full set of components, a mapping is
 * needed from the original component index to the mapped component index.
 *
 * The example above, provides components 0, 1 and 2. Slice uses component 0 and
 * glyph uses component 2. If glyph is made to use the output of Slice, then if
 * components 0 AND 2 are processed, they will be output as components 0 AND 1.
 *
 * --------------------------------------------------------------------------------
 * INPUT:
 *   &GROUP monitor_comps[] {
 *     selected - The Mappers selected component index to process from ALL components.
 *     mapped - The components index in output of a Mapper.
 *
 * OUTPUT:
 *   reqd_comps[] - Array containing the selected component indexes from ALL components.
 *
 * ----------------------------------------------------------------------
 * 
 * --------------------------------------------------------------------------------
 */

int
component_mapper::update( OMevent_mask event_mask,
			  int seq_num ) {

  xp_long ncomps;
  unsigned int nc_loop;
  unsigned int max_comp_selected;
  int *mapped;
  unsigned int nmapped = 0;
  unsigned int map_loop;
  unsigned int reqd_index = 0;
  unsigned int *reqd_comps_arr;
  
  // get number of components we are monitoring
  monitor_comps.get_array_size( &ncomps );

  cout << endl;

  cout << "Number of monitored component selections: " << ncomps << endl;

  // find the maximum component index in use
  max_comp_selected = -1;
  for( nc_loop=0; nc_loop < ncomps; nc_loop++ ) {
    // reset mapped index
    monitor_comps[nc_loop].mapped = -1;
    cout << "[" << nc_loop << "] selected: " << monitor_comps[nc_loop].selected << " mapped: " << monitor_comps[nc_loop].mapped << endl;
    if( (unsigned int)(monitor_comps[nc_loop].selected) > max_comp_selected )
      max_comp_selected = monitor_comps[nc_loop].selected;
  }

  // allocat an array to say if input component is mapped or not
  mapped = new int[max_comp_selected + 1];

  // for each possible mapped component, check if selected and if so set as original index; -1 is not mapped
  for( map_loop = 0; map_loop <= max_comp_selected; map_loop++ ) {
    mapped[map_loop] = -1; // initialize as not mapped
    
    for( nc_loop=0; nc_loop < ncomps; nc_loop++ ) {
      if( (int)(monitor_comps[nc_loop].selected) == map_loop )
	mapped[map_loop] = (int)(monitor_comps[nc_loop].selected);
    }
  }

  // count total number of mapped components
  for( map_loop = 0; map_loop <= max_comp_selected; map_loop++ ) {
    cout << "Component mapped: " << mapped[map_loop] << endl;
    if( mapped[map_loop] >= 0 ) {
      nmapped++;
    }
  }
  cout << "Number of mapped components: " << nmapped << endl;

  // allocated array for specifying mapped components
  reqd_comps.set_array_size( nmapped );
  reqd_comps_arr = (unsigned int *)reqd_comps.ret_array_ptr( OM_GET_ARRAY_WR, NULL, NULL );
  
  // populate array with indices of original selected component indices
  for( map_loop=0; map_loop < max_comp_selected; map_loop++ ) {
    if( mapped[map_loop] >= 0 ) {
      reqd_comps_arr[reqd_index] = mapped[map_loop];
      monitor_comps[map_loop].mapped = reqd_index++;
    }
  }

  cout << endl;
  
  ARRfree( reqd_comps_arr );
  
  return(1);
}
