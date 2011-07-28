/* --------------------------------------------------------------------------------
 * d2metin
 * VERSION: 1.0.1
 * --------------------------------------------------------------------------------
 * This program takes a ParaFEM input deck (.d file) and pre-processes it by
 * transforming it into the file format required for use by the METIS Partitioning
 * library and tools.
 *
 * USAGE:
 *   d2metin <input_model_filename> <output_model_filename>
 *
 *   e.g., d2metin helix.d helix.met
 *
 * Additional tools are then required to further process the output files e.g.,
 * partnmesh or partdmesh (from METIS).
 * --------------------------------------------------------------------------------
 * OVERVIEW:
 * (1) The code reads the ParaFEM input model file (*.d) and searches for the
 *     *ELEMENTS keyword. Only the ELEMENT connectivity data is required by METIS,
 *     so *NODES and all node coordinates are ignored.
 * (2) The dimensionality and number of nodes per element is read per node
 * --------------------------------------------------------------------------------
 * ISSUES:
 * (1) LML: Assumes non-mixed element type input; will fail if mixed is given; the
 *     "nod" value is set to whatever the last element read was.
 * (2) LML: Storage of connectivity data is overkill for this tool; the entire
 *     dynamic buffer code is not required; this code could easily read and write
 *     a single line at a time (assuming the non-mixed input holds true).
 * --------------------------------------------------------------------------------
 * The code was originally developed by Vendel Szeremi as part of his MSc
 * dissertation "Scalable Parallel Finite Element Method".
 *
 * AUTHOR:
 *   Vendel Szeremi
 *
 * MODIFICATIONS:
     Louise M. Lever
 * --------------------------------------------------------------------------------
 * CHANGES:
 * v1.0.1:
 *   Minor code formatting changes.
 *   Added support for partition mode flag IO (5th line of dat files).
 * --------------------------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define ELEM_BUF_SIZE 8192
typedef struct elemBuf_s {
  struct elemBuf_s *pNext;
  int numOcc;
  int elem[ELEM_BUF_SIZE];
} elemBuf;

int
main( int argc, char **argv ) {
  FILE *pIFH;
  FILE *pOFH;
  char buf[256];
  int i;
  elemBuf *pRootElemBuf;
  elemBuf *pElemBuf;
  int base;
  
  int iel;
  int ndim;
  int nod;
  
  /* Display USAGE if number of expected arguments incorrect 	 */
  if (argc != 3) {
    printf("Usage: dan2metin inputfile outputfile\n");
    printf(" Meshes with mixed element types not supported.\n");
    return(0);
  }

  /* open the existing ParaFEM .d model file and create the output METIS file	 */
  
  /* open input file - error and exit on fail */
  pIFH = fopen(argv[1], "r");
  if (pIFH == NULL) {
    printf("Could not open %s\n.", argv[1]);
    return(0);
  }
  
  /* open output file - error and exit on fail*/
  pOFH = fopen(argv[2], "w");
  if (pOFH == NULL) {
    printf("Could not open %s\n.", argv[2]);
    fclose(pIFH);
    return(0);
  }

  /* create first buffer root element - error and exit on malloc fail	 */

  pRootElemBuf = malloc(sizeof(elemBuf));
  if (pRootElemBuf == NULL) {
    printf("No memory.\n");
    fclose(pIFH);
    fclose(pOFH);
    return(0);
  }

  /* 	 */

  pElemBuf = pRootElemBuf;
  pElemBuf->pNext = NULL;
  pElemBuf->numOcc = 0;

  /* search for start of *ELEMENTS data - error and exit on unexpected EOF */
  
  while (1) {
    if( fgets(buf, 256, pIFH) == NULL ) {
      printf("Error in input file, elements not found.\n");
      fclose(pIFH);
      fclose(pOFH);
      return(0);
    }
    buf[strlen(buf) - 1] = 0;
    if (strcmp(buf, "*ELEMENTS") == 0) {
      break;
    }
  }

  /* loop over *ELEMENTS data - break on EOF or next label	 */
  /* LML: *DISPLACEMENTS label WILL NOT BE FOUND IN MODEL FILE	 */
  
  while (1) {
    if (fgets(buf, 256, pIFH) == NULL) {
      break;
    }
    if (strcmp(buf, "*DISPLACEMENTS") == 0) {
      break;
    }
    /* read ELEM-ID, NUM-DIMS and NODES-PER-ELEM tokens from line; skip the 4th token (1)	 */
    sscanf(strtok(buf, " "), "%d", &iel);
    sscanf(strtok(NULL, " "), "%d", &ndim);
    sscanf(strtok(NULL, " "), "%d", &nod);
    strtok(NULL, " ");

    /* grow the buffer array list if number of nodes exceeds remaining space in current buffer */
    if (pElemBuf->numOcc + nod > ELEM_BUF_SIZE) {
      /* allocate next buffer 	 */
      pElemBuf->pNext = malloc(sizeof(elemBuf));

      if (pElemBuf->pNext == NULL) {
	/* error and exit on malloc fail	 */
	printf("No memory.\n");
	fclose(pIFH);
	fclose(pOFH);
	return(0);
      }
      /* append new buffer to list	 */
      pElemBuf = pElemBuf->pNext;
      pElemBuf->numOcc = 0;
    }

    /* get current position in current buffer	 */
    base = pElemBuf->numOcc;

    /* read "nod" NODE INDICES from current line as required for element type and copy to current buffer */
    for (i=0; i<nod; i++) {
      sscanf(strtok(NULL, " "), "%d", &pElemBuf->elem[i+base]);
    }

    /* update buffer position to reflect new nodes added	 */
    pElemBuf->numOcc += nod;
  }


  
//VSZ
  printf("%d elements read / %d dim / %d nod\n", iel, ndim, nod);

  /* loop over buffer and output connectivity in METIS format	 */
  /* ONLY SUPPORTS Tet(nod=4) and Hex(nod=8) elements	 */

  /* Process HEX elements	 */
  if( nod == 8 ) {
    fprintf(pOFH, "%d 3\n", iel);
    pElemBuf = pRootElemBuf;
    base = 0;
    for (i=0; i<iel; i++) {
      fprintf(pOFH, "%7d %7d %7d %7d %7d %7d %7d %7d\n",
	      pElemBuf->elem[base+4], pElemBuf->elem[base+0],
	      pElemBuf->elem[base+3], pElemBuf->elem[base+7],
	      pElemBuf->elem[base+5], pElemBuf->elem[base+1],
	      pElemBuf->elem[base+2], pElemBuf->elem[base+6]);
      base += 8;
      if (pElemBuf->numOcc <= base) {
	base = 0;
	pElemBuf = pElemBuf->pNext;
      }
    }
  }
  /* Process TET elements	 */
  else if (nod == 4) {
    fprintf(pOFH, "%d 2\n", iel);
    pElemBuf = pRootElemBuf;
    base = 0;
    for (i=0; i<iel; i++) {
      fprintf(pOFH, "%7d %7d %7d %7d\n",
	      pElemBuf->elem[base+0], pElemBuf->elem[base+1],
	      pElemBuf->elem[base+1], pElemBuf->elem[base+2]);
      base += 4;
      if (pElemBuf->numOcc <= base) {
	base = 0;
	pElemBuf = pElemBuf->pNext;
      }
    }
  }
  /* Unsupported element type - error and exit */
  else {
    printf("%d nod not supported\n", nod);
    fclose(pIFH);
    fclose(pOFH);
    return(0);
  }

  /* Successful conversion - tody and exit	 */
  
  fclose(pIFH);
  fclose(pOFH);
  
  return(0);
}
