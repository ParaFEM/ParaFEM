/* ---------------------------------------------------------------
 * metout2pf
 * ---------------------------------------------------------------
 * This program takes an unpartitioned ParaFEM input deck
 * (.d, .bnd, .lds, .fix) and a METIS output partition file and
 * outputs a new ParaFEM input deck (.d, .bnd, .lds, .fix AND
 * .psize) that incorporates the partitioning information.
 *
 * NOTE: if the BND, LDS, or FIX files are not present, they will
 * not be processed but the program will exit after processing the
 * model files.
 *
 * NOTE: the METIS partition programs (partnmesh/partdmesh) will
 * output TWO partition files, with the suffix of .epart.# or
 * .npart.#, where # matches the number of partitions specified
 * when the partition programs were executed. The epart.# file
 * should be used for this program.
 *
 * USAGE:
 *   metout2pf <partition_file> <input_basename> <output_basename>
 *
 *   e.g., metout2pf helix.met.epart.3 helix helix_part
 *
 * NOTE: use the basename i.e., without the .d suffix.
 *
 * WORKFLOW:
 *   e.g., d2metin helix.d helix.met
 *         partnmesh helix.met 3
 *         metout2pf helix.met.epart.3 helix helix_part

 * This will generate helix_part.d and helix_part.psize files. It
 * will also generate helix_part.bnd, helix_part.lds, and
 * helix_part.fix files if those input files were present.
 * ---------------------------------------------------------------
 * The code was originally developed by Vendel Szeremi as part of
 * his MSc dissertation "Scalable Parallel Finite Element Method".
 *
 * AUTHOR:
 *   Vendel Szeremi
 *
 * MODIFICATIONS:
     Louise M. Lever
     Mark Filipiak
 * ---------------------------------------------------------------
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

typedef struct node_s {
  long nn;
  long rnn;
  float x1;
  float x2;
  float x3;
} node;

typedef struct lds_s {
  long nn;
  float x1;
  float x2;
  float x3;
} lds;

typedef struct fix_s {
  long nn;
  long x1;
  float x2;
} fix;

int icmp(const void *pI1, const void *pI2);
int ionecmp(const void *pI1, const void *pI2);
int ldscmp(const void *pL1, const void *pL2);
int fixcmp(const void *pF1, const void *pF2);


int
main( int argc, char **argv ) {
  FILE *pPFH;
  FILE *pIFH;
  FILE *pOFH;
  char buf[256];
  long i;
  long j;
  long base;
  long idx;

  long iel;
  long ndim;
  long nod;
  long nn;
  long nnnew;
  long pnum;

  char fnamebuf[256];

  long *pPart;
  long maxPartNumber;
  long *pPartSize;
  node *pNodes;
  long *pSortNodes;
  long *pElem;

  long numlds;
  long numbnd;
  long numfix;
  lds *pLds;
  long *pBnd;
  fix *pFix;

  if (argc != 4) {
    printf("Usage: metout2pf partitionfile basenamein basenameout\n");
    printf(" Meshes with mixed element types not supported.\n");
    return(0);
  }

  printf("\r                                                                           ");
  printf("\n %04d process partition file", __LINE__);
  fflush(stdout);

  /* process partition file */

  /* pass 1: determine number of elements */
  pPFH = fopen(argv[1], "r");
  if (pPFH == NULL) {
    printf("Could not open %s\n.", argv[1]);
    return(0);
  }
  /* count number of lines; number of elements	 */
  iel = 0;
  while (fgets(buf, 256, pPFH)) {
    iel++;
  }
  fclose(pPFH);

  /* allocate 2-item (partition,id) array for all elements	 */
  pPart = malloc(iel * 2 * sizeof(long));
  /* error and exit on malloc fail	 */
  if (pPart == NULL) {
    printf("No mem.\n");
    return(0);
  }

  /* pass 2: read and build (partition,id) array	 */
  pPFH = fopen(argv[1], "r");
  if (pPFH == NULL) {
    printf("Could not open %s\n.", argv[1]);
    return(0);
  }
  i = 0;
  /* store partition number from file with element id as a pair for sorting */
  while (fgets(buf, 256, pPFH)) {
    sscanf(buf, "%ld", &pPart[i*2]);
    pPart[i*2+1] = i+1;
    i++;
  }
  fclose(pPFH);

  printf("\r                                                                           ");
  printf("\r %04d sort elements by partition", __LINE__);
  fflush(stdout);

  /* sort elements by the partition they belong to */
  qsort(pPart, iel, 2*sizeof(long), icmp);



  for (i=0; i<20; i++) {
    printf("DBG pPart i=%ld %ld %ld\n", i, pPart[i*2], pPart[i*2+1]);
  }

  printf("\r                                                                           ");
  printf("\r %04d find largest partition number", __LINE__);
  fflush(stdout);



  /* find largest partition number */
  /* LML: could move this inside pass 2 above	 */
  /* LML: partition list sorted, so last entry has highest number?	 */
  maxPartNumber = 0;
  for (i=0; i<iel; i++) {
    if (maxPartNumber < pPart[i*2]) {
      maxPartNumber = pPart[i*2];
    }
  }

  /* find sizes of partitions */
  pPartSize = malloc((maxPartNumber+1) * sizeof(long));
  if (pPartSize == NULL) {
    printf("No mem.\n");
    return(0);
  }
  /* zero the counts for each partition	 */
  bzero(pPartSize, (maxPartNumber+1) * sizeof(long));
  /* traverse the list and count number of elements per partition	 */
  for (i=0; i<iel; i++) {
    pPartSize[ pPart[i*2] ]++;
  }


  printf("\r                                                                           ");
  printf("\r %04d read nodes", __LINE__);
  fflush(stdout);



  /* open original model file (basename + .d) and read nodes */

  /* pass 1 - count number of nodes	 */
  strcpy(fnamebuf, argv[2]);
  strcat(fnamebuf, ".d");
  pIFH = fopen(fnamebuf, "r");
  if (pIFH == NULL) {
    printf("Could not open %s\n.", fnamebuf);
    return(0);
  }
  /* read lines until *NODES keyword found	 */
  while (1) {
    if (fgets(buf, 256, pIFH) == NULL) {
      printf("Error in input file, nodes not found.\n");
      fclose(pIFH);
      return(0);
    }
    buf[strlen(buf) - 1] = 0;
    if (strcmp(buf, "*NODES") == 0) {
      break;
    }
  }
  /* read lines and count number of nodes until *ELEMENTS keyword found  */
  nn = 0;
  while (1) {
    if (fgets(buf, 256, pIFH) == NULL) {
      printf("Error in input file, elements not found.\n");
      fclose(pIFH);
      return(0);
    }
    buf[strlen(buf) - 1] = 0;
    if (strcmp(buf, "*ELEMENTS") == 0) {
      break;
    }
    nn++;
  }
  fclose(pIFH);

  /* pass 2 - allocate 'node' (struct) array for num nodes	 */
  /* read node index and coordinates into array	 */
  pNodes = malloc(nn * sizeof(node));
  if (pNodes == NULL) {
    printf("No mem.\n");
    return(0);
  }
  bzero(pNodes, nn*sizeof(node));
/* read lines until *NODES keyword is found	 */
  pIFH = fopen(fnamebuf, "r");
  if (pIFH == NULL) {
    printf("Could not open %s\n.", argv[1]);
    return(0);
  }
  while (1) {
    if (fgets(buf, 256, pIFH) == NULL) {
      printf("Error in input file, nodes not found.\n");
      fclose(pIFH);
      return(0);
    }
    buf[strlen(buf) - 1] = 0;
    if (strcmp(buf, "*NODES") == 0) {
      break;
    }
  }
  /* read node lines until *ELEMENTS keyword is found	 */
  i = 0;
  while (1) {
    if (fgets(buf, 256, pIFH) == NULL) {
      printf("Error in input file, elements not found.\n");
      fclose(pIFH);
      return(0);
    }
    buf[strlen(buf) - 1] = 0;
    if (strcmp(buf, "*ELEMENTS") == 0) {
      break;
    }
    /* get the node ID and coordinate values	 */
    sscanf(buf, "%ld%f%f%f", &j, &pNodes[i].x1, &pNodes[i].x2, &pNodes[i].x3);
    pNodes[i].nn = j;
    i++;
  }



  printf("\r                                                                           ");
  printf("\r %04d read elements", __LINE__);
  fflush(stdout);




  /* now read elements - already know how many */
  if (fgets(buf, 256, pIFH) == NULL) {
    printf("Error in input file, no elements.\n");
    fclose(pIFH);
    return(0);
  }
  /* get element id, dim and type of FIRST element - assume all same type	 */
  sscanf(buf, "%ld%ld%ld", &j, &ndim, &nod);
  /* allocate array to store element connectivity for all elements of retrieved type	 */
  pElem = malloc((iel * (nod + 2)) * sizeof(long));
  if (pElem == NULL) {
    printf("No mem.\n");
    fclose(pIFH);
    return(0);
  }
  base = 0;
  /* read element lines until *DISPLACEMENTS keyword found or EOF	 */
  /* LML: *DISPLACEMENTS will NOT BE FOUND in ParaFEM model .d files	 */
  while (1) {
    sscanf(strtok(buf, " "), "%ld", &j);
    sscanf(strtok(NULL, " "), "%ld", &ndim);
    sscanf(strtok(NULL, " "), "%ld", &nod);
    strtok(NULL, " ");
    pElem[base] = j;
    for (i=0; i<nod; i++) {
      sscanf(strtok(NULL, " "), "%ld", &pElem[i+base+2]);
    }
    base += nod + 2;
    if (fgets(buf, 256, pIFH) == NULL) {
      break;
    }
    if (strcmp(buf, "*DISPLACEMENTS") == 0) {
      break;
    }
  }
  fclose(pIFH);


  printf("\r                                                                           ");
  printf("\r %04d renumber nodes", __LINE__);
  fflush(stdout);

  /* renumber nodes - new order from indices encountered in connectivity lists */
  nnnew = 1;
  for (i=0; i<iel; i++) {
    for (j=0; j<nod; j++) {
      idx = pElem[((pPart[i*2+1]-1)*(nod+2))+j+2] - 1;
      if (pNodes[idx].rnn == 0) {
	pNodes[idx].rnn = nnnew;
	nnnew++;
      }
    }
  }



  printf("\nDBG nnnew %ld\n", nnnew);

  j = 0;
  for (i=0; i<nn; i++) {
    if (pNodes[i].nn != i+1) {
      printf("DBG .nn != i+1\n");
    }
    if (pNodes[i].rnn == 0) {
      printf("DBG .rnn == 0\n");
    }
    if (pNodes[i].rnn > j) {
      j = pNodes[i].rnn;
    }
  }

  printf("DBG .rnn max %ld\n", j);

  printf("\r                                                                           ");
  printf("\r %04d write renumbered nodes", __LINE__);
  fflush(stdout);




  /* write renumbered nodes - create the new .d file */
  strcpy(fnamebuf, argv[3]);
  strcat(fnamebuf, ".d");
  pOFH = fopen(fnamebuf, "w");
  if (pOFH == NULL) {
    printf("Could not open %s.\n", fnamebuf);
  }
  /* output headers	 */
  fprintf(pOFH, "*THREE_DIMENSIONAL\n");
  fprintf(pOFH, "*NODES\n");

  /* build (node id,new node id) list for sorting	 */
  pSortNodes = malloc(2*nn*sizeof(long));
  if (pSortNodes == NULL) {
    printf("No mem.\n");
    return(0);
  }
  for (i=0; i<nn; i++) {
    pSortNodes[i*2] = pNodes[i].nn;
    pSortNodes[i*2+1] = pNodes[i].rnn;
  }
  /* sort based on new node id	 */
  qsort(pSortNodes, nn, 2*sizeof(long), ionecmp);
  /* output all nodes in new node id order	 */
  for( i=0; i<nn; i++ ) {
    fprintf(pOFH, " %ld %.8E %.8E %.8E\n", pSortNodes[i*2+1],
	    pNodes[pSortNodes[i*2]-1].x1, pNodes[pSortNodes[i*2]-1].x2,
	    pNodes[pSortNodes[i*2]-1].x3);
  }



  printf("\r                                                                           ");
  printf("\r %04d write renumbered elements", __LINE__);
  fflush(stdout);



  /* write renumbered elements */
  fprintf(pOFH, "*ELEMENTS\n");
  pnum = 1;
  for (i=0; i<iel; i++) {
    /* output new node id, ndim, nod and fixed "1" label	 */
    fprintf(pOFH, "%ld %ld %ld %ld ", i+1, ndim, nod, 1);
    /* output the connectivity indices using renumbered node ids	 */
    for (j=0; j<nod; j++) {
      fprintf(pOFH, "%ld ", pNodes[ pElem[((pPart[i*2+1]-1)*(nod+2))+j+2]-1 ].rnn);
    }
    /* output partition number as the material ID	 */
    /* LML: THIS PREVENTS ANY USE OF MULTIPLE MATERIALS	 */
    fprintf(pOFH, "%ld\n", pnum);
    if (pPart[i*2] != pPart[(i+1)*2]) {
      pnum++;
    }
  }
  fclose(pOFH);



  printf("\r                                                                           ");
  printf("\r %04d write psize file", __LINE__);
  fflush(stdout);




  /* write psize file - list of number of elements per partition */
  strcpy(fnamebuf, argv[3]);
  strcat(fnamebuf, ".psize");
  pOFH = fopen(fnamebuf, "w");
  if (pOFH == NULL) {
    printf("Could not open %s.\n", fnamebuf);
  }
  /* sanity check of number of partitions - ignores any empty sets	 */
  j = 0;
  for (i=0; i<=maxPartNumber; i++) {
    if (pPartSize[i] > 0) {
      j++;
    }
  }
  /* output number of partitions (non-zero counted)	 */
  fprintf(pOFH, "%ld\n", j);
  /* output size of non-empty partitions	 */
  for (i=0; i<=maxPartNumber; i++) {
    if (pPartSize[i] > 0) {
      fprintf(pOFH, "%ld ", pPartSize[i]);
    }
  }
  fclose(pOFH);



 loads:
  printf("\r                                                                           ");
  printf("\r %04d renumber loads", __LINE__);
  fflush(stdout);


  /* <<<----------------- LML: HERE --------------------->>>	 */



  /* renumber loads */
  strcpy(fnamebuf, argv[2]);
  strcat(fnamebuf, ".lds");
  pIFH = fopen(fnamebuf, "r");
  if (pIFH == NULL) {
    printf("\nCould not open %s, skipping\n", fnamebuf);
    fflush(stdout);
    goto restraints;
  }
  numlds = 0;
  while (fgets(buf, 256, pIFH)) {
    numlds++;
  }
  fclose(pIFH);

  pIFH = fopen(fnamebuf, "r");

  strcpy(fnamebuf, argv[3]);
  strcat(fnamebuf, ".lds");
  pOFH = fopen(fnamebuf, "w");
  if (pOFH == NULL) {
    printf("Could not open %s.\n", fnamebuf);
    return(0);
  }

  pLds = malloc(numlds*sizeof(lds));
  i = 0;
  while (fgets(buf, 256, pIFH)) {
    sscanf(buf, "%ld%f%f%f", &pLds[i].nn,
			&pLds[i].x1, &pLds[i].x2, &pLds[i].x3);
    i++;
  }
  for (i=0; i<numlds; i++) {
    pLds[i].nn = pNodes[pLds[i].nn-1].rnn;
  }
  qsort(pLds, numlds, sizeof(lds), ldscmp);
  for (i=0; i<numlds; i++) {
    fprintf(pOFH, " %ld %.8E %.8E %.8E\n",
	    pLds[i].nn, pLds[i].x1, pLds[i].x2, pLds[i].x3);
  }
  fclose(pIFH);
  fclose(pOFH);


 restraints:
  printf("\r                                                                           ");
  printf("\r %04d renumber restraints", __LINE__);
  fflush(stdout);

  /* renumber restraints */
  strcpy(fnamebuf, argv[2]);
  strcat(fnamebuf, ".bnd");
  pIFH = fopen(fnamebuf, "r");
  if (pIFH == NULL) {
    printf("\nCould not open %s, skipping\n", fnamebuf);
    fflush(stdout);
    goto fixed;
  }
  numbnd = 0;
  while (fgets(buf, 256, pIFH)) {
    numbnd++;
  }
  fclose(pIFH);

  pIFH = fopen(fnamebuf, "r");

  strcpy(fnamebuf, argv[3]);
  strcat(fnamebuf, ".bnd");
  pOFH = fopen(fnamebuf, "w");
  if (pOFH == NULL) {
    printf("Could not open %s.\n", fnamebuf);
    return(0);
  }
  pBnd = malloc(numbnd * 4 * sizeof(long));
  if (pBnd == NULL) {
    printf("No mem.\n");
    return(0);
  }
  i = 0;
  while (fgets(buf, 256, pIFH)) {
    sscanf(buf, "%ld%ld%ld%ld", &pBnd[i*4], &pBnd[i*4+1],
	   &pBnd[i*4+2], &pBnd[i*4+3]);
    i++;
  }
  for (i=0; i<numbnd; i++) {
    pBnd[i*4] = pNodes[pBnd[i*4]-1].rnn;
  }
  qsort(pBnd, numbnd, 4*sizeof(long), icmp);

  for (i=0; i<numbnd; i++) {
    fprintf(pOFH, "%ld %ld %ld %ld\n",
	    pBnd[i*4], pBnd[i*4+1], pBnd[i*4+2], pBnd[i*4+3]);
  }
  fclose(pIFH);
  fclose(pOFH);

 fixed:
  printf("\r                                                                           ");
  printf("\r %04d renumber fixed", __LINE__);
  fflush(stdout);

  /* renumber fixed */
  strcpy(fnamebuf, argv[2]);
  strcat(fnamebuf, ".fix");
  pIFH = fopen(fnamebuf, "r");
  if (pIFH == NULL) {
    printf("\nCould not open %s, skipping\n", fnamebuf);
    fflush(stdout);
    goto done;
  }
  numfix = 0;
  while (fgets(buf, 256, pIFH)) {
    numfix++;
  }
  fclose(pIFH);

  pIFH = fopen(fnamebuf, "r");

  strcpy(fnamebuf, argv[3]);
  strcat(fnamebuf, ".fix");
  pOFH = fopen(fnamebuf, "w");
  if (pOFH == NULL) {
    printf("Could not open %s.\n", fnamebuf);
    return(0);
  }

  pFix = malloc(numfix*sizeof(fix));
  i = 0;
  while (fgets(buf, 256, pIFH)) {
    sscanf(buf, "%ld%ld%f", &pFix[i].nn, &pFix[i].x1, &pFix[i].x2);
    i++;
  }
  for (i=0; i<numfix; i++) {
    pFix[i].nn = pNodes[pFix[i].nn-1].rnn;
  }
  qsort(pFix, numfix, sizeof(fix), fixcmp);
  for (i=0; i<numfix; i++) {
    fprintf(pOFH, " %ld %ld %.8E\n",
	    pFix[i].nn, pFix[i].x1, pFix[i].x2);
  }
  fclose(pIFH);
  fclose(pOFH);

 done:
  printf("\rdone                                      \n");
}

int icmp(const void *pI1, const void *pI2) {
  if (*(long*)pI1 < *(long*)pI2) {
    return(-1);
  }
  else if (*(long*)pI1 > *(long*)pI2) {
    return(1);
  }
  else {
    return(0);
  }
}

int ionecmp(const void *pI1, const void *pI2) {
  long *p1 = (long *)pI1;
  long *p2 = (long *)pI2;
  p1++;
  p2++;
  if (*p1 < *p2) {
    return(-1);
  }
  else if (*p1 > *p2) {
    return(1);
  }
  else {
    return(0);
  }
}

int ldscmp(const void *pL1, const void *pL2) {
  lds *p1 = (lds*)pL1;
  lds *p2 = (lds*)pL2;

  if (p1->nn < p2->nn) {
    return(-1);
  }
  else if (p1->nn > p2->nn) {
    return(1);
  }
  else {
    return(0);
  }
}

int fixcmp(const void *pF1, const void *pF2) {
  fix *p1 = (fix*)pF1;
  fix *p2 = (fix*)pF2;

  if (p1->nn < p2->nn) {
    return(-1);
  }
  else if (p1->nn > p2->nn) {
    return(1);
  }
  else {
    return(0);
  }
}
