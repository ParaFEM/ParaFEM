/* ---------------------------------------------------------------
 * metout2d
 * ---------------------------------------------------------------
 * This program takes an unpartitioned ParaFEM input deck
 * (.d, .bnd, .lds) and a METIS output partition file and outputs
 * a new ParaFEM input deck (.d, .bnd, .lds AND .psize) that
 * incorporates the partitioning information.
 *
 * NOTE: if the BND or LDS files are not present, they will not be
 * processed but the program will exit after processing the model
 * files.
 *
 * NOTE: the METIS partition programs (partnmesh/partdmesh) will
 * output TWO partition files, with the suffix of .epart.# or
 * .npart.#, where # matches the number of partitions specified
 * when the partition programs were executed. The epart.# file
 * should be used for this program.
 *
 * NOTE: ParaFEM .fix files are not supported.
 * 
 * USAGE:
 *   metout2d <partition_file> <input_basename> <output_basename>
 *
 *   e.g., metout2d helix.met.epart.3 helix helix_part
 *
 * NOTE: use the basename i.e., without the .d suffix.
 * 
 * WORKFLOW:
 *   e.g., d2metin helix.d helix.met
 *         partnmesh helix.met 3
 *         metout2d helix.met.epart.3 helix helix_part

 * This will generate helix_part.d and helix_part.psize files. It
 * will also generate helix_part.bnd and helix_part.lds if those
 * input files were present.
 * ---------------------------------------------------------------
 * The code was originally developed by Vendel Szeremi as part of
 * his MSc dissertation "Scalable Parallel Finite Element Method".
 *
 * AUTHOR:
 *   Vendel Szeremi
 *
 * MODIFICATIONS:
     Louise M. Lever
 * ---------------------------------------------------------------
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

typedef struct node_s {
	int nn;
	int rnn;
	float x1;
	float x2;
	float x3;
} node;

typedef struct lds_s {
	int nn;
	float x1;
	float x2;
	float x3;
} lds;


int icmp(const void *pI1, const void *pI2);
int ionecmp(const void *pI1, const void *pI2);
int ldscmp(const void *pL1, const void *pL2);


int main(int argc, char **argv) {
	FILE *pPFH;
	FILE *pIFH;
	FILE *pOFH;
	char buf[256];
	int i;
	int j;
	int base;
	int idx;

	int iel;
	int ndim;
	int nod;
	int nn;
	int nnnew;
	int pnum;

	char fnamebuf[256];

	int *pPart;
	int maxPartNumber;
	int *pPartSize;
	node *pNodes;
	int *pSortNodes;
	int *pElem;

	int numlds;
	int numbnd;
	lds *pLds;
	int *pBnd;

	if (argc != 4) {
		printf("Usage: metout2dan partitionfile basenamein basenameout\n");
		printf(" Meshes with mixed element types not supported.\n");
		return(0);
	}

printf("\r                                                                           ");
printf("\n %04d process partition file", __LINE__);
fflush(stdout);
	/* process partition file */
	pPFH = fopen(argv[1], "r");
	if (pPFH == NULL) {
		printf("Could not open %s\n.", argv[1]);
		return(0);
	}
	iel = 0;
	while (fgets(buf, 256, pPFH)) {
		iel++;
	}
	fclose(pPFH);
	pPart = malloc(iel * 2 * sizeof(int));
	if (pPart == NULL) {
		printf("No mem.\n");
		return(0);
	}
	pPFH = fopen(argv[1], "r");
	if (pPFH == NULL) {
		printf("Could not open %s\n.", argv[1]);
		return(0);
	}
	i = 0;
	while (fgets(buf, 256, pPFH)) {
		sscanf(buf, "%d", &pPart[i*2]);
		pPart[i*2+1] = i+1;
		i++;
	}
	fclose(pPFH);

printf("\r                                                                           ");
printf("\r %04d sort elements by partition", __LINE__);
fflush(stdout);
	/* sort elements by the partition they belong to */
	qsort(pPart, iel, 2*sizeof(int), icmp);

for (i=0; i<20; i++) {
	printf("DBG pPart i=%d %d %d\n", i, pPart[i*2], pPart[i*2+1]);
}

printf("\r                                                                           ");
printf("\r %04d find largest partition number", __LINE__);
fflush(stdout);
	/* find largest partition number */
	maxPartNumber = 0;
	for (i=0; i<iel; i++) {
		if (maxPartNumber < pPart[i*2]) {
			maxPartNumber = pPart[i*2];
		}
	}

	/* find sizes of partitions */
	pPartSize = malloc((maxPartNumber+1) * sizeof(int));
	if (pPartSize == NULL) {
		printf("No mem.\n");
		return(0);
	}
	bzero(pPartSize, (maxPartNumber+1) * sizeof(int));
	for (i=0; i<iel; i++) {
		pPartSize[ pPart[i*2] ]++;
	}


printf("\r                                                                           ");
printf("\r %04d read nodes", __LINE__);
fflush(stdout);
	/* read nodes */
	strcpy(fnamebuf, argv[2]);
	strcat(fnamebuf, ".d");
	pIFH = fopen(fnamebuf, "r");
	if (pIFH == NULL) {
		printf("Could not open %s\n.", fnamebuf);
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
	pNodes = malloc(nn * sizeof(node));
	if (pNodes == NULL) {
		printf("No mem.\n");
		return(0);
	}
	bzero(pNodes, nn*sizeof(node));
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
		sscanf(buf, "%d%f%f%f", &j, &pNodes[i].x1, &pNodes[i].x2, &pNodes[i].x3);
		pNodes[i].nn = j;
		i++;
	}

printf("\r                                                                           ");
printf("\r %04d read elements", __LINE__);
fflush(stdout);
	/* read elements */
	if (fgets(buf, 256, pIFH) == NULL) {
		printf("Error in input file, no elements.\n");
		fclose(pIFH);
		return(0);
	}
	sscanf(buf, "%d%d%d", &j, &ndim, &nod);
	pElem = malloc((iel * (nod + 2)) * sizeof(int));
	if (pElem == NULL) {
		printf("No mem.\n");
		fclose(pIFH);
		return(0);
	}
	base = 0;
	while (1) {
		sscanf(strtok(buf, " "), "%d", &j);
		sscanf(strtok(NULL, " "), "%d", &ndim);
		sscanf(strtok(NULL, " "), "%d", &nod);
		strtok(NULL, " ");
		pElem[base] = j;
		for (i=0; i<nod; i++) {
			sscanf(strtok(NULL, " "), "%d", &pElem[i+base+2]);
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
	/* renumber nodes */
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

printf("\nDBG nnnew %d\n", nnnew);
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
printf("DBG .rnn max %d\n", j);

printf("\r                                                                           ");
printf("\r %04d write renumbered nodes", __LINE__);
fflush(stdout);
	/* write renumbered nodes */
	strcpy(fnamebuf, argv[3]);
	strcat(fnamebuf, ".d");
	pOFH = fopen(fnamebuf, "w");
	if (pOFH == NULL) {
		printf("Could not open %s.\n", fnamebuf);
	}
	fprintf(pOFH, "*THREE_DIMENSIONAL\n");
	fprintf(pOFH, "*NODES\n");
	pSortNodes = malloc(2*nn*sizeof(int));
	if (pSortNodes == NULL) {
		printf("No mem.\n");
		return(0);
	}
	for (i=0; i<nn; i++) {
		pSortNodes[i*2] = pNodes[i].nn;
		pSortNodes[i*2+1] = pNodes[i].rnn;
	}
	qsort(pSortNodes, nn, 2*sizeof(int), ionecmp);
	for (i=0; i<nn; i++) {
		fprintf(pOFH, " %4d  % 1.4E  % 1.4E  % 1.4E\n", pSortNodes[i*2+1],
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
		fprintf(pOFH, "%4d %4d %4d %4d ", i+1, ndim, nod, 1);
		for (j=0; j<nod; j++) {
			fprintf(pOFH, "%4d ", pNodes[ pElem[((pPart[i*2+1]-1)*(nod+2))+j+2]-1 ].rnn);
		}
		fprintf(pOFH, "%4d\n", pnum);
		if (pPart[i*2] != pPart[(i+1)*2]) {
			pnum++;
		}
	}
	fclose(pOFH);

printf("\r                                                                           ");
printf("\r %04d write psize file", __LINE__);
fflush(stdout);
	/* write psize file */
	strcpy(fnamebuf, argv[3]);
	strcat(fnamebuf, ".psize");
	pOFH = fopen(fnamebuf, "w");
	if (pOFH == NULL) {
		printf("Could not open %s.\n", fnamebuf);
	}
	j = 0;
	for (i=0; i<=maxPartNumber; i++) {
		if (pPartSize[i] > 0) {
			j++;
		}
	}
	fprintf(pOFH, "%d\n", j);
	for (i=0; i<=maxPartNumber; i++) {
		if (pPartSize[i] > 0) {
			fprintf(pOFH, "%d ", pPartSize[i]);
		}
	}
	fclose(pOFH);


printf("\r                                                                           ");
printf("\r %04d renumber loads", __LINE__);
fflush(stdout);
	/* renumber loads */
	strcpy(fnamebuf, argv[2]);
	strcat(fnamebuf, ".lds");
	pIFH = fopen(fnamebuf, "r");
	if (pIFH == NULL) {
		printf("Could not open %s.\n", fnamebuf);
		return(0);
	}
	numlds = 0;
	while (fgets(buf, 256, pIFH)) {
		numlds++;
	}
	fclose(pIFH);

	pIFH = fopen(fnamebuf, "r");
	if (pIFH == NULL) {
		printf("Could not open %s.\n", fnamebuf);
		return(0);
	}
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
		sscanf(buf, "%d%f%f%f", &pLds[i].nn,
			&pLds[i].x1, &pLds[i].x2, &pLds[i].x3);
		i++;
	}
	for (i=0; i<numlds; i++) {
		pLds[i].nn = pNodes[pLds[i].nn-1].rnn;
	}
	qsort(pLds, numlds, sizeof(lds), ldscmp);
	for (i=0; i<numlds; i++) {
		fprintf(pOFH, " %4d  % 1.4E  % 1.4E  % 1.4E\n",
			pLds[i].nn, pLds[i].x1, pLds[i].x2, pLds[i].x3);
	}
	fclose(pIFH);
	fclose(pOFH);


printf("\r                                                                           ");
printf("\r %04d renumber restraints", __LINE__);
fflush(stdout);
	/* renumber restraints */
	strcpy(fnamebuf, argv[2]);
	strcat(fnamebuf, ".bnd");
	pIFH = fopen(fnamebuf, "r");
	if (pIFH == NULL) {
		printf("Could not open %s.\n", fnamebuf);
		return(0);
	}
	numbnd = 0;
	while (fgets(buf, 256, pIFH)) {
		numbnd++;
	}
	fclose(pIFH);

	pIFH = fopen(fnamebuf, "r");
	if (pIFH == NULL) {
		printf("Could not open %s.\n", fnamebuf);
		return(0);
	}
	strcpy(fnamebuf, argv[3]);
	strcat(fnamebuf, ".bnd");
	pOFH = fopen(fnamebuf, "w");
	if (pOFH == NULL) {
		printf("Could not open %s.\n", fnamebuf);
		return(0);
	}
	pBnd = malloc(numbnd * 4 * sizeof(int));
	if (pBnd == NULL) {
		printf("No mem.\n");
		return(0);
	}
	i = 0;
	while (fgets(buf, 256, pIFH)) {
		sscanf(buf, "%d%d%d%d", &pBnd[i*4], &pBnd[i*4+1],
			&pBnd[i*4+2], &pBnd[i*4+3]);
		i++;
	}
	for (i=0; i<numbnd; i++) {
		pBnd[i*4] = pNodes[pBnd[i*4]-1].rnn;
	}
	qsort(pBnd, numbnd, 4*sizeof(int), icmp);

	for (i=0; i<numbnd; i++) {
		fprintf(pOFH, "%4d   %d   %d   %d\n",
			pBnd[i*4], pBnd[i*4+1], pBnd[i*4+2], pBnd[i*4+3]);
	}

	fclose(pIFH);
	fclose(pOFH);

printf("\rdone                                      \n");
}

int icmp(const void *pI1, const void *pI2) {
	if (*(int*)pI1 < *(int*)pI2) {
		return(-1);
	}
	else if (*(int*)pI1 > *(int*)pI2) {
		return(1);
	}
	else {
		return(0);
	}
}

int ionecmp(const void *pI1, const void *pI2) {
	int *p1 = (int *)pI1;
	int *p2 = (int *)pI2;
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
