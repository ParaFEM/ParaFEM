/* ---------------------------------------------------------------
 * d2metin
 * ---------------------------------------------------------------
 * This program takes a ParaFEM input deck (.d file) and
 * pre-processes it by transforming it into the file format
 * required for use by the METIS Partitioning library and tools.
 *
 * USAGE:
 *   d2metin <input_model_filename> <output_model_filename>
 *
 *   e.g., d2metin helix.d helix.met
 *
 * Additional tools are then required to further process the
 * output files e.g., partnmesh or partdmesh (from METIS).
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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define ELEM_BUF_SIZE 8192
typedef struct elemBuf_s {
	struct elemBuf_s *pNext;
	int numOcc;
	int elem[ELEM_BUF_SIZE];
} elemBuf;



int main(int argc, char **argv) {
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


	if (argc != 3) {
		printf("Usage: dan2metin inputfile outputfile\n");
		printf(" Meshes with mixed element types not supported.\n");
		return(0);
	}

	/* open input file */
	pIFH = fopen(argv[1], "r");
	if (pIFH == NULL) {
		printf("Could not open %s\n.", argv[1]);
		return(0);
	}

	/* open output file */
	pOFH = fopen(argv[2], "w");
	if (pOFH == NULL) {
		printf("Could not open %s\n.", argv[2]);
		fclose(pIFH);
		return(0);
	}

	pRootElemBuf = malloc(sizeof(elemBuf));
	if (pRootElemBuf == NULL) {
		printf("No memory.\n");
		fclose(pIFH);
		fclose(pOFH);
		return(0);
	}
	pElemBuf = pRootElemBuf;
	pElemBuf->pNext = NULL;
	pElemBuf->numOcc = 0;

	while (1) {
		if (fgets(buf, 256, pIFH) == NULL) {
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

	while (1) {
		if (fgets(buf, 256, pIFH) == NULL) {
			break;
		}
		if (strcmp(buf, "*DISPLACEMENTS") == 0) {
			break;
		}
		sscanf(strtok(buf, " "), "%d", &iel);
		sscanf(strtok(NULL, " "), "%d", &ndim);
		sscanf(strtok(NULL, " "), "%d", &nod);
		strtok(NULL, " ");
		if (pElemBuf->numOcc + nod > ELEM_BUF_SIZE) {
			pElemBuf->pNext = malloc(sizeof(elemBuf));
			if (pElemBuf->pNext == NULL) {
				printf("No memory.\n");
				fclose(pIFH);
				fclose(pOFH);
				return(0);
			}
			pElemBuf = pElemBuf->pNext;
			pElemBuf->numOcc = 0;
		}
		base = pElemBuf->numOcc;
		for (i=0; i<nod; i++) {
			sscanf(strtok(NULL, " "), "%d", &pElemBuf->elem[i+base]);
		}
		pElemBuf->numOcc += nod;
	}

//VSZ
printf("%d elements read / %d dim / %d nod\n", iel, ndim, nod);

	if (nod == 8) {
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
	else {
		printf("%d nod not supported\n", nod);
		fclose(pIFH);
		fclose(pOFH);
		return(0);
	}

	fclose(pIFH);
	fclose(pOFH);

	return(0);

}
