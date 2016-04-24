#include <stdio.h>
#include <stdlib.h>
#include "global_variables.h"
#include "read.h"

void ReadInput() {
	FILE *input;
	char *line = NULL;
	size_t len;
	int pos;
	int count = 0;
	int flg1 = 0;
	int flg2 = 0;
	int ints[6];
	double dous[14];

	input = fopen("input.txt", "r");
	if (!input) {
		printf("Cannot open file: input.txt!\n");
		exit(EXIT_FAILURE);
	}

	while (getline(&line, &len, input) != -1) {
		if (line[0] != '\n') { // Skip blank line
			for (pos = 0; pos < len; ++pos) {
				if (line[pos] == ':')
					break;
			}
			if (pos != len) {
				count++;
				int s, e, i;
				for (s = pos + 1; s < len; ++s) {
					if (line[s] != ' ' && line[s] != '\t')
						break;
				}
				for (e = s; e < len; ++e) {
					if (line[e] == ' ' || line[e] == '\t' \
						|| line[e] == '\n' || line[e] == '\0')
						break;
				}
				char buffer[50];
				if (count == 1 || count == 2 || count == 5 || \
					count == 11 || count == 13 || count == 15) {
					for (i = s; i < e; ++i) {
						buffer[i-s] = line[i];
					}
					buffer[i-s] = '\0';
					ints[flg1++] = atoi(buffer);
				} else {
					for (i = s; i < e; ++i) {
						buffer[i-s] = line[i];
						if (buffer[i-s] == 'd') {
							buffer[i-s] = 'e';
						}
					}
					buffer[i-s] = '\0';
					dous[flg2++] = atof(buffer);
				}
			}
		}	
	}
	fclose(input);
	if (line) free(line);

	mode = ints[0];
    isbinary = ints[1];
    tol = dous[0];
    penalty = dous[1];
    maxit = ints[2];
    gravity[0] = dous[2];
    gravity[1] = dous[3];
    gravity[2] = dous[4];
    inistep = dous[5];
    adjust = dous[6];
    nsteps = ints[3];
    dt = dous[7];
    nprint = ints[4];
    damp = dous[8];
    materialtype = ints[5];
    materialprops[0] = dous[9];
    materialprops[1] = dous[10];
    materialprops[2] = dous[11];
    materialprops[3] = dous[12];
    materialprops[4] = dous[13];
}

void ReadMesh() {
	int i, j;
	// coords
	FILE *incoords;
	incoords = fopen("coords.txt", "r");
	if (!incoords) {
		printf("Can not open file: coords.txt!\n");
		exit(EXIT_FAILURE);
	}
	fscanf(incoords, "%d %d", &nsd, &nn);

	coords = (double**)calloc(nn, sizeof(double*));
	for (i = 0; i < nn; ++i) coords[i] = (double*)calloc(nsd, sizeof(double));

	for (i = 0; i < nn; ++i) {
		for (j = 0; j < nsd; ++j) {
			fscanf(incoords, "%lf", &coords[i][j]);
		}
	}
	fclose(incoords);
	
	// connect
	FILE *inconnect;
	inconnect = fopen("connect.txt", "r");
	if (!inconnect) {
		printf("Can not open file: connect.txt!\n");
		exit(EXIT_FAILURE);
	}
	fscanf(inconnect, "%d %d", &nel, &nen);
	
	connect = (int**)calloc(nel, sizeof(int*));
	for (i = 0; i < nel; ++i) connect[i] = (int*)calloc(nen, sizeof(int));

	for (i = 0; i < nel; ++i) {
		for (j = 0; j < nen; ++j) {
			fscanf(inconnect, "%d", &connect[i][j]);
		}
	}
	fclose(inconnect);
	
	// bc
	FILE *inbc;
	inbc = fopen("bc.txt", "r");
	if (!inbc) {
		printf("Can not open file: bc.txt!\n");
		exit(EXIT_FAILURE);
	}
	fscanf(inbc, "%d", &bc_size);

	bc_num = (int**)calloc(bc_size, sizeof(int*));
	for (i = 0; i < bc_size; ++i) bc_num[i] = (int*)calloc(2, sizeof(int));
	bc_val = (double*)calloc(bc_size, sizeof(double));
	
	for (i = 0; i < bc_size; ++i) {
		fscanf(inbc, "%d %d %lf", &bc_num[i][0], &bc_num[i][1], &bc_val[i]);
	}
	fclose(inbc);

	// load
	FILE *inload;
	inload = fopen("load.txt", "r");
	if (!inload) {
		printf("Can not open file: load.txt!\n");
		exit(EXIT_FAILURE);
	}
	fscanf(inload, "%d %d", &load_size, &load_type);
	load_type -= 2;
	
	load_num = (int**)calloc(load_size, sizeof(int*));
	load_val = (double**)calloc(load_size, sizeof(double*));
	for (i = 0; i < load_size; ++i) {
		load_num[i] = (int*)calloc(2, sizeof(int));
		load_val[i] = (double*)calloc(load_type, sizeof(double));
	}
	
	for (i = 0; i < load_size; ++i) {
		for (j = 0; j < load_type + 2; ++j) {
			if (j < 2) fscanf(inload, "%d", &load_num[i][j]);
			else fscanf(inload, "%lf", &load_val[i][j-2]);
		}
	}
	fclose(inload);
	
	// share
	share = (int*)calloc(nn, sizeof(int));
	for (i = 0; i < nel; ++i) {
		for (j = 0; j < nen; ++j) {
			share[connect[i][j] - 1]++;
		}
	}
}

void Echo() {
	int i, j;

	// Echo the inputs
    printf("mode = %d\n", mode);
    printf("isbinary = %d\n", isbinary);
    printf("tol = %12.4e\n", tol);
    printf("penalty = %12.4e\n", penalty);
    printf("maxit = %d\n", maxit);
    printf("gravity = [%12.4e %12.4e %12.4e]\n", gravity[0], gravity[1], gravity[2]);
    printf("inistep = %12.4e\n", inistep);
    printf("adjust = %12.4e\n", adjust);
    printf("nsteps = %d\n", nsteps);
    printf("dt = %12.4e\n", dt);
    printf("nprint = %d\n", nprint);
    printf("damp = %12.4e\n", damp);
    printf("materialtype = %d\n", materialtype);
    printf("materialprops =[%12.4e %12.4e %12.4e %12.4e %12.4e]\n", \
     materialprops[0], materialprops[1], materialprops[2], materialprops[3], materialprops[4]);
    printf("nsd = %d\n", nsd);
    printf("nn = %d\n", nn);
    printf("nel = %d\n", nel);
    printf("nen = %d\n", nen);
    printf("bc_size = %d\n", bc_size);
	printf("load_size = %d\n", load_size);
	printf("load_type = %d\n", load_type);

	printf("Coordinates of the first 1 nodes\n");
	for (i = 0; i < 1; ++i) {
		for (j = 0; j < nsd; ++j) {
			printf("%f ", coords[i][j]);
		}
		printf("\n");
	}
	printf("\n");
	
	printf("Connect of the first 1 elements\n");
	for (i = 0; i < 1; ++i) {
		for (j = 0; j < nen; ++j) {
			printf("%d ", connect[i][j]);
		}
		printf("\n");
	}
	printf("\n");	

	printf("First 1 boundary conditions\n");
	for (i = 0; i < 1; ++i) {
		printf("%d %d %lf\n", bc_num[i][0], bc_num[i][1], bc_val[i]);
	}
	printf("\n");
	
	printf("First 1 loads\n");
	for (i = 0; i < 1; ++i) {
		for (j = 0; j < load_type + 2; ++j) {
			if (j < 2) printf("%d ", load_num[i][j]);
			else printf("%lf ", load_val[i][j-2]);
		}
		printf("\n");
	}
	printf("\n");
}

