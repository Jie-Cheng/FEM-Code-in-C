#include <stdio.h>
#include <math.h>
#include "read.h"
#include "helpers.h"
#include "material.h"
#include "internal_force.h"
#include "external_force.h"
#include "tangent_stiffness.h"

void Statics(int maxit, int isbinary, double inistep, double adjust, double tol, \
	double penalty, int materialtype, double materialprops[5], double gravity[3], \
	int nsd, int nen, int nn, int nel, double coords[nn][nsd], \
	int connect[nel][nen], int bc_size, int bc_num[bc_size][2], double bc_val[bc_size], \
	int load_size, int load_type, int load_num[load_size][2], \
	double load_val[load_size][load_type], int share[nn]);

extern void ReadInput();
extern void ReadMesh();
extern void ma57ds_(int* ndofs, double A[*ndofs][*ndofs], double b[*ndofs], double x[*ndofs]);

int main() {
	// Initialization
	int mode, maxit, nsteps, nprint, isbinary;
	int nsd, nen, nn, nel, materialtype, bc_size, load_size, load_type;
	double inistep, adjust, tol, dt, damp, penalty;
	double materialprops[5], gravity[3];
    ReadInput(&mode, &maxit, &nsteps, &nprint, &isbinary, &inistep, &adjust, &tol, &dt, \
    	&damp, &penalty, &materialtype, materialprops, gravity);
    double coords[nn][nsd];
    int connect[nel][nen];
    int bc_num[bc_size][2];
    double bc_val[bc_size];
    int load_num[load_size][2];
    double load_val[load_size][load_type];
    int share[nn];
    ReadMesh(&nsd, &nn, &nel, &nen, coords, connect, &bc_size, bc_num, bc_val, &load_size, &load_type, \
    	load_num, load_val, share);

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

    if (mode == 0) {
    	Statics(maxit, isbinary, inistep, adjust, tol, penalty, materialtype, materialprops, gravity, \
    		nsd, nen, nn, nel, coords, connect, bc_size, bc_num, bc_val, load_size, load_type, \
    		load_num, load_val, share);
    }

    return 0;
}


void Statics(int maxit, int isbinary, double inistep, double adjust, double tol, \
	double penalty, int materialtype, double materialprops[5], double gravity[3], \
	int nsd, int nen, int nn, int nel, double coords[nn][nsd], \
	int connect[nel][nen], int bc_size, int bc_num[bc_size][2], double bc_val[bc_size], \
	int load_size, int load_type, int load_num[load_size][2], \
	double load_val[load_size][load_type], int share[nn]) {

	double Fint[nn*nsd+nel];
	double Fext[nn*nsd+nel];
	double R[nn*nsd+nel];
	double w[nn*nsd+nel];
	double w1[nn*nsd+nel];
	double dw[nn*nsd+nel];
	double A[nn*nsd+nel][nn*nsd+nel];

	double constraint[bc_size];
	double der_constraint[bc_size][nn*nsd];

	int i, j, row, col;
	for (i = 0; i < nn*nsd+nel; ++i) w[i] = 0.0;
	for (i = 0; i < bc_size; ++i) {
		constraint[i] = 0.0;
		for (j = 0; j < nn*nsd; ++j) {
			der_constraint[i][j] = 0.0;
		}
	}

	double load_factor = 0.0;
	double increment = inistep;
	int step = 0;
	// output

	if (load_type != 1) // Traction
		ExternalTraction(nsd, nn, nel, nen, coords, connect, load_size, load_num, load_type, load_val, Fext);

	while (load_factor < 1.0) {
		step++;
		if (load_factor + increment > 1.0) {
			increment = 1.0 - load_factor;
		}
		for (i = 0; i < nn*nsd+nel; ++i) w1[i] = w[i];
		double err1 = 1.0;
		double err2 = 1.0;
		int nit = 0;
		printf("Step = %8d, LoadFactor = %12.4e\n", step, load_factor);
		while ((err1 > tol || err2 > tol) && nit < maxit) {
			nit++;
			InternalForce(nsd, nn, nel, nen, w, coords, connect, materialtype, materialprops, Fint);
			if (load_type == 1) // Pressure
				ExternalPressure(nsd, nn, nel, nen, w, coords, connect, load_size, load_num, load_type, load_val, Fext);
			InternalTangent(nsd, nn, nel, nen, w, coords, connect, materialtype, materialprops, A);
			// R = Fint - load_factor*Fext
			for (i = 0; i < nn*nsd+nel; ++i) {
				R[i] = Fint[i] - load_factor*Fext[i];
			} 
			for (i = 0; i < bc_size; ++i) {
				row = nsd*(bc_num[i][0] - 1) + bc_num[i][1] - 1;
				constraint[i] = w[row] - bc_val[i];
				der_constraint[i][row] = 1.0;
				A[row][row] += penalty*der_constraint[i][row]*der_constraint[i][row];
				R[row] += penalty*der_constraint[i][row]*constraint[i];
			}
			// MA57DS(nn*nsd+nel, A, -R, dw)
			for (i = 0; i < nn*nsd+nel; ++i) {
				R[i] = -R[i];
			}
			int n = nn*nsd + nel;
			ma57ds_(&n, A, R, dw);
			for (i = 0; i < nn*nsd+nel; ++i) w[i] += dw[i];
			err1 = sqrt(DotProduct(nn*nsd+nel, dw, dw)/DotProduct(nn*nsd+nel, w, w));
			err2 = sqrt(DotProduct(nn*nsd+nel, R, R)/(nsd*nn+nel));
			printf("Iteration = %8d     Err1 = %12.4e,     Err2 = %12.4e\n", nit, err1, err2);
		}
		if (nit == maxit) {
			for (i = 0; i < nn*nsd+nel; ++i) w[i] = w1[i];
			load_factor = load_factor - increment;
			increment /= 2.;
		} else if (nit > 6) {
			increment *= (1.0 - adjust);
		} else if (nit < 6) {
			increment *= (1.0 + adjust);
		}
	}
	// output
}