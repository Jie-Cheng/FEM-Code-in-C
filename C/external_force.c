#include <math.h>
#include <petscksp.h>
#include "global_variables.h"
#include "external_force.h"
#include "helpers.h"
#include "shape.h"
#include "face.h"
#include "integration.h"

int ExternalPressure(Vec* dofs, Vec* fglo) {
	PetscErrorCode ierr;
	int i, j, k, intpt, row;
	// Initialize
	ierr = VecZeroEntries(*fglo); CHKERRQ(ierr);
	int nfacenodes = FaceNumNodes(nsd, nen); // The number of nodes on a face
	int npt = IntNum(nsd-1, nfacenodes, 0); // The number of integration points on a face
	double xilist[npt][nsd-1]; // The natural coordinates of the integration points
	IntPoints(nsd-1, nfacenodes, npt, xilist);
	double weights[npt]; // The weights of the integration points
	IntWeights(nsd-1, nfacenodes, npt, weights);

	double facecoord[nfacenodes][nsd];
	double facedof[nfacenodes][nsd];
	double fele[nsd*nfacenodes];

	for (k = 0; k < load_size; ++k) {
		for (i = 0; i < nsd*nfacenodes; ++i) fele[i] = 0.;
		int ele = load_num[k][0];
		int face = load_num[k][1];
		int nodelist[nfacenodes]; // The number (of order) of the nodes
		FaceNodes(nsd, nen, nfacenodes, face, nodelist);
		for (i = 0; i < nfacenodes; ++i) {
			for (j = 0; j < nsd; ++j) {
				int n = connect[ele-1][nodelist[i]-1]; // Global node number
				facecoord[i][j] = coords[n-1][j];
				int id = nsd*(n-1)+j;
				ierr = VecGetValues(*dofs, 1, &id, &facedof[i][j]); CHKERRQ(ierr);
				//facedof[i][j] = dofs[nsd*(n-1)+j];
			}
		}
		double external_pressure = load_val[k][0];
		for (intpt = 0; intpt < npt; ++intpt) {
			double dNdxi[nfacenodes][nsd-1];
			double N[nfacenodes];
			double dydxi[nsd][nsd-1];
			ShapeFun(nsd-1, xilist[intpt], nfacenodes, N);
			ShapeDer(nsd-1, xilist[intpt], nfacenodes, dNdxi);
			double temp[nfacenodes][nsd];
			MatMatAdd_(nfacenodes, nsd, facecoord, facedof, temp);
			MatTransposeMatMult_(nsd, nfacenodes, nsd-1, temp, dNdxi, dydxi); // (facecoord+facedof)^T * dNdxi
			double traction[nsd];
			if (nsd == 3) {
				double a[nsd], b[nsd];
				for (i = 0; i < nsd; ++i) {
					a[i] = dydxi[i][0];
					b[i] = dydxi[i][1];
				}
				Cross3d(a, b, traction);
			} else if (nsd == 2) {
				traction[0] =  dydxi[1][0];
				traction[1] = -dydxi[0][0];
			}
			// Compute the force
			for (i = 0; i < nfacenodes; ++i) {
				for (j = 0; j < nsd; ++j){
					row = nsd*i + j;
					fele[row] += external_pressure*traction[j]*N[i]*weights[intpt];
				}
			}
		}
		// scatter
		for (i = 0; i < nfacenodes; ++i) {
			for (j = 0; j < nsd; ++j) {
				int n = connect[ele-1][nodelist[i]-1];
				row = nsd*(n-1)+j;
				ierr = VecSetValue(*fglo, row, fele[nsd*i+j], ADD_VALUES); CHKERRQ(ierr);
				//fglo[row] += fele[nsd*i+j];
			}
		}
	}
	ierr = VecAssemblyBegin(*fglo); CHKERRQ(ierr);
	ierr = VecAssemblyEnd(*fglo); CHKERRQ(ierr);
	return 0;
}

int ExternalTraction(Vec* fglo) {
	PetscErrorCode ierr;
	int i, j, k, intpt, row;
	// Initialize
	ierr = VecZeroEntries(*fglo); CHKERRQ(ierr);

	int nfacenodes = FaceNumNodes(nsd, nen); // The number of nodes on a face
	int npt = IntNum(nsd-1, nfacenodes, 0); // The number of integration points on a face
	double xilist[npt][nsd-1]; // The natural coordinates of the integration points
	IntPoints(nsd-1, nfacenodes, npt, xilist);
	double weights[npt]; // The weights of the integration points
	IntWeights(nsd-1, nfacenodes, npt, weights);

	double facecoord[nfacenodes][nsd];
	double fele[nsd*nfacenodes];

	for (k = 0; k < load_size; ++k) {
		for (i = 0; i < nsd*nfacenodes; ++i) fele[i] = 0.;
		int ele = load_num[k][0];
		int face = load_num[k][1];
		int nodelist[nfacenodes]; // The number (of order) of the nodes
		FaceNodes(nsd, nen, nfacenodes, face, nodelist);
		for (i = 0; i < nfacenodes; ++i) {
			for (j = 0; j < nsd; ++j) {
				int n = connect[ele-1][nodelist[i]-1]; // Global node number
				facecoord[i][j] = coords[n-1][j];
			}
		}
		double traction[nsd];
		for (i = 0; i < nsd; ++i) traction[i] = load_val[k][i];
		for (intpt = 0; intpt < npt; ++intpt) {
			double dNdxi[nfacenodes][nsd-1];
			double N[nfacenodes];
			double dxdxi[nsd][nsd-1];
			ShapeFun(nsd-1, xilist[intpt], nfacenodes, N);
			ShapeDer(nsd-1, xilist[intpt], nfacenodes, dNdxi);
			MatTransposeMatMult_(nsd, nfacenodes, nsd-1, facecoord, dNdxi, dxdxi); // facecoords^T * dNdxi
			double det;
			if (nsd == 3) {
				det = sqrt( pow(dxdxi[1][0]*dxdxi[2][1] - dxdxi[1][1]*dxdxi[2][0], 2.) \
					+ pow(dxdxi[0][0]*dxdxi[2][1] - dxdxi[0][1]*dxdxi[2][0], 2.) \
					+ pow(dxdxi[0][0]*dxdxi[1][1] - dxdxi[0][1]*dxdxi[1][0], 2.));
			} else if (nsd == 2) {
				det = sqrt(pow(dxdxi[0][0], 2.) + pow(dxdxi[1][0], 2.));
			}
			// Compute the force
			for (i = 0; i < nfacenodes; ++i) {
				for (j = 0; j < nsd; ++j){
					row = nsd*i + j;
					fele[row] += traction[j]*N[i]*weights[intpt]*det;
				}
			}
		}
		// scatter
		for (i = 0; i < nfacenodes; ++i) {
			for (j = 0; j < nsd; ++j) {
				int n = connect[ele-1][nodelist[i]-1];
				row = nsd*(n-1)+j;
				ierr = VecSetValue(*fglo, row, fele[nsd*i+j], ADD_VALUES); CHKERRQ(ierr);
				//fglo[row] += fele[nsd*i+j];
			}
		}
	}
	ierr = VecAssemblyBegin(*fglo); CHKERRQ(ierr);
	ierr = VecAssemblyEnd(*fglo); CHKERRQ(ierr);
	return 0;
}

int ExternalGravity(Vec* fglo) {
	PetscErrorCode ierr;
	int i, j, ele, intpt, row;
	// Initialize
	ierr = VecZeroEntries(*fglo); CHKERRQ(ierr);
	double rho = materialprops[0];
	int npt = IntNum(nsd, nen, 0); // The number of integration points in an element
	double xilist[npt][nsd]; // The natural coordinates of the integration points
	IntPoints(nsd, nen, npt, xilist);
	double weights[npt]; // The weights of the integration points
	IntWeights(nsd, nen, npt, weights);

	double elecoord[nen][nsd];
	double fele[nsd*nen];

	for (ele = 0; ele < nel; ++ele) {
		for (i = 0; i < nsd*nen; ++i) fele[i] = 0.;
		for (i = 0; i < nen; ++i) {
			for (j = 0; j < nsd; ++j) {
				elecoord[i][j] = coords[connect[ele][i]-1][j];
			}
		}
		for (intpt = 0; intpt < npt; ++intpt) {
			double dNdxi[nen][nsd];
			double N[nen];
			double dxdxi[nsd][nsd];
			ShapeFun(nsd, xilist[intpt], nen, N);
			ShapeDer(nsd, xilist[intpt], nen, dNdxi);
			MatTransposeMatMult_(nsd, nen, nsd, elecoord, dNdxi, dxdxi); // elecoord^T * dNdxi
			double det = Determinant(nsd, dxdxi);
			// Compute the force
			for (i = 0; i < nen; ++i) {
				for (j = 0; j < nsd; ++j){
					row = nsd*i + j;
					fele[row] += gravity[j]*N[i]*weights[intpt]*det*rho;
				}
			}
		}
		// scatter
		for (i = 0; i < nen; ++i) {
			for (j = 0; j < nsd; ++j) {
				row = nsd*(connect[ele][i] - 1) + j;
				ierr = VecSetValue(*fglo, row, fele[nsd*i+j], ADD_VALUES); CHKERRQ(ierr);
				//fglo[row] += fele[nsd*i+j];
			}
		}
	}
	ierr = VecAssemblyBegin(*fglo); CHKERRQ(ierr);
	ierr = VecAssemblyEnd(*fglo); CHKERRQ(ierr);
	return 0;
}
