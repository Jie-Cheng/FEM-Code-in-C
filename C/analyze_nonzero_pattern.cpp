#include <map>
#include <stdlib.h>
#include "analyze_nonzero_pattern.h"
#include "global_variables.h"

void AnalyzeNonzerosPattern() {
	int ele, a, b, i, j, row, col, pos, ndofs;
	std::map<int, int> data;
	ndofs = nn*nsd + nel;
	nnz = 0;
	nonzeros = (int*)calloc(ndofs, sizeof(int));

	for (ele = 0; ele < nel; ++ele) {
		for (a = 0; a < nen; ++a) {
			for (i = 0; i < nsd; ++i) {
				row = nsd*(connect[ele][a] - 1) + i;
				for ( b = 0; b < nen; ++b) {
					for (j = 0; j < nsd; ++j) {
						col = nsd*(connect[ele][b] - 1) + j;
						pos = row*ndofs + col;
						data[pos]++;
						if (data[pos] == 1) nonzeros[row]++;
					}
				}
				col = nsd*nn + ele;
				pos = row*ndofs + col;
				data[pos]++;
				if (data[pos] == 1) nonzeros[row]++;
				pos = col*ndofs + row;
				data[pos]++;
				if (data[pos] == 1) nonzeros[col]++;
			}
		}
		row = nsd*nn + ele;
		col = nsd*nn + ele;
		pos = row*ndofs + col;
		data[pos]++;
		if (data[pos] == 1) nonzeros[row]++;
	}
	
	nnz = data.size();
}