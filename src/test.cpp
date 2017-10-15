#include <iostream> 
#include <iomanip>
#include "smatrix.h"
#include "scalapack.h"


typedef complex<double> complex_d;

using namespace std;

const int SIZE = 3; 

SType *EV;

// void mult(SMatrix **A, SMatrix **B, SMatrix **C) {
// 	SMatrix	*u = new SMatrix(*this, nRows, nMin),
// 		*vt = new SMatrix(*this, nMin, nCols);
// 	SReal *S = new SReal[nMin];
// 	SType *SC = new SType[nMin];
//         int	*desc = getDesc(),
// 		*descU = u->getDesc(),
// 		*descVT = vt->getDesc();

// 	SType *dataCopy = new SType[myProcSize];
// 	memcpy(dataCopy, data, myProcSize * sizeof(SType));

// 	SType *work = new SType[10 * nRows * nCols];

// 	int ione = 1, info, lwork = 10 * nRows * nCols;

// #ifdef USE_COMPLEX
// 	SType *rwork = new SType[1 + 4 * max(nRows, nCols)];
// 	int lrwork = 1 + 4 * max(nRows, nCols) * max(nRows, nCols);
// #endif

// 	pzheev_((char*) "V", (char*) "U", &nRows, dataCopy, &ione, &ione, 
// 		desc, S, u->data, &ione, &ione, descU, work, &lwork, rwork, 
// 		&lrwork, &info);
	
// #ifdef USE_COMPLEX
// 	lwork = (int) work[0].real();
// #else
// 	lwork = (int) work[0];
// #endif
	
// 	for (int i = 0; i < nMin; i++)
// 	{
// 		SC[i] = exp(S[i] * deltat * SType(0,-1));
// 	}

// 	delete[] work;
// 	work = new SType[lwork];
	
// 	delete[] work, dataCopy;
// #ifdef USE_COMPLEX
// 	delete[] rwork;
// #endif
// 	delete[] desc, descU, descVT;

// 	if (SSize != NULL) {
// 		*SSize = nMin;
// 	}

// 	if (U != NULL) {
// 		*U = u;
// 	} else {
// 		delete u;
// 	}

// 	if (VT != NULL) {
// 		*VT = vt;
// 	} else {
// 		delete vt;
// 	}

// 	*ev = SC;

// 	return S;
// }

SType fillerZero(int i, int j) {
	return 0;
}
SType fillerEV(int i, int j) {
	if (i == j) {
		return EV[i];
	}
	return 0;
}

SType filler1(int i, int j)
{
	//SType mat[] = {SType(1,1), SType(2,0), SType(-1.5,3), SType(2,0), SType(2,-1),
		 //SType(4,-0.76), SType(-1.5,-3), SType(4,0.76), SType(-1,-1)};
		 // SType mat[] = {SType(1,0), SType(2,0), SType(3,0), SType(2,0), SType(-5,0),
		 // SType(4,0), SType(3,0), SType(4,0), SType(0,0)};
	SType mat[] = {SType(2,0), SType(0,0), SType(0,0), SType(0,0), SType(1,0),
		 SType(3,0), SType(0,0), SType(3,0), SType(0,0)};
	return mat[i*3+j];
}

#ifdef USE_COMPLEX
SType filler0(int i, int j) {
	i++;
	j++;
	return SType(10 * i + j, 10 * j + i);
}
#endif

/// SVD demo
void printSVD(SMatrix &x, bool isRoot) {
	int SSize;
	SMatrix *U, *VT;
	SReal *S = x.calculateSVD(&SSize, &U, &VT, isRoot, &EV);
	if (isRoot) {
		cout << endl;
		cout << "Producing SVD:" << endl;
		cout << " =" << endl;
	}
	cout << *U;
	if (isRoot) cout << " *" << endl;
	SMatrix eigenM(*U,SIZE,SIZE);
	// eigenM.setIdentity();
	eigenM.populate(&fillerEV);
	cout << eigenM;
	if (isRoot) {
		// cout << " *" << endl;
		// cout << "  diag([";
		// for (int i = 0; i < SSize; i++) {
		// 	if (i > 0) {
		// 		cout << " ";
		// 	}
		// 	cout << setprecision(8) << S[i];
		// }
		// cout << "])" << endl;
		// cout << " *" << endl;
		// cout << EV << endl;
		// for (int i = 0; i < SSize; i++)
		// {
		// 	cout << EV[i] << " ";
		// }
		// cout << endl;
	}
	if (isRoot) cout << " *" << endl;
	// cout << *VT;
	
	SMatrix res(*U,SIZE,SIZE);
	// res.setIdentity();
	res.populate(&fillerZero);

	int size = SIZE, intone = 1;
	complex_d one = complex_d(1,0), zero = complex_d(0,0);

	int *desca = U->getDesc(), *descb = eigenM.getDesc(), *descc = res.getDesc();
	if (isRoot) cout << desca[6] << " " << descb[6] << " " << descc[6] << endl;
	pzgemm_((char *) "N", (char *) "N", &size, &size, &size, &one, U->data, &intone, &intone, desca, eigenM.data, 
		&intone, &intone, descb, &zero, res.data, &intone, &intone, descc);

	eigenM.populate(&fillerZero);

	pzgemm_((char *) "N", (char *) "N", &size, &size, &size, &one, res.data, &intone, &intone, descc, U->data, &intone, &intone, desca,
		&zero, eigenM.data, &intone, &intone, descb);

	cout << eigenM;

	delete U, VT, SSize; 
	delete[] S;
}

int main(int argc, char **argv) {
	MPI_Init(&argc, &argv);
	int size,rank;
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	
	bool isRoot = SMatrix::init(0, 0, 5);
	
	{
		SMatrix x(SIZE, SIZE);
		x.setIdentity();
		x.populate(&filler1);
		cout << x;
		printSVD(x, isRoot);
	}

	#ifdef USE_COMPLEX
	//SType we[] ={SType(2,1)};
	//cout << we[0];
	#endif
	SMatrix::exit();
	return 0;
}
