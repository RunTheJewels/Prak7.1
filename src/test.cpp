#include <iostream> 
#include <iomanip>
#include "smatrix.h"
#include "scalapack.h"


typedef complex<double> complex_d;

using namespace std;

const int SIZE = 3; 

const int N = 20;

SType *EV;

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
	
	SType mat[] = {SType(2,0), SType(0,0), SType(0,0),
		 		   SType(0,0), SType(1,0), SType(3,0), 
		 		   SType(0,0), SType(3,0), SType(0,0)};
	
	return mat[i*3+j];
}

SType fillerRo(int i, int j) {
	if (i == 1 and j == 1) return 1;
	//if (i==0 or i ==2 or j==0 or j==2) return 0.5;
	return 0;
}

/// SVD demo
void printSVD(SMatrix &x, bool isRoot, int rank) {
	
	int SSize;
	
	SMatrix *U, *VT;
	SReal *S = x.calculateSVD(&SSize, &U, &VT, isRoot, &EV);
	
	if (isRoot) {
		cout << endl;
		cout << "Матрица собственных векторов:" << endl;
		cout << " =" << endl;
	}
	
	cout << *U;
	
	SMatrix eigenM(*U,SIZE,SIZE);
	eigenM.populate(&fillerEV);
	
	if (isRoot) {
		cout << endl;
		cout << "Матрица собственных значений:" << endl;
		cout << " =" << endl;
	}
	
	cout << eigenM;

	if (isRoot) cout << endl;
	
	SMatrix res(*U,SIZE,SIZE);
	res.populate(&fillerZero);

	int size = SIZE, intone = 1;
	complex_d one = complex_d(1,0), zero = complex_d(0,0);

	int *desca = U->getDesc(), *descb = eigenM.getDesc(), *descc = res.getDesc();

	pzgemm_((char *) "N", (char *) "N", &size, &size, &size, &one, U->data, &intone, &intone, desca, eigenM.data, 
		&intone, &intone, descb, &zero, res.data, &intone, &intone, descc);

	eigenM.populate(&fillerZero);

	pzgemm_((char *) "N", (char *) "C", &size, &size, &size, &one, res.data, &intone, &intone, descc, U->data, &intone, &intone, desca,
		&zero, eigenM.data, &intone, &intone, descb);


	SMatrix Ro(*U,SIZE,SIZE);
	int *descRo = Ro.getDesc();

	Ro.populate(&fillerRo);
	if (isRoot)	cout << "iter = 0:\n" << Ro;

	for (int i = 0; i < N; i++)
	{
		Ro.barrier();
		res.barrier();
		eigenM.barrier();
		MPI_Barrier(MPI_COMM_WORLD);
		
		res.populate(&fillerZero);
		
		pzgemm_((char *) "N", (char *) "N", &size, &size, &size, &one, eigenM.data, &intone, &intone, descb, Ro.data, 
		&intone, &intone, descRo, &zero, res.data, &intone, &intone, descc);
		
		Ro.barrier();
		res.barrier();
		eigenM.barrier();
		MPI_Barrier(MPI_COMM_WORLD);
		
		pzgemm_((char *) "N", (char *) "C", &size, &size, &size, &one, res.data, &intone, &intone, descc, eigenM.data, 
		&intone, &intone, descb, &zero, Ro.data, &intone, &intone, descRo);
		
		if (isRoot)	cout << "iter = " << i+1 << ":\n";
		
		
		Ro.barrier();
		res.barrier();
		eigenM.barrier();
		MPI_Barrier(MPI_COMM_WORLD);
		
		cout << Ro;
		
		Ro.barrier();
		res.barrier();
		eigenM.barrier();
		
		MPI_Barrier(MPI_COMM_WORLD);
	}
	
		cout << rank << endl;
	
		Ro.barrier();
		res.barrier();
		eigenM.barrier();
		MPI_Barrier(MPI_COMM_WORLD);
	
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
		printSVD(x, isRoot,rank);
	}
	
	SMatrix::exit();
	return 0;
}
