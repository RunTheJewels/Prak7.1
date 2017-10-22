#include <iostream> 
#include <iomanip>
#include "smatrix.h"
#include "scalapack.h"


typedef complex<double> complex_d;

using namespace std;

const int SIZE = 3; //Размер матрицы

const int N = 20; //количество итераций

const SType deltat = SType(1,0); //Дельта-t, делённая на постоянную Планка

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
	if (i == 0 and j == 0) return 1;
	// if (i==0 or i ==2 or j==0 or j==2) return 0.5;
	return 0;
}

void printEvol(SMatrix &x, bool isRoot) {
	
	int SSize;
	
	SMatrix *U;
	SReal *S = x.calculateEigen(&SSize, &U, isRoot, &EV, deltat);
	
	if (isRoot) {
		cout << endl;
		cout << "Матрица собственных векторов = " << endl;
	}
	
	cout << *U;
	
	SMatrix eigenM(*U,SIZE,SIZE);
	eigenM.populate(&fillerEV); // Изначально eigenM хранит exp(-i * deltat * собств. значения) на диаголнальных элементах
	
	if (isRoot) {
		cout << endl;
		cout << "Матрица собственных значений = " << endl;
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
		&zero, eigenM.data, &intone, &intone, descb); // Теперь и далее eigenM хранит exp(-i * deltat * гамильтониан)


	SMatrix Ro(*U,SIZE,SIZE); // Матрица плотности
	int *descRo = Ro.getDesc();

	Ro.populate(&fillerRo);
	if (isRoot)	cout << "iter = 0:\n";
	cout << Ro;

	for (int i = 0; i < N; i++)
	{		
		res.populate(&fillerZero);
		
		pzgemm_((char *) "N", (char *) "N", &size, &size, &size, &one, eigenM.data, &intone, &intone, descb, Ro.data, 
		&intone, &intone, descRo, &zero, res.data, &intone, &intone, descc);
		
		pzgemm_((char *) "N", (char *) "C", &size, &size, &size, &one, res.data, &intone, &intone, descc, eigenM.data, 
		&intone, &intone, descb, &zero, Ro.data, &intone, &intone, descRo);
		
		if (isRoot)	cout << "iter = " << i+1 << ":\n";
		
		cout << Ro;
	}
	
	delete U, SSize; 
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
		printEvol(x, isRoot);
	}
	
	SMatrix::exit();
	return 0;
}
