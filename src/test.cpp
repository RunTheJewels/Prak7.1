#include <iostream> 
#include <iomanip>
#include "smatrix.h"
#include "scalapack.h"


typedef complex<double> complex_d;

using namespace std;

SType filler(int i, int j) {
	i++;
	j++;
	return (i == j) ? i : (0.1*j+0.001*i);
}
SType fillerTest(int i, int j) {
	if (i == 0 && j == 0) {
		return 1;
	} else if (i == 0 && j == 4) {
		return 2;
	} else if (i == 1 && j == 2) {
		return 3;
	} else if (i == 3 && j == 1) {
		return 4;
	}
	return 0;
}

SType filler1(int i, int j)
{
	//SType mat[] = {SType(1,1), SType(2,0), SType(-1.5,3), SType(2,0), SType(2,-1),
		 //SType(4,-0.76), SType(-1.5,-3), SType(4,0.76), SType(-1,-1)};
		 //SType mat[] = {SType(1,0), SType(2,0), SType(3,0), SType(2,0), SType(-5,0),
		 //SType(4,0), SType(3,0), SType(4,0), SType(0,0)};
	SType mat[] = {SType(1.34,0), SType(0,0), SType(0,0), SType(0,0), SType(2,0),
		 SType(0,0), SType(0,0), SType(0,0), SType(-1,0)};
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
	SReal *S = x.calculateSVD(&SSize, &U, &VT);

	if (isRoot) {
		cout << endl;
		cout << "Producing SVD:" << endl;
		cout << " =" << endl;
	}
	cout << *U;
	if (isRoot) {
		cout << " *" << endl;
		cout << "  diag([";
		for (int i = 0; i < SSize; i++) {
			if (i > 0) {
				cout << " ";
			}
			cout << setprecision(8) << S[i];
		}
		cout << "])" << endl;
		cout << " *" << endl;
	}
	cout << *VT;
	
	delete U, VT, SSize; 
	delete[] S;
}

int main(int argc, char **argv) {
	MPI_Init(&argc, &argv);
	int size,rank;
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	
	bool isRoot = SMatrix::init(0, 0, 5);
	
	//////const int N_CONST = 3;
	//////int N = N_CONST;
	////////const int R = 3, C = 3;
	////////int ROWS = R, COLS = C;
	//////complex_d A[N_CONST*N_CONST] = {1,0,0,0,1,0,0,0,1}; //{1,2,3,2,-5,4,3,4,0};
	//////int ind = 0, LWORK=10*N*N, LRWORK=40*N, INFO;
	//////complex_d Z[N*N], WORK[LWORK], RWORK[LRWORK];
	//////double W[N];
	
	//////int brank,bsize,ictxt,myrow,mycol,nprow=1,npcol=1,zero=0,one=1,desc[9],info;
	//////Cblacs_pinfo(&brank, &bsize);
	//////Cblacs_get(-1,0,&ictxt);
	//////Cblacs_gridinit(&ictxt,(char*)"R",nprow,npcol);
	//////Cblacs_gridinfo(ictxt,&nprow,&npcol,&myrow,&mycol);
	//////int NB=(N-1)/nprow+1,MB=(N-1)/npcol+1,MXLLD=MB*NB;
	////////cout << NB << " " << MB << " " << zero << " " << ictxt << " " << MXLLD << endl;
	//////descinit_(desc,&N,&N,&NB,&MB,&zero,&zero,&ictxt,&MXLLD,&info);
	//////if (!rank) for (int i = 0; i < 9; i++)
	//////{
		//////cout << A[i] << " ";
		//////if (i % 3 == 2) cout << endl;
	//////}
	//////pzheev_((char*)"V",(char*)"U",&N,A,&one,&one,desc,W,Z,&one,&one,desc,WORK,&LWORK,RWORK,&LRWORK,&INFO);
	//////if (rank == 3) for (int i = 0; i < 9; i++)
	//////{
		//////cout << Z[i] << "\t";
		//////if (i % 3 == 2) cout << endl;
	//////}
	//////if (rank == 3) {for (int i = 0; i < 3; i++) cout << W[i] << " ";
	//cout << endl; cout << INFO << endl;}
	{
		SMatrix x(3, 3);
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
