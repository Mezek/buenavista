/**
 * $Date: 2013-05-30 14:45:21 +0200 (Thu, 30 May 2013) $
 * $Revision: 361 $
 * $Author: bartos $
 * $HeadURL: http://triglav/repos/BuenaVista/Nucleons/NucleonsDebugCS.cpp $
 * $Id: NucleonsDebugCS.cpp 361 2013-05-30 12:45:21Z bartos $
 *
 * @file
 * @brief	Debug file, v. with cross section data.
 */

/**
 * Some useful functions.
 */

/// Show values of FF in the point

int FFinPoint ( char* p, TComplex x )
{
	FFactor debugFF(12);
    debugFF.LoadParameters(p);
	std::cout << "\n t = " << x << std::endl;
	std::cout << "F1S: " << "\t" << debugFF.ScalarOne(x) << std::endl;
	std::cout << "F1V: " << "\t" << debugFF.VectorOne(x) << std::endl;
	std::cout << "F2S: " << "\t" << debugFF.ScalarTwo(x) << std::endl;
	std::cout << "F2V: " << "\t" << debugFF.VectorTwo(x) << std::endl;
	return 0;
}

using namespace ROOT::Minuit2;

/// Perform debug operations.

void performDebug ( char* p, char* f ) {

	std::cout << "\n> Debug parameters:   '" << p << "'" << std::endl;
	std::cout << "> Debug data output:  '" << f << "'" << std::endl;	
	std::cout << "\nDebug stuff:" << std::endl;


	/**
	std::cout << "> Values of data: " << std::endl;
	char dataFile[] = "dataNucleonsWithMainzGCyril.dat";	
	ExperimentalData A;
	A.ReadData(dataFile);
	A.ShowData();
	*/
	
	/**
	std::cout << "> Values of FF: " << std::endl;
	TComplex cA(-10.,0.);
	TComplex cB(-1.,0.);
	TComplex cC(10.,0.);
	FFinPoint ( p, cA );
	FFinPoint ( p, cB );
	FFinPoint ( p, k0 );
	FFinPoint ( p, k1 );
	FFinPoint ( p, cC );

	ofstream os (debugFile);

	time_t rawtime;
	time ( &rawtime );
	os << ctime (&rawtime) << "Performed debugging of values" << std::endl;

	os.close();
	*/
	
	/**
	std::cout << "> Cholesky decomposition test: " << std::endl;

	int nPar = 4;

	//Double_t a1[] = {18.,22.,54.,42.};
	//Double_t a2[] = {22.,70.,86.,62.};
	//Double_t a3[] = {54.,86.,174.,134.};
	//Double_t a4[] = {42.,62.,134.,106.};
	
	Double_t a1[] = {1.0,0.7,0.8,0.9};
	Double_t a2[] = {0.7,1.0,0.8,0.75};
	Double_t a3[] = {0.8,0.8,1.0,0.85};
	Double_t a4[] = {0.9,0.75,0.85,1.0};

	TVectorD b1; b1.Use(nPar,a1);
	TVectorD b2; b2.Use(nPar,a2);
	TVectorD b3; b3.Use(nPar,a3);
	TVectorD b4; b4.Use(nPar,a4);
	
	TMatrixD X(nPar,nPar);
	TMatrixDRow(X,0) = b1;
	TMatrixDRow(X,1) = b2;
	TMatrixDRow(X,2) = b3;
	TMatrixDRow(X,3) = b4;
	//X.Print();
	
	TMatrixDSym Y(nPar);
	Y.Use(0,3,a1);
	Y.Print();
	
	TDecompChol A(nPar,nPar);
	A.SetMatrix(Y);
	A.Decompose();

	TMatrixD B(nPar,nPar);
	B=A.GetU();
	B.Transpose(B);
	B.Print();
	*/

	std::cout << "\nEnd of debugging." << std::endl;

}
