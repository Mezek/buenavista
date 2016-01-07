/**
 * $Date$
 * $Revision$
 * $Author$
 *
 * @file
 * @brief	Debug file.
 */

using namespace ROOT::Minuit2;

/// Perform debug operations.

void performDebug ( char* p, char* f ) {

	std::cout << "\nDebug stuff:" << std::endl;
	std::cout << "> Debug parameters:   `" << p << "'" << std::endl;
	std::cout << "> Debug data output:  `" << f << "'" << std::endl;	

	/*
	// Test argument parsing

	std::cout << "Arguments:" << std::endl;
	std::cout << "-f: " << globalArgs.fit << std::endl;
	std::cout << "-p: " << globalArgs.parameters << std::endl;
	std::cout << "-o: " << globalArgs.output << std::endl;
	std::cout << "-d: " << globalArgs.data << std::endl;
	std::cout << "-v: " << globalArgs.verbose << std::endl;*/

	// Convert parameters
	double q1 = 2.2181;
	double q2 = 2.2203;
	double q3 = 6.0527;
	double q4 = 5.3767;
	std::cout << ">> 1s : "<< (q1*q1+1.)*t0s << std::endl;
	std::cout << ">> 1v : "<< (q3*q3+1.)*t0v << std::endl;
	std::cout << ">> 2s : "<< (q2*q2+1.)*t0s << std::endl;
	std::cout << ">> 2v : "<< (q4*q4+1.)*t0v << std::endl;
	

	/*// Inverse Convert parameters
	double t1 = 2.01907;
	double t2 = 2.52493;
	double t3 = 1.83342;
	double t4 = 2.31511;
	std::cout << ">> 1s : "<< sqrt(t1/t0s-1.) << std::endl;
	std::cout << ">> 1v : "<< sqrt(t2/t0v-1.) << std::endl;
	std::cout << ">> 2s : "<< sqrt(t3/t0s-1.) << std::endl;
	std::cout << ">> 2v : "<< sqrt(t4/t0v-1.) << std::endl;*/

	FFactor dM(12);
    dM.LoadParameters(p);

	double r = 0.1;
	for (int i = 0; i < 6; ++i) {
		double s = 1.;
		for (int j = 0; j < i+1; ++j) {
			s = s*r;
		}
		std::cout << "> Proton radius: " << dM.RadiusEP(s) << " step: " << s << std::endl;
	}

	//std::cout << dM.Derive(0, 0., 0.001) << std::endl;
	std::cout << dM.DeriveOld(0., 0.001) << std::endl;
	//std::cout << dM.Derive(1, 0., 0.001) << std::endl;


	/*double tD = -1.;
	for (int i = 0; i < 30; ++i) {
		std::cout << "\n>> Debug: t = " << tD << std::endl;
		std::cout << tD << " " << dM.ScalarOne(tD) << std::endl;
		std::cout << tD << " " << dM.VectorOne(tD) << std::endl;
		std::cout << tD << " " << dM.ScalarTwo(tD) << std::endl;
		std::cout << tD << " " << dM.VectorTwo(tD) << std::endl;
		std::cout << tD << " " << dM.AbsGEP(tD) << std::endl;
		std::cout << tD << " " << dM.AbsGMP(tD) << std::endl;
		std::cout << tD << " " << dM.AbsGEN(tD) << std::endl;
		std::cout << tD << " " << dM.AbsGMN(tD) << std::endl;
		tD = tD + .05;
	}*/

	/*TComplex a(1.,0.);
	TComplex t0s(t0s,0.);
	TComplex ct(1.05,0.5);
	std::cout << "\n>> Debug: t = " << ct << std::endl;
	std::cout << ct << " " << dM.ScalarOne(ct) << std::endl;*/

	/*ct(1.05,0.);
	std::cout << "\n>> Debug: t = " << ct << std::endl;
	std::cout << ct << " " << dM.ScalarOne(ct) << std::endl;

	ct(1.05,-0.5);
	std::cout << "\n>> Debug: t = " << ct << std::endl;
	std::cout << ct << " " << dM.ScalarOne(ct) << std::endl;*/

	
/*	ofstream myOutput (f);
	const int dPoints = 10000;
	double tMin = -10.;
	double tMax = 10.;
	double tStep = (tMax-tMin)/dPoints;
	TComplex z,gA,gB,gC,gD;

	for (int i = 0; i < dPoints; ++i) {
		z = tMin + i*tStep;
	
		gA = dM.ScalarOne(z);
		gB = dM.VectorOne(z);
		gC = dM.ScalarTwo(z);
		gD = dM.VectorTwo(z);

		myOutput << "\n\t t=" << z.Re() << std::endl;
		myOutput << z.Re() << "\t" << gA.Re() << "\t" << gA.Im() << std::endl;
		myOutput << z.Re() << "\t" << gB.Re() << "\t" << gB.Im() << std::endl;
		myOutput << z.Re() << "\t" << gC.Re() << "\t" << gC.Im() << std::endl;
		myOutput << z.Re() << "\t" << gD.Re() << "\t" << gD.Im() << std::endl;
    }
	myOutput.close();*/
}
