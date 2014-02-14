/**
 * $Date$
 * $Revision$
 * $Author$
 *
 * @file
 * @brief	Debug file, v. hyperon.
 */

using namespace ROOT::Minuit2;

void performDebug ( char* p, char* f ) {

	std::cout << "\n> Debug parameters:   '" << p << "'" << std::endl;
	std::cout << "> Debug data output:  '" << f << "'" << std::endl;	
	std::cout << "\nDebug stuff:" << std::endl;

	/*
	// Test argument parsing

	std::cout << "Arguments:" << std::endl;
	std::cout << "-f: " << globalArgs.fit << std::endl;
	std::cout << "-p: " << globalArgs.parameters << std::endl;
	std::cout << "-o: " << globalArgs.output << std::endl;
	std::cout << "-d: " << globalArgs.data << std::endl;
	std::cout << "-v: " << globalArgs.verbose << std::endl;*/

	/*// Convert parameters
	double q1 = 3.365414;
	double q2 = 3.375674;
	double q3 = 5.971356;
	double q4 = 5.384470;
	std::cout << ">> 1s : "<< (q1*q1+1.)*t0s << std::endl;
	std::cout << ">> 1v : "<< (q3*q3+1.)*t0v << std::endl;
	std::cout << ">> 2s : "<< (q2*q2+1.)*t0s << std::endl;
	std::cout << ">> 2v : "<< (q4*q4+1.)*t0v << std::endl;*/
	

	/*// Inverse Convert parameters
	double t1 = 2.01907;
	double t2 = 2.52493;
	double t3 = 1.83342;
	double t4 = 2.31511;
	std::cout << ">> 1s : "<< sqrt(t1/t0s-1.) << std::endl;
	std::cout << ">> 1v : "<< sqrt(t2/t0v-1.) << std::endl;
	std::cout << ">> 2s : "<< sqrt(t3/t0s-1.) << std::endl;
	std::cout << ">> 2v : "<< sqrt(t4/t0v-1.) << std::endl;*/

	FFactor dM;
    dM.LoadParameters(p);
	double r = 0.1;
	for (int i = 0; i < 6; ++i) {
		double s = 1.;
		for (int j = 0; j < i+1; ++j) {
			s = s*r;
		}
		std::cout << "> Proton radius: " << dM.RadiusE1(s) << " step: " << s << std::endl;
	}
	
	/*
	double tD = -0.5;
	for (int i = 0; i < 50; ++i) {
		std::cout << "\n>> Debug:" << std::endl;
		std::cout << tD << " " << dM.ScalarOne(tD) << std::endl;
		std::cout << tD << " " << dM.VectorOne(tD) << std::endl;
		std::cout << tD << " " << dM.ScalarTwo(tD) << std::endl;
		std::cout << tD << " " << dM.VectorTwo(tD) << std::endl;
		std::cout << tD << " " << dM.AbsGEp(tD) << std::endl;
		std::cout << tD << " " << dM.AbsGMp(tD) << std::endl;
		std::cout << tD << " " << dM.AbsGEn(tD) << std::endl;
		std::cout << tD << " " << dM.AbsGMn(tD) << std::endl;
		tD = tD + .05;
	}*/

	/*
	ofstream myOutput (f);
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
	myOutput.close();
	*/
}
