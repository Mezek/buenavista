/**
 * $Date: 2013-05-22 09:56:01 +0200 (Wed, 22 May 2013) $
 * $Revision: 347 $
 * $Author: bartos $
 * $HeadURL: http://triglav/repos/BuenaVista/Hyperons/HyperonsRadii.cpp $
 * $Id: HyperonsRadii.cpp 347 2013-05-22 07:56:01Z bartos $
 *
 * @file
 * @brief	Calculate value of hyperon radii.
 */

using namespace ROOT::Minuit2;

/*struct FF {
	//FF() : name(), nor(), mesons(), tresh() {} // automatic initialization with constructor
	FFactor B(1);
};*/

void performRadii ( char* p, char* f ) {

	std::cout << "\n> Radii parameters:   '" << p << "'" << std::endl;
	std::cout << "> Output:  '" << f << "'" << std::endl;	

	FFactor rBar[6];

	/*FFactor rLam(2);
	FFactor rSpm(3);
	FFactor rSig(4);
	FFactor rXzm(5);*/

	/*// Debug
	FFactor rNic;
	rNic.LoadParameters(p);
	rNic.PrintParameters();

	FFactor rDef;
	rDef.SetParticle(0);
	rDef.LoadParameters(p);
	rDef.PrintParameters();
	
	FFactor rNuc(3);
	*/

	double r = 0.1;
	for ( int i = 1; i < 6; ++i ) {

		rBar[i].SetParticle(i);
		rBar[i].LoadParameters(p);
		//rBar[i].PrintParameters();
		
		std::cout << "GE1N:  " << rBar[i].GE1N() << std::endl;
		std::cout << "GM1N:  " << rBar[i].GM1N() << std::endl;
		std::cout << "GE2N:  " << rBar[i].GE2N() << std::endl;
		std::cout << "GM2N:  " << rBar[i].GM2N() << std::endl;
		
		std::cout << "\n> " << rBar[i].particleName1 << std::endl;
		std::cout << "GE1:  " << rBar[i].RadiusE1(0.00001) << std::endl;
		std::cout << "GM1:  " << rBar[i].RadiusM1(0.01) << std::endl;
		
		std::cout << "\n> " << rBar[i].particleName2 << std::endl;
		std::cout << "GE2:  " << rBar[i].RadiusE2(0.00001) << std::endl;
		std::cout << "GM2:  " << rBar[i].RadiusM2(0.01) << std::endl;

		/*double s = 1.;
		for (int j = 0; j < 6; ++j) {
			s = TMath::Power(0.1,j);
			std::cout << "> Magnetic radius: " << rBar[i].RadiusM1(s) << " " << rBar[i].RadiusM2(s)<< " step: " << s << std::endl;
		}*/

	}

	for (int i = 0; i < 6; ++i) {
	}
}
