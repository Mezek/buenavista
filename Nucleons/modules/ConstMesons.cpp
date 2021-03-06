/**
 * $Date$
 * $Revision$
 * $Author$
 * $HeadURL: http://triglav/repos/BuenaVista/Nucleons/modules/ConstMesons.cpp $
 * $Id$
 *
 * @file
 * @brief	Meson constants for numerical calculations.
 */

	/// Masses and widths are in GeV
	static const double massPi0 = 0.1349766; // Pi0	
	static const double widthPi0 = 0.0006;
	static const double massPiPm = 0.13957018;
	static const double massPi = massPiPm;
	static const double massPi2 = massPi*massPi;

	static const double t0s = 9.*massPi2;                 // isoscalar squared branch point		
	static const double t0v = 4.*massPi2;                 // isovector squared branch point	

	// Om
	static const double massOm = 0.78265;
	static const double widthOm = 0.00849;                // PDG value
	static const double widthEEOm = 0.0000006;            // widthOm*7.28*10-5;

	// Om1P
	static const double massOm1P = 1.425;
	static const double widthOm1P = 0.215;                // PDG range 0.180-0.250
	// Donnachie value A
	//static const double massOm1P = 1.391;                 // +/- 0.018
	//static const double widthOm1P = 0.224;                // +/- 0.049
	static const double widthEEOm1P = 0.137*0.000001;     // +/- 0.057*0.000001 

	// Om2P
	static const double massOm2P = 1.670;
	static const double widthOm2P = 0.315;
	// Donnachie value A
	//static const double massOm2P = 1.594;                 // +/- 0.012
	//static const double widthOm2P = 0.100;                // +/- 0.030
	static const double widthEEOm2P = 0.152*0.000001;     // +/- 0.047*0.000001 

	// Phi
	static const double massPhi = 1.019461;
	static const double widthPhi = 0.004266;              // PDG value
	// Donachie value
	static const double widthEEPhi = 0.000001258;         // widthPhi*0.0002954;

	// Phi1P
	static const double massPhi1P = 1.680;
	//static const double massPhi1P = 1.378;                // D. value
	static const double widthPhi1P = 0.150;

	// Phi2P
	static const double massPhi2P = 2.175;
	static const double widthPhi2P = 0.061;
	//static const double massPhi2P = 1.680;                // D. value
	//static const double widthPhi2P = 0.49;                // D. value

	// Rho
	static const double massRho = 0.77549;
	static const double widthRho = 0.1491;
	// Donachie value
    static const double widthEERho = 0.00000704;          // widthRho*4.72*10-5;

	// Rho1P
	static const double massRho1P = 1.465;
	static const double widthRho1P = 0.400;
	// Donachie value
	static const double widthEERho1P = 9.*0.195*0.000001; // +/- 9.*0.027*0.000001 

	// Rho2P
	static const double massRho2P = 1.720;
	//static const double massRho2P = 1.570;                // D. value  
	static const double widthRho2P = 0.250;
	// Donachie value
	static const double widthEERho2P = 9.*0.084*0.000001; // +/- 9.*0.017*0.000001 

	// Rho3P
	static const double massRho3P = 2.1693;
	static const double widthRho3P = 0.25;

	// Rho4P
	static const double massRho4P = 2.455;                // +/-53 by D. fit
	static const double widthRho4P = 0.728;               // +/-2 by D. fit

