/**
 * $Date$
 * $Revision$
 * $Author$
 * $HeadURL: http://triglav/repos/BuenaVista/Nucleons/modules/ConstMesons-D.cpp $
 * $Id$
 *
 * @file
 * @brief	Meson constants, multiple of mass of &pi; meson, v. D.
 */

	/// Masses and widths are in GeV
	static const double massPi = 0.14; // 0.13957018;
	static const double massPi0 = massPi;
	static const double massPi2 = massPi*massPi;

	static const double t0s = 9.*massPi2;                 // isoscalar squared branch point		
	static const double t0v = 4.*massPi2;                 // isovector squared branch point	

	// Om
	static const double massOm = massPi*5.590357;         // 
	static const double widthOm = massPi*0.0642857;       // 
	static const double widthEEOm = 0.0000006;            // 

	// Om1P
	static const double massOm1P = massPi*10.178571;      // 
	static const double widthOm1P = massPi*1.535714;      // 
	// Donnachie value A
	//static const double massOm1P = 1.391;               // +/- 0.018
	//static const double widthOm1P = 0.224;              // +/- 0.049
	static const double widthEEOm1P = 0.137*0.000001;     // +/- 0.057*0.000001 

	// Om2P
	static const double massOm2P = massPi*11.92857;       // 
	static const double widthOm2P = massPi*2.25;          // 
	// Donnachie value A
	//static const double massOm2P = 1.594;               // +/- 0.012
	//static const double widthOm2P = 0.100;              // +/- 0.030
	static const double widthEEOm2P = 0.152*0.000001;     // +/- 0.047*0.000001 

	// Phi
	static const double massPhi = massPi*7.281821;        // 
	static const double widthPhi = massPi*0.0304286;      // 
	// Donachie value
	static const double widthEEPhi = 0.000001258;         // widthPhi*0.0002954;

	// Phi1P
	static const double massPhi1P = massPi*12.;           // 
	static const double widthPhi1P = massPi*1.07143;      // 

	// Phi2P
	static const double massPhi2P = massPi*15.535714;     // 
	static const double widthPhi2P = massPi*0.435714;     // 

	// Rho
    static const double massRho = massPi*5.539214;        // 
	static const double widthRho = massPi*1.065;          // 
	// Donachie value
    static const double widthEERho = 0.00000704;          // widthRho*4.72*10-5;

	// Rho1P
	static const double massRho1P = massPi*10.464285;     // 
	static const double widthRho1P = massPi*2.857143;     // 
	// Donachie value
	static const double widthEERho1P = 9.*0.195*0.000001; // +/- 9.*0.027*0.000001 

	// Rho2P
	static const double massRho2P = massPi*12.285714;     // 
	static const double widthRho2P = massPi*1.785714;     // 
	// Donachie value
	static const double widthEERho2P = 9.*0.084*0.000001; // +/- 9.*0.017*0.000001 

	// Rho3P
	static const double massRho3P = massPi*15.495;        // 
	static const double widthRho3P = massPi*1.785714;     // 

	// Rho4P
	static const double massRho4P = massPi*17.535714;     // D. fit
	static const double widthRho4P = massPi*5.2;          // D. fit

