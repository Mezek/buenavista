/**
 * $Date: 2013-05-20 14:30:44 +0200 (Mon, 20 May 2013) $
 * $Revision: 338 $
 * $Author: bartos $
 * $HeadURL: http://triglav/repos/BuenaVista/Nucleons/modules/ConstBasic.cpp $
 * $Id: ConstBasic.cpp 338 2013-05-20 12:30:44Z bartos $
 *
 * @file
 * @brief	Basic constants for numerical calculations.
 */

	static const TComplex kI(0,1);
	static const TComplex k0(0,0);
	static const TComplex k1(1,0);
	
	static const double alpha = 1.0/137.036;              // 1.0/137.035999074;
	static const double hTransC = 197.3269718;            // [MeV*fm]
	static const double hTransC2 = .389379338;            // [GeV^2*mbarn]

	/// Masses and widths are in GeV
	static const double massP = 0.93827203;
	static const double massN = 0.93956536;
	static const double massL = 1.115683;
	static const double massSp = 1.18937;
	static const double massS0 = 1.192642;
	static const double massSm = 1.197449; 
	static const double massX0 = 1.31486;
	static const double massXm = 1.32171;

	static const double massMu = 0.105658369;
	static const double massTau = 1.77699;

	static const double trashP = 4.*massP*massP;
	static const double trashN = 4.*massN*massN;

	static const double muP = 2.792847351;
	static const double muN = -1.9130427;
	static const double muL = -0.613;
	static const double muSp = 2.458;
	static const double muSm = -1.16;
	static const double muX0 = -1.25;
	static const double muXm = -0.6507;

	static const double ammP = 1.792847351;
	static const double ammN = -1.9130427;

