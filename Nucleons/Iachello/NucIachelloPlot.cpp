/**
 * $Date: 2013-10-08 11:56:38 +0200 (Tue, 08 Oct 2013) $
 * $Revision: 390 $
 * $Author: bartos $
 * $HeadURL: http://triglav/repos/BuenaVista/Nucleons/Iachello/NucIachelloPlot.cpp $
 * $Id: NucIachelloPlot.cpp 390 2013-10-08 09:56:38Z bartos $
 *
 * @file
 * @brief	Plot of Iachello form factors.
 */

using namespace ROOT::Minuit2;

/// Perform plot of form factors.

void performPlot ( char* p, char* f1, char* f2 ) {

	std::cout << "\n> Plotting:" << std::endl;
	std::cout << "> Plotted parameters:          `" << p << "'" << std::endl;
	std::cout << "> Plotted form factor data:    `" << f1 << "'" << std::endl;
	std::cout << "> Plotted cross section data:  `" << f2 << "'" << std::endl;

	const int nPoints = 2500;
	double tMin;
	double tMax;
	double tStep;
	double tA;

	ExperimentalData W;
	W.ReadData(f1);
	W.ReadDataCrossSection(f2);
	int numG = W.size();
	std::cout << "\nNumber of plotted points:      " << numG << std::endl;
	std::vector<int> type = W.Type();
	std::vector<double> x = W.X();
	std::vector<double> val = W.Val();
	std::vector<double> errUp = W.ErrUp();
	std::vector<double> errDown = W.ErrDown();
	std::vector<double> theta = W.Theta();
	std::vector<double> energy = W.Energy();
	/*
	FFactor nPlot(12);
	nPlot.LoadParameters(p);
	nPlot.PrintParameters();
	//nPlot.CheckFormFactor("all", 0.01);
	
	double plotX[nPoints], plotY0[nPoints], plotY1[nPoints], plotY2[nPoints], plotY3[nPoints];
	int numA = W.size();
	double dataX0[numA];
	double dataY0[numA];
	double dataU0[numA];
	double dataD0[numA];

	tMin = .1;
	tMax = 1.;
	tStep = (tMax-tMin)/nPoints;
	tA = tMin;
	double res;
	TComplex z0(tMin,0.000001);
	TComplex z,zY,zW,gA,gB,gC,gD;

	for (int i = 0; i < nPoints; ++i) {
		tA = tMin + i*tStep;
		plotX[i] = tA;
 		z = tMin + i*tStep;

		//gA = nPlot.ScalarOne(z);
		//gB = nPlot.VectorOne(z);
		//gC = nPlot.ScalarTwo(z);
		//gD = nPlot.VectorTwo(z);

		TComplex tau = z/(4*massP*massP);
		//TComplex eps = 1./(1.+2.*(1.+tau)*tan(th/2.)*tan(th/2.));
		TComplex a = a.Abs(nPlot.GMp(z));
		TComplex b = b.Abs(nPlot.GEp(z));
		TComplex k = 4./3.*TMath::Pi()*alpha*alpha/z;
		TComplex r = r.Sqrt(1.-4*massP*massP/z);
		TComplex v = k*r*(a*a+2.*massP*massP/z*b*b);	

		//std::cout << plotY0[i] << std::endl;
    }

	for (int i = 0; i < numG; i++) {
		dataX0[i] = W.x[i];
		dataY0[i] = W.val[i];
		dataD0[i] = W.errDown[i];
		dataU0[i] = W.errUp[i];
	}
			
	PlotGraph graf1,graf2;
	Char_t const* title;
	title = "#sigma_{exp}^{p}";

	//graf1.view(numE,dataX0,dataY0,title);
	graf2.viewPlusDataE(nPoints,plotX,plotY0,numG,dataX0,dataY0,0,dataD0,title);
	*/

	// Plot ratios
	FFactor pPlot(5);
	pPlot.LoadParameters(p);

	double plotRX[nPoints], plotRpY[nPoints], plotRnY[nPoints];
	int numE = W.typeN(9);
	int numF = W.typeN(10);
	double dataX4[numE], dataX5[numF];
	double dataY4[numE], dataY5[numF];
	double dataU4[numE], dataU5[numF];
	double dataD4[numE], dataD5[numF];
	double tZ = 0.;
	
	tMin = 0.0;
	tMax = 15.0;
	tStep = (tMax-tMin)/nPoints;
	tA = tMin;
	TComplex w,hA,hB,hC,hD;

	int k4 = 0 , k5 = 0;
	for (int i = 0; i < W.size(); i++) {
		if (type[i] == 9) {
			dataX4[k4] = -x[i];
			dataY4[k4] = val[i];
			dataD4[k4] = errDown[i];
			dataU4[k4] = errUp[i];
			++k4;
		}
		if (type[i] == 10) {
			dataX5[k5] = -x[i];
			dataY5[k5] = val[i];
			dataD5[k5] = errDown[i];
			dataU5[k5] = errUp[i];
			++k5;
		}
		if (type[i] == 0) {
			std::cout << ">> Warning: Some data have undefined type!" << std::endl;
		}
	}

	for (int i = 0; i < nPoints; ++i) {
		tA = tMin + i*tStep;
		plotRX[i] = tA;
 		w = tMin + i*tStep;

		hA = pPlot.GEP(-w);
		hB = pPlot.GMP(-w);
		hC = pPlot.GEN(-w);
		hD = pPlot.GMN(-w);

		plotRpY[i] = (1+ammP)*hA/hB;
		plotRnY[i] = (ammN)*hC/hD;
		//std::cout << i << ". " << plotRpY[i] << " " << tA <<  std::endl;
		if (plotRpY[i] > 0.) { tZ = tA; }
    }
	std::cout << "Zero point at [ GeV^2 ] :      " << tZ << std::endl;

	PlotGraph grafR(2);
	Char_t const* titleR;
	titleR = "mu_p*G_E^p/G_M^p";
	grafR.viewPlusDataAE(nPoints,plotRX,plotRpY,numE,dataX4,dataY4,0,0,dataU4,dataD4,titleR);
	titleR = "mu_n*G_E^n/G_M^n";
	//grafR.viewPlusDataAE(nPoints,plotRX,plotRnY,numF,dataX5,dataY5,0,0,dataU5,dataD5,titleR);

}
