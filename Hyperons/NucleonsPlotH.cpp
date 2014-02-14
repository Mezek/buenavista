/**
 * $Date$
 * $Revision$
 * $Author$
 *
 * @file
 * @brief	Plot of nucleon form factors, v. hyperon.
 */

using namespace ROOT::Minuit2;

/// Perform plot of form factors.

void performPlot ( char* p, char* f ) {

	std::cout << "\n> Plotting:" << std::endl;
	std::cout << "> Plotted parameters: `" << p << "'" << std::endl;
	std::cout << "> Plotted data:       `" << f << "'" << std::endl;

	const int nPoints = 2500;
	double tMin;
	double tMax;
	double tStep;
	double tA;

	ExperimentalData W;
	W.ReadData(f);
	int numG = W.size();
	std::cout << "\nNumber of plotted points:      " << numG << std::endl;
	std::vector<int> type = W.Type();
	std::vector<double> x = W.X();
	std::vector<double> val = W.Val();
	std::vector<double> errUp = W.ErrUp();
	std::vector<double> errDown = W.ErrDown();
	std::vector<double> theta = W.Theta();
	std::vector<double> energy = W.Energy();

	FFactor nPlot;
	nPlot.LoadParameters(p);
	//nPlot.PrintParameters();
	nPlot.CheckFormFactor("all", 0.01);
	
	double plotX[nPoints], plotY0[nPoints], plotY1[nPoints], plotY2[nPoints], plotY3[nPoints];
	int numA = W.typeN(1)+W.typeN(2);
	int numB = W.typeN(3)+W.typeN(4);
	int numC = W.typeN(5)+W.typeN(6);
	int numD = W.typeN(7)+W.typeN(8);
	double dataX0[numA], dataX1[numB], dataX2[numC], dataX3[numD];
	double dataY0[numA], dataY1[numB], dataY2[numC], dataY3[numD];
	double dataU0[numA], dataU1[numB], dataU2[numC], dataU3[numD];
	double dataD0[numA], dataD1[numB], dataD2[numC], dataD3[numD];

	//std::fill(dataX3, dataX3+num3, 0);
	tMin = -25.0;
	tMax = 25.0;
	tStep = (tMax-tMin)/nPoints;
	tA = tMin;
	double res;
	TComplex z0(tMin,0.000001);
	TComplex z,zY,zW,gA,gB,gC,gD;

	for (int i = 0; i < nPoints; ++i) {
		tA = tMin + i*tStep;
		plotX[i] = tA;
		//~???
 		//~if (tA < t0v) { 
 		z = tMin + i*tStep;
 		//~}
 		//~else {
 			//~z = z0 + i*shag;
 		//~}

		/*
		gA = nPlot.ScalarOne(z);
		gB = nPlot.VectorOne(z);
		gC = nPlot.ScalarTwo(z);
		gD = nPlot.VectorTwo(z);

		plotY0[i] = TComplex::Abs(gA+gB + z/(4*massP*massP)*(gC+gD));
		plotY1[i] = TComplex::Abs(gA+gB + gC+gD);
		plotY2[i] = TComplex::Abs(gA-gB + z/(4*massN*massN)*(gC-gD));
		plotY3[i] = TComplex::Abs(gA-gB + gC-gD);*/

		plotY0[i] = nPlot.AbsGE1(z);
		plotY1[i] = nPlot.AbsGM1(z);
		plotY2[i] = nPlot.AbsGE2(z);
		plotY3[i] = nPlot.AbsGM2(z);

    }
	
	int k0 = 0 , k1 = 0, k2 = 0, k3 = 0;
	for (int i = 0; i < W.size(); i++) {
		//std::cout << type[i] << ": " << x[i] << " " << val[i] << std::endl;		
		// "protonElectric"
		if ((type[i] == 1) || (type[i] == 2)) {
			dataX0[k0] = x[i];
			dataY0[k0] = val[i];
			dataD0[k0] = errDown[i];
			dataU0[k0] = errUp[i];
			//std::cout << i << ": " << x[i] << " " << val[i] << std::endl;
			++k0;
		}
		// "protonMagnetic"
		if ((type[i] == 3) || (type[i] == 4)) {
			dataX1[k1] = x[i];
			dataY1[k1] = val[i];
			dataD1[k1] = errDown[i];
			dataU1[k1] = errUp[i];
			//if (type[i] == 4) { dataY1[k1] = 0.; }
			//std::cout << " " << k1 << ". " << dataX1[k1] << " " << dataY1[k1] << std::endl;
			++k1;
		}
		// "neutronElectric"
		if ((type[i] == 5) || (type[i] == 6)) {
			dataX2[k2] = x[i];
			dataY2[k2] = val[i];
			dataD2[k1] = errDown[i];
			dataU2[k1] = errUp[i];
			//std::cout << W.type[i] << " " << k2 << ". " << dataX2[k2] << " " << dataY2[k2] << std::endl;
			++k2;
		}
		// "neutronMagnetic"
		if ((type[i] == 7) || (type[i] == 8)) {
			dataX3[k3] = x[i];
			dataD3[k3] = errDown[i];
			dataU3[k3] = errUp[i];
			if (type[i] == 7) { dataY3[k3] = -val[i]; }
				else { dataY3[k3] = val[i]; }
			//std::cout << " " << k3 << ". " << dataX3[k3] << " " << dataY3[k3] << std::endl;
			++k3;
		}
		if (type[i] == 9) {
		}
		if (type[i] == 10) {
		}
		if (type[i] == 0) {
			std::cout << ">> Warning: Some data have undefined type!" << std::endl;
		}
	}

	// Small check of values
	if ( k0 != numA ) {	std::cout << "> Warning: data GEp: " << k0 << " " << numA << std::endl; }
	if ( k1 != numB ) {	std::cout << "> Warning: data GMp: " << k1 << " " << numB << std::endl; }
	if ( k2 != numC ) {	std::cout << "> Warning: data GEn: " << k2 << " " << numC << std::endl; }
	if ( k3 != numD ) {	std::cout << "> Warning: data GMn: " << k3 << " " << numD << std::endl; }
	
	PlotGraph graf(13);
	Char_t const* title;
	title = "|G_{E}^{p}|";
	//graf.viewPlusData(nPoints,plotX,plotY0,k0,dataX0,dataY0,title);
	graf.viewPlusDataAE(nPoints,plotX,plotY0,k0,dataX0,dataY0,0,0,dataD0,dataU0,title);
	title = "|G_{M}^{p}|";
	//graf.viewPlusData(nPoints,plotX,plotY1,k1,dataX1,dataY1,title);
	graf.viewPlusDataAE(nPoints,plotX,plotY1,k1,dataX1,dataY1,0,0,dataD1,dataU1,title);
	title = "|G_{E}^{n}|";
	//graf.viewPlusData(nPoints,plotX,plotY2,k2,dataX2,dataY2,title);
	graf.viewPlusDataAE(nPoints,plotX,plotY2,k2,dataX2,dataY2,0,0,dataD2,dataU2,title);
	title = "|G_{M}^{n}|";
	//graf.viewPlusData(nPoints,plotX,plotY3,k3,dataX3,dataY3,title);
	graf.viewPlusDataAE(nPoints,plotX,plotY3,k3,dataX3,dataY3,0,0,dataD3,dataU3,title);
	
	//graf.view4(nPoints,plotX,plotY0,plotY1,plotY2,plotY3);
	//graf.view4Exp(nPoints,plotX,plotY0,plotY1,plotY2,plotY3,k0,dataX0,dataY0,k1,dataX1,dataY1,k2,dataX2,dataY2,k3,dataX3,dataY3);

	
	// Plot ratios
	FFactor pPlot;
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
	tMax = 19.0;
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

		hA = pPlot.GE1(-w);
		hB = pPlot.GM1(-w);
		hC = pPlot.GE2(-w);
		hD = pPlot.GM2(-w);

		plotRpY[i] = muP*hA/hB;
		plotRnY[i] = muN*hC/hD;
		//std::cout << i << ". " << w << " " << hA <<  std::endl;
		if (plotRpY[i] > 0.) { tZ = tA; }
    }
	std::cout << "Zero point at [ GeV^2 ] :      " << tZ << std::endl;

	PlotGraph grafR(2);
	Char_t const* titleR;
	titleR = "#mu_{p}G_{E}^{p}/G_{M}^{p}";
	graf.viewPlusDataAE(nPoints,plotRX,plotRpY,numE,dataX4,dataY4,0,0,dataU4,dataD4,titleR);
	titleR = "#mu_{p}G_{E}^{n}/G_{M}^{n}";
	//graf.viewPlusDataAE(nPoints,plotRX,plotRnY,numF,dataX5,dataY5,0,0,dataU5,dataD5,titleR);

	/*// Plot isoFF
	FFactor oPlot(12);
	oPlot.LoadParameters(outputFile);
	oPlot.PerformCheck("all", 0.001);
	
	double plotIX[nPoints], plotIY[nPoints];
	tMin = -0.5;
	tMax = 0.0;
	tStep = (tMax-tMin)/nPoints;
	tA = tMin;
	TComplex hA;
	
	for (int i = 0; i < nPoints; ++i) {
		tA = tMin + i*tStep;
		plotIX[i] = tA;
 		z = tMin + i*tStep;

		//hA = oPlot.ScalarOne(z);
		//hA = nPlot.VectorOne(z);
		hA = nPlot.ScalarTwo(z);
		//hA = nPlot.VectorTwo(z);

		plotIY[i] = TComplex::Abs(hA);
		//std::cout << i << ". " << z << " " << hA <<  std::endl;
    }

	PlotGraph grafI(1);
	const Char_t* tit = "IsoFF";
	grafI.view(nPoints,plotIX,plotIY,tit);
	*/

}
