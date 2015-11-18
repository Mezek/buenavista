/**
 * $Date$
 * $Revision$
 * $Author$
 *
 * @file
 * @brief	Plot of hyperon form factors.
 */

using namespace ROOT::Minuit2;

void performPlot ( char* p, char* f ) {

	std::cout << "\n> Plotted parameters: '" << p << "'" << std::endl;
	std::cout << "> Plotted data:       '" << f << "'" << std::endl;

	FFactor nPlot(0);
	nPlot.LoadParameters(p);
	//nPlot.PrintParameters();
	nPlot.CheckFormFactor("all", 0.01);

	std::cout << nPlot.ScalarOne(0.) << std::endl;
	std::cout << nPlot.ScalarTwo(0.) << std::endl;
	std::cout << nPlot.VectorOne(0.) << std::endl;
	std::cout << nPlot.VectorTwo(0.) << std::endl;
	std::cout << nPlot.GE1(0.) << std::endl;
	std::cout << nPlot.GM1(0.) << std::endl;
	std::cout << nPlot.GE2(0.) << std::endl;
	std::cout << nPlot.GM2(0.) << std::endl;

	// Plot graphs
	const int nPoints = 2500;
	double plotX[nPoints], plotY0[nPoints], plotY1[nPoints], plotY2[nPoints], plotY3[nPoints];
	double plotY4[nPoints], plotY5[nPoints];
	//double dataX0[numA], dataX1[numB], dataX2[numC], dataX3[numD];
	//double dataY0[numA], dataY1[numB], dataY2[numC], dataY3[numD];
	//double dataU0[numA], dataU1[numB], dataU2[numC], dataU3[numD];
	//double dataD0[numA], dataD1[numB], dataD2[numC], dataD3[numD];
	
	double tMin = -10.;
	double tMax = 15.;
	double tStep = (tMax-tMin)/nPoints;
	double tA = tMin;
	double res;
	TComplex z0(tMin,0.000001);
	TComplex z,zY,zW,gA,gB,gC,gD,hA,hB,hC,hD;

	for (int i = 0; i < nPoints; i++) {
		tA = tMin + i*tStep;
		plotX[i] = tA;
 		z = tMin + i*tStep;

		plotY0[i] = nPlot.AbsGE1(z);
		plotY1[i] = nPlot.AbsGM1(z);
		plotY2[i] = nPlot.AbsGE2(z);
		plotY3[i] = nPlot.AbsGM2(z);
	}


	PlotGraph graf(13);
	PlotGraph graf2(10);
	
	graf.view4(nPoints,plotX,plotY0,plotY1,plotY2,plotY3);

	//graf.view(nPoints,plotX,plotY0,"|G_{E}^{P1}|");
	//graf.view(nPoints,plotX,plotY1,"|G_{M}^{P2}|");
	//graf.view(nPoints,plotX,plotY2,"|G_{E}^{P1}|");
	//graf.view(nPoints,plotX,plotY3,"|G_{M}^{P2}|");

	graf.view(nPoints,plotX,plotY0,"|G_{E}^{#Sigma^{+}}|");
	graf.view(nPoints,plotX,plotY1,"|G_{M}^{#Sigma^{+}}|");
	graf.view(nPoints,plotX,plotY2,"|G_{E}^{#Sigma^{-}}|");
	graf.view(nPoints,plotX,plotY3,"|G_{M}^{#Sigma^{-}}|");

	//graf.view(nPoints,plotX,plotY0,"|G_{E}^{#Xi^{-}}|");
	//graf.view(nPoints,plotX,plotY1,"|G_{M}^{#Xi^{-}}|");
	//graf.view(nPoints,plotX,plotY2,"|G_{E}^{#Xi^{0}}|");
	//graf.view(nPoints,plotX,plotY3,"|G_{M}^{#Xi^{0}}|");

	//graf.view(nPoints,plotX,plotY0,"|G_{E}^{#Sigma^{0}}|");
	//graf.view(nPoints,plotX,plotY1,"|G_{M}^{#Sigma^{0}}|");

	//graf.view(nPoints,plotX,plotY0,"|G_{E}^{#Lambda}|");
	//graf.view(nPoints,plotX,plotY1,"|G_{M}^{#Lambda}|");

/*
	// Plot cross sections
	tMin = 0.;
	tMax = 15.;
	tStep = (tMax-tMin)/nPoints;
	tA = tMin;

	for (int i = 0; i < nPoints; i++) {
		tA = tMin + i*tStep;
		plotX[i] = tA;
 		z = tMin + i*tStep;

		if (z.Re() < 4.*nPlot.massH1*nPlot.massH1) {
			plotY4[i] = 0.;
		} else {
			plotY4[i] = TComplex::Abs(nPlot.SigmaTotal1(z));
		}

		if (z.Re() < 4.*nPlot.massH2*nPlot.massH2) {
			plotY5[i] = 0.;
		} else {
			plotY5[i] = TComplex::Abs(nPlot.SigmaTotal2(z));
		}

	}

	const Char_t* title = "";

    TCanvas *c1 = new TCanvas ("c1", "Graph", 10, 10, 600, 500);
	c1->SetLogy();
    TGraph *gr1 = new TGraph (nPoints, plotX, plotY4);
	
    gr1->Draw("AL");
    gr1->GetXaxis()->CenterTitle();
    gr1->GetYaxis()->CenterTitle();

	gr1->GetXaxis()->SetTitle("t [GeV^{2}]");
	gr1->GetYaxis()->SetTitle("sigma_{tot}^{1}");
	gr1->SetTitle(title);

	TCanvas *c2 = new TCanvas ("c2", "Graph", 10, 10, 600, 500);
	c2->SetLogy();
    TGraph *gr2 = new TGraph (nPoints, plotX, plotY5);
	
    gr2->Draw("AL");
    gr2->GetXaxis()->CenterTitle();
    gr2->GetYaxis()->CenterTitle();

	gr2->GetXaxis()->SetTitle("t [GeV^{2}]");
	gr2->GetYaxis()->SetTitle("sigma_{tot}^{2}");
	gr2->SetTitle(title);
*/
}
