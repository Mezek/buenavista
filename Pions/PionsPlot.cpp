/**
 * $Date$
 * $Revision$
 * $Author$
 *
 * @file
 * @brief	Plot of pion form factors.
 */

using namespace ROOT::Minuit2;

/// Perform plot of form factors.

void performPlot ( char* p, char* f ) {

	std::cout << "\n> Plotted parameters: '" << p << "'" << std::endl;
	std::cout << "> Plotted data:       '" << f << "'" << std::endl;

	FFactor nPlot(11);
	nPlot.LoadParameters(p);
	//nPlot.PrintParameters();
	
	ExperimentalData W;
	W.ReadData(f);
	std::cout << "\nNumber of plotted points:      " << W.size() << std::endl;
	std::vector<double> x = W.X();
	std::vector<double> val = W.Val();
	std::vector<double> errUp = W.ErrUp();
	std::vector<double> errDown = W.ErrDown();
	std::vector<int> type = W.Type();	
	
	const int nPoints = 2500;
	double tMin;
	double tMax;
	double tStep;
	double tA;

	double plotX[nPoints], plotY0[nPoints], plotY1[nPoints], plotY2[nPoints], plotY3[nPoints];
	int num = x.size();
	double dataX0[num], dataX1[num];
	double dataY0[num], dataY1[num];
	double dataU0[num], dataU1[num];
	double dataD0[num], dataD1[num];

	for (int i = 0; i < num; ++i) {
		dataX0[i] = x[i];
		dataY0[i] = val[i];
    }

	tMin = -10.0;
	tMax = 15.0;
	tStep = (tMax-tMin)/nPoints;
	tA = tMin;
	double res;
	TComplex z0(tMin,0.000001);
	TComplex z,zY,zW,gA,gB,gC,gD;

	for (int i = 0; i < nPoints; ++i) {
		tA = tMin + i*tStep;
		plotX[i] = tA;
 		z = tMin + i*tStep;
		plotY0[i] = nPlot.AbsValue(z);
    }
	
	
	PlotGraph graf1(12);
	Char_t const* title;
	title = "|F_{pi}|";
	graf1.viewPlusData(nPoints,plotX,plotY0,num,dataX0,dataY0,title);

}
