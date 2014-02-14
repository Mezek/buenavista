/**
 * $Date$
 * $Revision$
 * $Author$
 * 
 * @file
 * @brief       Nucleon form factors usage, v. with cross section data.
 */

int val;
int tab;
int deb;
struct globalArgs_t {
	int fit;
	int tabulate;
	char *parameters;
	char *output;
	char *ffdata;
	char *csdata;
	int verbose;
} globalArgs;	

static const char *optString = "p:o:d:c:h?";

static struct option longOpts[] =
{
	{"fit", no_argument, &val, 1},
	{"tabulate", no_argument, &tab, 1},
	{"parameters", required_argument, NULL, 'p'},
	{"output", required_argument, NULL, 'o'},
	{"ffdata", required_argument, NULL, 'd'},
	{"csdata", required_argument, NULL, 'c'},
	{"help", no_argument, NULL, 'h'},
	{"verbose", no_argument, &deb, 1},
	{NULL, no_argument, NULL, 0}
};

/// Display program usage and exit.

void displayUsage( void )
{
	std::cout << "\nUsage: ./NucleonsM [options]" << std::endl;
	std::cout << "\nOptions:" << std::endl;	
	std::cout << "     --fit          set fit,              default= no fit" << std::endl;
	std::cout << "     --tabulate     numerical tabulation, default= no tabulation" << std::endl;
	std::cout << " -p, --parameters   parametersFile,       default='" << parametersFile << "'" << std::endl;
	std::cout << " -o, --output       outputFile,           default='" << outputFile << "'" << std::endl;
	std::cout << " -d, --ffdata       form factor data,     default='" << dataFile1 << "'" << std::endl;
	std::cout << " -c, --csdata       cross section data,   default='" << dataFile2 << "'" << std::endl;
	std::cout << " -h, --help         display usage" << std::endl;
	std::cout << "     --verbose      print some debug stuff\n" << std::endl;

	exit( EXIT_FAILURE );
}

/// Process global and optional arguments.

void performUsage ( int argc, char **argv ) {

	// Argument parsing

	int opt = 0;
	int longIndex = 0;
	
	// Initialize globalArgs before we get to work
	globalArgs.fit = 0;
	globalArgs.tabulate = 0;
	globalArgs.parameters = parametersFile;
	globalArgs.output = outputFile;
	globalArgs.ffdata = dataFile1;
	globalArgs.csdata = dataFile2;

	// Process the arguments with getopt_long()
	
	opt = getopt_long( argc, argv, optString, longOpts, &longIndex );

	while( opt != -1 ) {

		switch( opt ) {
			case 'p':
				globalArgs.parameters = optarg;
				break;
				
			case 'o':
				globalArgs.output = optarg;
				break;
				
			case 'd':
				globalArgs.ffdata = optarg;
				break;

			case 'c':
				globalArgs.csdata = optarg;
				break;
				
			case 'h':
				displayUsage();
				break;
				
			case '?':
				displayUsage();
				break;

			default:
				//displayUsage();
				break;
		}
		
		opt = getopt_long( argc, argv, optString, longOpts, &longIndex );
	}
	if ( val==1 ) { globalArgs.fit = 1; }
	if ( tab==1 ) { globalArgs.tabulate = 1; }
	if ( deb==1 ) { globalArgs.verbose = 1; }
	
}
