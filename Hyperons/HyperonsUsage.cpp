/**
 * $Date$
 * $Revision$
 * $Author$
 *
 * @file
 * @brief	Hyperon form factors usage.
 */

int deb;
struct globalArgs_t {
	char *parameters;
	char *output;
	char *data;
	int verbose;
} globalArgs;	

static const char *optString = "p:o:d:h?";

static struct option longOpts[] =
{
	{"parameters", required_argument, NULL, 'p'},
	{"output", required_argument, NULL, 'o'},
	{"data", required_argument, NULL, 'd'},
	{"help", no_argument, NULL, 'h'},
	{"verbose", no_argument, &deb, 1},
	{NULL, no_argument, NULL, 0}
};

// Display program usage, and exit.
void displayUsage( void )
{
	std::cout << "\nUsage: ./Hyperons [options]" << std::endl;
	std::cout << "\nOptions:" << std::endl;	
	std::cout << " -p, --parameters       parametersFile, default='" << dataFile << "'" << std::endl;
	std::cout << " -o, --output           outputFile,     default='" << outputFile << "'" << std::endl;
	std::cout << " -d, --data             dataFile,       default='" << parametersFile << "'" << std::endl;
	std::cout << " -h, --help             display usage" << std::endl;
	std::cout << "     --verbose          print some debug stuff\n" << std::endl;

	exit( EXIT_FAILURE );
}

void performUsage ( int argc, char **argv ) {

	// Argument parsing

	int opt = 0;
	int longIndex = 0;
	
	// Initialize globalArgs before we get to work
	globalArgs.parameters = parametersFile;
	globalArgs.output = outputFile;
	globalArgs.data = dataFile;

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
				globalArgs.data = optarg;
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
	if ( deb==1 ) { globalArgs.verbose = 1; }
	
}
