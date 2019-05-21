#ifndef DEF_GAUSS
#define DEF_GAUSS

// Header of the class

class Gauss
{
	public :
	
	// Static attributes
	double x_gauss ;
	double y_gauss ;
	double weight ;
	double jacobian ;
	int    number_influencing_nodes ;
	vector<int> influencing_nodes ;
	vector<double> shape_functions ;
	vector<double> shape_xderiv ;
	vector<double> shape_yderiv ;
	
	// Constructor and Destructor
	Gauss(double x, double y, double w, double j, int n) ;
	~Gauss() ;
	
	// Accessors
	
	// Modifiers
	
	// Methods
} ;




// Content of the class

// Constructor and Destructor

Gauss::Gauss (double x, double y, double w, double j, int n)
{
	x_gauss = x ;
	y_gauss = y ;
	weight = w ;
	jacobian = j ;
	number_influencing_nodes = n ;
	//vector<int> influencing_nodes ;
	//vector<double> shape_functions ;
	//vector<double> shape_xderiv ;
	//vector<double> shape_yderiv ;
}

Gauss::~Gauss()
{
}


// Accessors

// Modifiers

// Methods


#endif
