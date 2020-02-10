#ifndef DEF_SPY
#define DEF_SPY

// Header of the class

class Spy
{
public :

    // Attributes
    string filename ;
    int    nb_quantities ;
    double period ;
    vector<vector<double>> quantities ;
    double next_time ;

    // Constructor and Destructor
    Spy(string filename, int nb_quantities, double period, vector<vector<double>> quantities) ;
    ~Spy() ;

    // Accessors

    // Modifiers

    // Methods
} ;




// Content of the class

// Constructor and Destructor
Spy::Spy(string f, int nq, double p, vector<vector<double>> q)
{
    filename = f ;
    nb_quantities = nq ;
    period = p ;
    quantities = q ;
    next_time = 0. ; //p ;
}

Spy::~Spy()
{
}

// Accessors

// Modifiers

// Methods

#endif
