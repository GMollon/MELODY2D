#ifndef DEF_CONTACT_LAW
#define DEF_CONTACT_LAW

// Header of the class

class Contact_law
{
public :

    // Attributes
    int index ;
    string material1 ;
    string material2 ;
    string type ;
    string length_evolution ;
    vector<double> parameters ;

    // Constructor and Destructor
    Contact_law(int index, string material1, string material2, string type, string length_evolution, vector<double> parameters) ;
    ~Contact_law() ;

    // Accessors

    // Modifiers

    // Methods

} ;



// Content of the class

// Constructor and Destructor

Contact_law::Contact_law(int i, string m1, string m2, string t, string l, vector<double> p)
{
    index = i ;
    material1 = m1 ;
    material2 = m2 ;
    type = t ;
    length_evolution = l ;
    parameters = p ;
}

Contact_law::~Contact_law()
{
}

// Accessors

// Modifiers

// Methods


#endif
