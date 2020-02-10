#ifndef DEF_MATERIAL
#define DEF_MATERIAL

// Header of the class

class Material
{
public :

    // Attributes

    int index ;
    string name ;
    string type ;
    vector<double> parameters ;

    // Constructor and Destructor
    Material(int index, string name, string type, vector<double> parameters) ;
    ~Material() ;

    // Accessors

    // Modifiers

    // Methods

} ;



// Content of the class

// Constructor and Destructor

Material::Material(int i, string n, string t, vector<double> p)
{
    index = i ;
    name = n ;
    type = t ;
    parameters = p ;
}

Material::~Material()
{
}

// Accessors

// Modifiers

// Methods

#endif
