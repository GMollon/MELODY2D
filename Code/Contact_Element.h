#ifndef DEF_CONTACT_ELEM
#define DEF_CONTACT_ELEM

// Header of the class

class Contact_element
{
public :

    // Attributes
    int    borderS ;
    int    border_nodeS ;
    int    nodeS ;
    int    bodyM ;
    int    borderM ;
    int    shiftM ;
    int    segmentM ;
    int    nodeM0 ;
    int    shiftM0 ;
    double shapeM0 ;
    int    nodeM1 ;
    int    shiftM1 ;
    double shapeM1 ;
    int    nodeM2 ;
    int    shiftM2 ;
    double shapeM2 ;
    int    nodeM3 ;
    int    shiftM3 ;
    double shapeM3 ;
    double gapn ;
    double vgapn ;
    double gapt ;
    double vgapt ;
    double xsi ;
    double xnorm ;
    double ynorm ;
    double xtan ;
    double ytan ;
    double effective_mass ;
    double length ;
    double fx ;
    double fy ;
    int nb_internal ;
    vector<double> internal ;

    // Constructor and Destructor
    Contact_element() ;
    ~Contact_element() ;

    // Accessors

    // Modifiers

    // Methods
} ;




// Content of the class

// Constructor and Destructor
Contact_element::Contact_element()
{
    borderS = 0 ;
    border_nodeS = 0 ;
    nodeS = 0 ;
    bodyM = 0 ;
    borderM = 0 ;
    shiftM = 0 ;
    segmentM = 0 ;
    nodeM0 = 0 ;
    shiftM0 = 0 ;
    shapeM0 = 0. ;
    nodeM1 = 0 ;
    shiftM1 = 0 ;
    shapeM1 = 0. ;
    nodeM2 = 0 ;
    shiftM2 = 0 ;
    shapeM2 = 0. ;
    nodeM3 = 0 ;
    shiftM3 = 0 ;
    shapeM3 = 0. ;
    gapn = 0. ;
    vgapn = 0. ;
    gapt = 0. ;
    vgapt = 0. ;
    xsi = 0. ;
    xnorm = 0. ;
    ynorm = 0. ;
    xtan = 0. ;
    ytan = 0. ;
    effective_mass = 0. ;
    length = 0. ;
    fx = 0. ;
    fy = 0. ;
    nb_internal = 0 ;
}

Contact_element::~Contact_element()
{
}

// Accessors

// Modifiers

// Methods


#endif
