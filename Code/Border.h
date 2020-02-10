#ifndef DEF_BORDER
#define DEF_BORDER

// Header of the class

class Border
{
public :

    // Static attributes
    int            number_border_nodes ;
    vector<int>    border_nodes ;
    string         periodicity ;
    string         interpolant ;
    int			   body ;
    int			   index ;
    int            number_xdirichlet_values ;
    vector<double> xdirichlet_instants ;
    vector<double> xdirichlet_values ;
    vector<double> xdirichlet_parameters ;
    int            number_ydirichlet_values ;
    vector<double> ydirichlet_instants ;
    vector<double> ydirichlet_values ;
    vector<double> ydirichlet_parameters ;
    int            number_xneumann_values ;
    vector<double> xneumann_instants ;
    vector<double> xneumann_values ;
    int            number_yneumann_values ;
    vector<double> yneumann_instants ;
    vector<double> yneumann_values ;

    // Dynamic attributes
    string         x_bc_type ;
    double         x_bc_value ;
    double         x_bc_velocity ;
    string         y_bc_type ;
    double         y_bc_value ;
    double         y_bc_velocity ;
    vector<double> length ;
    vector<double> initial_length ;
    vector<int>    node0 ;
    vector<int>    shift0 ;
    vector<int>    node1 ;
    vector<int>    shift1 ;
    vector<int>    node2 ;
    vector<int>    shift2 ;
    vector<int>    node3 ;
    vector<int>    shift3 ;
    vector<double> x_contact_pressure ;
    vector<double> y_contact_pressure ;
    vector<double> x_bc_pressure ;
    vector<double> y_bc_pressure ;
    double		   xmin_prox ;
    double		   xmax_prox ;
    double		   ymin_prox ;
    double		   ymax_prox ;
    int			   shiftmin ;
    int			   shiftmax ;
    int            border_before ;
    int            border_after ;

    // Constructor and Destructor
    Border(int n, string p, int b, int i) ;
    ~Border() ;

    // Accessors

    // Modifiers

    // Methods
    void Update_bc(double t) ;
} ;




// Content of the class

// Constructor and Destructor
Border::Border (int n, string p, int b, int i)
{
    number_border_nodes = n ;
    periodicity = p ;
    body = b ;
    index = i ;
    interpolant = "Linear" ;//"Spline" ;//
    xmin_prox = 0. ;
    xmax_prox = 0. ;
    ymin_prox = 0. ;
    ymax_prox = 0. ;
    shiftmin = 0 ;
    shiftmax = 0 ;
    x_contact_pressure = {0.} ;
    y_contact_pressure = {0.} ;
    for (int i(0) ; i<n ; i++)
    {
        x_contact_pressure.push_back(0.) ;
        y_contact_pressure.push_back(0.) ;
    }
    border_before = 0 ;
    border_after = 0 ;
}

Border::~Border()
{
}

// Accessors

// Modifiers

// Methods



//********************************************//
//** UPDATE BC *******************************//
//********************************************//

void Border::Update_bc(double Time)
{
    double t0, t1, p0, p1 ;

    // Setting current bc along x (type and value)


    //
    if (x_bc_type != "driven")
    {
        if (x_bc_type == "soft")
        {
            int flag_exit = 0 ;
            int n = 0 ;
            x_bc_value = 0. ;
            while ( flag_exit == 0 )
            {
                t0 = xdirichlet_instants[n] ;
                p0 = xdirichlet_values[n] ;
                t1 = xdirichlet_instants[n+1] ;
                p1 = xdirichlet_values[n+1] ;
                if ( Time >= t1 )
                {
                    x_bc_value = x_bc_value + ( t1 - t0 ) * ( p1 + p0 ) * 0.5 ;
                }
                else
                {
                    x_bc_value = x_bc_value + ( Time - t0 ) * ( p0 + ( p1 - p0 ) / ( t1 - t0 ) * ( ( Time - t0 ) * 0.5 ) ) ;
                    x_bc_velocity = p0 + ( p1 - p0 ) / ( t1 - t0 ) * ( Time - t0 ) ;
                    flag_exit = 1 ;
                }
                n++ ;
            }
        }
        else if (x_bc_type == "oriented" || x_bc_type == "following")
        {
            for (int n(0) ; n<number_xneumann_values ; n++)
            {
                t0 = xneumann_instants[n] ;
                p0 = xneumann_values[n] ;
                t1 = xneumann_instants[n+1] ;
                p1 = xneumann_values[n+1] ;
                if ( (t0 <= Time) & (Time <= t1) )
                {
                    x_bc_value = p0 + ( p1 - p0 ) / ( t1 - t0 ) * ( Time - t0 ) ;
                    break ;
                }
            }
        }
        else if (x_bc_type == "none")
        {
            x_bc_value = 0. ;
        }
    }

    if (y_bc_type != "driven")
    {
        if (y_bc_type == "soft")
        {
            int flag_exit = 0 ;
            int n = 0 ;
            y_bc_value = 0. ;
            while ( flag_exit == 0 )
            {
                t0 = ydirichlet_instants[n] ;
                p0 = ydirichlet_values[n] ;
                t1 = ydirichlet_instants[n+1] ;
                p1 = ydirichlet_values[n+1] ;
                if ( Time >= t1 )
                {
                    y_bc_value = y_bc_value + ( t1 - t0 ) * ( p1 + p0 ) * 0.5 ;
                }
                else
                {
                    y_bc_value = y_bc_value + ( Time - t0 ) * ( p0 + ( p1 - p0 ) / ( t1 - t0 ) * ( ( Time - t0 ) * 0.5 ) ) ;
                    y_bc_velocity = p0 + ( p1 - p0 ) / ( t1 - t0 ) * ( Time - t0 ) ;
                    flag_exit = 1 ;
                }
                n++ ;
            }
        }
        else if (y_bc_type == "oriented" || y_bc_type == "following")
        {
            for (int n(0) ; n<number_yneumann_values ; n++)
            {
                t0 = yneumann_instants[n] ;
                p0 = yneumann_values[n] ;
                t1 = yneumann_instants[n+1] ;
                p1 = yneumann_values[n+1] ;
                if ( (t0 <= Time) & (Time <= t1) )
                {
                    y_bc_value = p0 + ( p1 - p0 ) / ( t1 - t0 ) * ( Time - t0 ) ;
                    break ;
                }
            }
        }
        else if (y_bc_type == "none")
        {
            y_bc_value = 0. ;
        }
    }
    /*
    if (x_bc_type != "driven")
    {
        if ( number_xneumann_values != 0 )
        {
            x_bc_type = "neumann" ;
            for (int n(0) ; n<number_xneumann_values ; n++)
            {
                t0 = xneumann_instants[n] ;
                p0 = xneumann_values[n] ;
                t1 = xneumann_instants[n+1] ;
                p1 = xneumann_values[n+1] ;
                if ( (t0 <= Time) & (Time <= t1) )
                {
                    x_bc_value = p0 + ( p1 - p0 ) / ( t1 - t0 ) * ( Time - t0 ) ;
                    break ;
                }
            }
        }
        else if ( number_xdirichlet_values != 0 )
        {
            x_bc_type = "dirichlet" ;
            int flag_exit = 0 ;
            int n = 0 ;
            x_bc_value = 0. ;
            while ( flag_exit == 0 )
            {
                t0 = xdirichlet_instants[n] ;
                p0 = xdirichlet_values[n] ;
                t1 = xdirichlet_instants[n+1] ;
                p1 = xdirichlet_values[n+1] ;
                if ( Time >= t1 )
                {
                    x_bc_value = x_bc_value + ( t1 - t0 ) * ( p1 + p0 ) * 0.5 ;
                }
                else
                {
                    x_bc_value = x_bc_value + ( Time - t0 ) * ( p0 + ( p1 - p0 ) / ( t1 - t0 ) * ( ( Time - t0 ) * 0.5 ) ) ;
                    flag_exit = 1 ;
                }
                n++ ;
            }
        }
        else
        {
            x_bc_type = "none" ;
            x_bc_value = 0. ;
        }
    }

    // Setting current bc along y (type and value)
    if (y_bc_type != "driven")
    {
        if ( number_yneumann_values != 0 )
        {
            y_bc_type = "neumann" ;
            for (int n(0) ; n<number_yneumann_values ; n++)
            {
                t0 = yneumann_instants[n] ;
                p0 = yneumann_values[n] ;
                t1 = yneumann_instants[n+1] ;
                p1 = yneumann_values[n+1] ;
                if ( (t0 <= Time) & (Time <= t1) )
                {
                    y_bc_value = p0 + ( p1 - p0 ) / ( t1 - t0 ) * ( Time - t0 ) ;
                    break ;
                }
            }
        }
        else if ( number_ydirichlet_values != 0 )
        {
            if (y_bc_type != "driven") y_bc_type = "dirichlet" ;
            int flag_exit = 0 ;
            int n = 0 ;
            y_bc_value = 0. ;
            while ( flag_exit == 0 )
            {
                t0 = ydirichlet_instants[n] ;
                p0 = ydirichlet_values[n] ;
                t1 = ydirichlet_instants[n+1] ;
                p1 = ydirichlet_values[n+1] ;
                if ( Time >= t1 )
                {
                    y_bc_value = y_bc_value + ( t1 - t0 ) * ( p1 + p0 ) * 0.5 ;
                }
                else
                {
                    y_bc_value = y_bc_value + ( Time - t0 ) * ( p0 + ( p1 - p0 ) / ( t1 - t0 ) * ( ( Time - t0 ) * 0.5 ) ) ;
                    flag_exit = 1 ;
                }
                n++ ;
            }
        }
        else
        {
            y_bc_type = "none" ;
            y_bc_value = 0. ;
        }
    }
    */
}

#endif
