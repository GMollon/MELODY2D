#ifndef DEF_IO
#define DEF_IO
#include <sstream>



//********************************************//
//** READ SPY LINE ***************************//
//********************************************//

void Read_spy_line(
    vector<double>& New_Read,
    ifstream& Static_Control_file)
{
    string Keyword ;
    double b ;
    double bo ;
    double n ;
    double xmin, xmax, ymin, ymax ;

    Static_Control_file >> Keyword ;
    if (Keyword == "Position")
    {
        Static_Control_file >> Keyword ;
        if (Keyword == "X")
        {
            Static_Control_file >> b >> n ;
            New_Read = { 0, 0, b, n } ;
        }
        else if (Keyword == "Y")
        {
            Static_Control_file >> b >> n ;
            New_Read = { 0, 1, b, n } ;
        }
        else if (Keyword == "R")
        {
            Static_Control_file >> b >> n ;
            New_Read = { 0, 2, b, n } ;
        }
        else if (Keyword == "Length")
        {
            Static_Control_file >> b ;
            New_Read = { 0, 3, b } ;
        }
        else if (Keyword == "LengthInContact")
        {
            Static_Control_file >> b ;
            New_Read = { 0, 4, b } ;
        }
        else if (Keyword == "Area")
        {
            Static_Control_file >> b ;
            New_Read = { 0, 5, b } ;
        }
        else if (Keyword == "AllX")
        {
            New_Read = { 0, 6 } ;
        }
        else if (Keyword == "AllY")
        {
            New_Read = { 0, 7 } ;
        }
    }
    else if (Keyword == "Displacement")
    {
        Static_Control_file >> Keyword ;
        if (Keyword == "X")
        {
            Static_Control_file >> b >> n ;
            New_Read = { 1, 0, b, n } ;
        }
        else if (Keyword == "Y")
        {
            Static_Control_file >> b >> n ;
            New_Read = { 1, 1, b, n } ;
        }
        else if (Keyword == "Mag")
        {
            Static_Control_file >> b >> n ;
            New_Read = { 1, 2, b, n } ;
        }
        else if (Keyword == "R")
        {
            Static_Control_file >> b >> n ;
            New_Read = { 1, 3, b, n } ;
        }
    }
    else if (Keyword == "Velocity")
    {
        Static_Control_file >> Keyword ;
        if (Keyword == "X")
        {
            Static_Control_file >> b >> n ;
            New_Read = { 2, 0, b, n } ;
        }
        else if (Keyword == "Y")
        {
            Static_Control_file >> b >> n ;
            New_Read = { 2, 1, b, n } ;
        }
        else if (Keyword == "Mag")
        {
            Static_Control_file >> b >> n ;
            New_Read = { 2, 2, b, n } ;
        }
        else if (Keyword == "R")
        {
            Static_Control_file >> b >> n ;
            New_Read = { 2, 3, b, n } ;
        }
    }
    else if (Keyword == "Acceleration")
    {
        Static_Control_file >> Keyword ;
        if (Keyword == "X")
        {
            Static_Control_file >> b >> n ;
            New_Read = { 3, 0, b, n } ;
        }
        else if (Keyword == "Y")
        {
            Static_Control_file >> b >> n ;
            New_Read = { 3, 1, b, n } ;
        }
        else if (Keyword == "Mag")
        {
            Static_Control_file >> b >> n ;
            New_Read = { 3, 2, b, n } ;
        }
        else if (Keyword == "R")
        {
            Static_Control_file >> b >> n ;
            New_Read = { 3, 3, b, n } ;
        }
    }
    else if (Keyword == "Force")
    {
        Static_Control_file >> Keyword ;
        if (Keyword == "X")
        {
            Static_Control_file >> Keyword ;
            if (Keyword == "Total")
            {
                Static_Control_file >> b >> n ;
                New_Read = { 4, 0, 0, b, n } ;
            }
            else if (Keyword == "Internal")
            {
                Static_Control_file >> b >> n ;
                New_Read = { 4, 0, 1, b, n } ;
            }
            else if (Keyword == "Contact")
            {
                Static_Control_file >> b >> n ;
                New_Read = { 4, 0, 2, b, n } ;
            }
            else if (Keyword == "Body")
            {
                Static_Control_file >> b >> n ;
                New_Read = { 4, 0, 3, b, n } ;
            }
            else if (Keyword == "Dirichlet")
            {
                Static_Control_file >> b >> n ;
                New_Read = { 4, 0, 4, b, n } ;
            }
            else if (Keyword == "Neumann")
            {
                Static_Control_file >> b >> n ;
                New_Read = { 4, 0, 5, b, n } ;
            }
            else if (Keyword == "Damping")
            {
                Static_Control_file >> b >> n ;
                New_Read = { 4, 0, 6, b, n } ;
            }
            else if (Keyword == "Alid")
            {
                Static_Control_file >> b >> n ;
                New_Read = { 4, 0, 7, b, n } ;
            }
        }
        else if (Keyword == "Y")
        {
            Static_Control_file >> Keyword ;
            if (Keyword == "Total")
            {
                Static_Control_file >> b >> n ;
                New_Read = { 4, 1, 0, b, n } ;
            }
            else if (Keyword == "Internal")
            {
                Static_Control_file >> b >> n ;
                New_Read = { 4, 1, 1, b, n } ;
            }
            else if (Keyword == "Contact")
            {
                Static_Control_file >> b >> n ;
                New_Read = { 4, 1, 2, b, n } ;
            }
            else if (Keyword == "Body")
            {
                Static_Control_file >> b >> n ;
                New_Read = { 4, 1, 3, b, n } ;
            }
            else if (Keyword == "Dirichlet")
            {
                Static_Control_file >> b >> n ;
                New_Read = { 4, 1, 4, b, n } ;
            }
            else if (Keyword == "Neumann")
            {
                Static_Control_file >> b >> n ;
                New_Read = { 4, 1, 5, b, n } ;
            }
            else if (Keyword == "Damping")
            {
                Static_Control_file >> b >> n ;
                New_Read = { 4, 1, 6, b, n } ;
            }
            else if (Keyword == "Alid")
            {
                Static_Control_file >> b >> n ;
                New_Read = { 4, 1, 7, b, n } ;
            }
        }
        else if (Keyword == "Mag")
        {
            Static_Control_file >> Keyword ;
            if (Keyword == "Total")
            {
                Static_Control_file >> b >> n ;
                New_Read = { 4, 2, 0, b, n } ;
            }
            else if (Keyword == "Internal")
            {
                Static_Control_file >> b >> n ;
                New_Read = { 4, 2, 1, b, n } ;
            }
            else if (Keyword == "Contact")
            {
                Static_Control_file >> b >> n ;
                New_Read = { 4, 2, 2, b, n } ;
            }
            else if (Keyword == "Body")
            {
                Static_Control_file >> b >> n ;
                New_Read = { 4, 2, 3, b, n } ;
            }
            else if (Keyword == "Dirichlet")
            {
                Static_Control_file >> b >> n ;
                New_Read = { 4, 2, 4, b, n } ;
            }
            else if (Keyword == "Neumann")
            {
                Static_Control_file >> b >> n ;
                New_Read = { 4, 2, 5, b, n } ;
            }
            else if (Keyword == "Damping")
            {
                Static_Control_file >> b >> n ;
                New_Read = { 4, 2, 6, b, n } ;
            }
            else if (Keyword == "Alid")
            {
                Static_Control_file >> b >> n ;
                New_Read = { 4, 2, 7, b, n } ;
            }
        }
        else if (Keyword == "R")
        {
            Static_Control_file >> Keyword ;
            if (Keyword == "Total")
            {
                Static_Control_file >> b >> n ;
                New_Read = { 4, 3, 0, b, n } ;
            }
            else if (Keyword == "Internal")
            {
                Static_Control_file >> b >> n ;
                New_Read = { 4, 3, 1, b, n } ;
            }
            else if (Keyword == "Contact")
            {
                Static_Control_file >> b >> n ;
                New_Read = { 4, 3, 2, b, n } ;
            }
            else if (Keyword == "Body")
            {
                Static_Control_file >> b >> n ;
                New_Read = { 4, 3, 3, b, n } ;
            }
            else if (Keyword == "Dirichlet")
            {
                Static_Control_file >> b >> n ;
                New_Read = { 4, 3, 4, b, n } ;
            }
            else if (Keyword == "Neumann")
            {
                Static_Control_file >> b >> n ;
                New_Read = { 4, 3, 5, b, n } ;
            }
            else if (Keyword == "Damping")
            {
                Static_Control_file >> b >> n ;
                New_Read = { 4, 3, 6, b, n } ;
            }
            else if (Keyword == "Alid")
            {
                Static_Control_file >> b >> n ;
                New_Read = { 4, 3, 7, b, n } ;
            }
        }
    }
    else if (Keyword == "Jacobian")
    {
        Static_Control_file >> b >> n ;
        New_Read = { 5, b, n } ;
    }
    else if (Keyword == "Stress")
    {
        Static_Control_file >> Keyword ;
        if (Keyword == "XX")
        {
            Static_Control_file >> b >> n ;
            New_Read = { 6, 0, b, n } ;
        }
        else if (Keyword == "YY")
        {
            Static_Control_file >> b >> n ;
            New_Read = { 6, 1, b, n } ;
        }
        else if (Keyword == "XY")
        {
            Static_Control_file >> b >> n ;
            New_Read = { 6, 2, b, n } ;
        }
        else if (Keyword == "ZZ")
        {
            Static_Control_file >> b >> n ;
            New_Read = { 6, 3, b, n } ;
        }
        else if (Keyword == "Tresca")
        {
            Static_Control_file >> b >> n ;
            New_Read = { 6, 4, b, n } ;
        }
        else if (Keyword == "VM")
        {
            Static_Control_file >> b >> n ;
            New_Read = { 6, 5, b, n } ;
        }
        else if (Keyword == "Major")
        {
            Static_Control_file >> b >> n ;
            New_Read = { 6, 6, b, n } ;
        }
        else if (Keyword == "Intermediate")
        {
            Static_Control_file >> b >> n ;
            New_Read = { 6, 7, b, n } ;
        }
        else if (Keyword == "Minor")
        {
            Static_Control_file >> b >> n ;
            New_Read = { 6, 8, b, n } ;
        }
        else if (Keyword == "Spherical")
        {
            Static_Control_file >> b >> n ;
            New_Read = { 6, 9, b, n } ;
        }
    }
    else if (Keyword == "Contact")
    {
        Static_Control_file >> Keyword ;
        if (Keyword == "Gapn")
        {
            Static_Control_file >> b >> bo >> n ;
            New_Read = { 7, 0, b, bo, n } ;
        }
        else if (Keyword == "Gapt")
        {
            Static_Control_file >> b >> bo >> n ;
            New_Read = { 7, 1, b, bo, n } ;
        }
        else if (Keyword == "Xsi")
        {
            Static_Control_file >> b >> bo >> n ;
            New_Read = { 7, 2, b, bo, n } ;
        }
        else if (Keyword == "Xnorm")
        {
            Static_Control_file >> b >> bo >> n ;
            New_Read = { 7, 3, b, bo, n } ;
        }
        else if (Keyword == "Ynorm")
        {
            Static_Control_file >> b >> bo >> n ;
            New_Read = { 7, 4, b, bo, n } ;
        }
        else if (Keyword == "Damage")
        {
            Static_Control_file >> b >> bo >> n ;
            New_Read = { 7, 5, b, bo, n } ;
        }
        else if (Keyword == "Fx")
        {
            Static_Control_file >> b >> bo >> n ;
            New_Read = { 7, 6, b, bo, n } ;
        }
        else if (Keyword == "Fy")
        {
            Static_Control_file >> b >> bo >> n ;
            New_Read = { 7, 7, b, bo, n } ;
        }
        else if (Keyword == "Length")
        {
            Static_Control_file >> b >> bo >> n ;
            New_Read = { 7, 8, b, bo, n } ;
        }
        else if (Keyword == "Master")
        {
            Static_Control_file >> b >> bo >> n ;
            New_Read = { 7, 9, b, bo, n } ;
        }
        else if (Keyword == "Px")
        {
            Static_Control_file >> b >> bo >> n ;
            New_Read = { 7, 10, b, bo, n } ;
        }
        else if (Keyword == "Py")
        {
            Static_Control_file >> b >> bo >> n ;
            New_Read = { 7, 11, b, bo, n } ;
        }
        else if (Keyword == "Slave")
        {
            Static_Control_file >> b >> bo >> n ;
            New_Read = { 7, 12, b, bo, n } ;
        }
        else if (Keyword == "Xslave")
        {
            Static_Control_file >> b >> bo >> n ;
            New_Read = { 7, 13, b, bo, n } ;
        }
        else if (Keyword == "Yslave")
        {
            Static_Control_file >> b >> bo >> n ;
            New_Read = { 7, 14, b, bo, n } ;
        }
    }
    else if (Keyword == "Damage")
    {
        Static_Control_file >> b ;
        New_Read = { 8, b } ;
    }
    else if (Keyword == "Energy")
    {
        Static_Control_file >> Keyword >> b ;
        if (Keyword == "Kinetic")
        {
            New_Read = { 9, 0, b } ;
        }
        else if (Keyword == "Deformation")
        {
            New_Read = { 9, 1, b } ;
        }
    }
    else if (Keyword == "Work")
    {
        Static_Control_file >> Keyword >> b ;
        if (Keyword == "Internal")
        {
            New_Read = { 10, 0, b } ;
        }
        else if (Keyword == "Contact")
        {
            New_Read = { 10, 1, b } ;
        }
        else if (Keyword == "Body")
        {
            New_Read = { 10, 2, b } ;
        }
        else if (Keyword == "Dirichlet")
        {
            New_Read = { 10, 3, b } ;
        }
        else if (Keyword == "Neumann")
        {
            New_Read = { 10, 4, b } ;
        }
        else if (Keyword == "Damping")
        {
            New_Read = { 10, 5, b } ;
        }
        else if (Keyword == "Alid")
        {
            New_Read = { 10, 6, b } ;
        }
    }
    else if (Keyword == "Error")
    {
        Static_Control_file >> Keyword ;
        if (Keyword == "X")
        {
            Static_Control_file >> b >> n ;
            New_Read = { 11, 0, b, n } ;
        }
        else if (Keyword == "Y")
        {
            Static_Control_file >> b >> n ;
            New_Read = { 11, 1, b, n } ;
        }
        else if (Keyword == "Mag")
        {
            Static_Control_file >> b >> n ;
            New_Read = { 11, 2, b, n } ;
        }
        else if (Keyword == "Norm")
        {
            Static_Control_file >> b >> n ;
            New_Read = { 11, 3, b, n } ;
        }
        else if (Keyword == "Max")
        {
            Static_Control_file >> b >> n ;
            New_Read = { 11, 4, b, n } ;
        }
    }
    else if (Keyword == "MassScaling")
    {
        Static_Control_file >> Keyword ;
        if (Keyword == "X")
        {
            Static_Control_file >> b >> n ;
            New_Read = { 12, 0, b, n } ;
        }
        else if (Keyword == "Y")
        {
            Static_Control_file >> b >> n ;
            New_Read = { 12, 1, b, n } ;
        }
        else if (Keyword == "Moy")
        {
            Static_Control_file >> b >> n ;
            New_Read = { 12, 2, b, n } ;
        }
        else if (Keyword == "dX")
        {
            Static_Control_file >> b >> n ;
            New_Read = { 12, 3, b, n } ;
        }
        else if (Keyword == "dY")
        {
            Static_Control_file >> b >> n ;
            New_Read = { 12, 4, b, n } ;
        }
        else if (Keyword == "dMoy")
        {
            Static_Control_file >> b >> n ;
            New_Read = { 12, 5, b, n } ;
        }
        else if (Keyword == "Xmass")
        {
            Static_Control_file >> b >> n ;
            New_Read = { 12, 6, b, n } ;
        }
        else if (Keyword == "Ymass")
        {
            Static_Control_file >> b >> n ;
            New_Read = { 12, 7, b, n } ;
        }
        else if (Keyword == "Mass")
        {
            Static_Control_file >> b >> n ;
            New_Read = { 12, 8, b, n } ;
        }

        else if (Keyword == "Polar")
        {
            Static_Control_file >> Keyword >> n >> xmin >> xmax >> ymin >> ymax ;
            if (Keyword == "Contact_Normal")
            {
                New_Read = { 11, 0, n, xmin, xmax, ymin, ymax } ;
            }
        }
    }
}



//********************************************//
//** LOAD STATIC *****************************//
//********************************************//

void Load_static(
    string& Simulation_name,
    int& Nb_materials,
    vector<Material>& Materials,
    int& Nb_contact_laws,
    vector<Contact_law>& Contact_laws,
    vector<vector<int>>& Contacts_Table,
    string& Solver,
    double& Tini,
    double& Deltat,
    double& Tend,
    double& Target_error,
    double& Inv_Target_error,
    double& Control_parameter,
    double& Accepted_ratio,
    double& Max_mass_scaling,
    double& Control_parameter_mass_scaling,
    double& Error_factor_mass_scaling,
    double& Decrease_factor_mass_scaling,
    double& Save_period,
    double& Print_period,
    double& Contact_update_period,
    double& Xmin_period,
    double& Xmax_period,
    double& Penalty,
    double& Xgravity,
    double& Ygravity,
    int& Activate_plot,
    double& Xmin_plot,
    double& Xmax_plot,
    double& Ymin_plot,
    double& Ymax_plot,
    int& Nb_monitored,
    vector<vector<double>>& Monitored,
    int& Nb_deactivated,
    vector<vector<double>>& Deactivated,
    int& Nb_spies,
    vector<Spy>& Spies,
    int& Nb_regions,
    vector<vector<int>>& Regions,
    int& Nb_bodies,
    vector<Body>& Bodies,
    vector<int>& To_Plot )
{
    cout << endl ;
    cout << "Loading static control" << endl ;
    int index_body = 0 ;
    int nb(0), ndef(0), nr(0), nn(0), nd(0), ng(0) ;
    ifstream Static_Control_file ("STATIC_CONTROL.asc") ;
    string line ;
    string token ;
    int flag_body_rigid(0) ;
    while(getline(Static_Control_file, line))
    {
        if (line.substr(0,15)=="SIMULATION_NAME")
        {
            Static_Control_file >> Simulation_name ;
            getline(Static_Control_file, token) ;
            cout << "Simulation name " << Simulation_name << endl ;
        }

        if (line.substr(0,9)=="MATERIALS")
        {
            Static_Control_file >> Nb_materials ;
            getline(Static_Control_file, token) ;
            getline(Static_Control_file, token) ;
            for (int i(0) ; i < Nb_materials ; i++)
            {
                string name ;
                Static_Control_file >> name ;
                getline(Static_Control_file, token) ;
                string type ;
                Static_Control_file >> type ;
                getline(Static_Control_file, token) ;
                vector<double> parameters({0}) ;
                if (type=="ElasticLinear")
                {
                    double p1, p2, p3, p4, p5, p6 ;
                    Static_Control_file >> p1 >> p2 >> p3 >> p4 >> p5 >> p6 ;
                    getline(Static_Control_file, token) ;
                    vector<double> parameters({p1, p2, p3, p4 / ( 2. * (1. + p5 )), p4 * p5 / (( 1. + p5 ) * ( 1. - 2. * p5 )), p6}) ;
                    Material mat ( i, name, type, parameters ) ;
                    Materials.push_back(mat) ;
                }
                else if (type=="NeoHookean")
                {
                    double p1, p2, p3, p4, p5, p6 ;
                    Static_Control_file >> p1 >> p2 >> p3 >> p4 >> p5 >> p6 ;
                    getline(Static_Control_file, token) ;
                    vector<double> parameters({p1, p2, p3, p4 / ( 2. * (1. + p5 )), p4 / ( 3. * (1. - 2. * p5 )), p6}) ;
                    Material mat ( i, name, type, parameters ) ;
                    Materials.push_back(mat) ;
                }
                // OTHER MATERIALS ??
                getline(Static_Control_file, token) ;
            }
            cout << "Materials " << Nb_materials << endl ;
            for (int i=0 ; i<Nb_materials ; i++)
            {
                cout << i << ' ' << Materials[i].name << ' ' << Materials[i].type << ' ' << Materials[i].parameters[0] << endl ;
            }
        }

        if (line.substr(0,12)=="CONTACT_LAWS")
        {
            Static_Control_file >> Nb_contact_laws ;
            getline(Static_Control_file, token) ;
            getline(Static_Control_file, token) ;
            for (int i(0) ; i < Nb_contact_laws ; i++)
            {
                string material1 ;
                Static_Control_file >> material1 ;
                getline(Static_Control_file, token) ;
                string material2 ;
                Static_Control_file >> material2 ;
                getline(Static_Control_file, token) ;
                string type ;
                Static_Control_file >> type ;
                getline(Static_Control_file, token) ;
                string length_evolution ;
                Static_Control_file >> length_evolution ;
                getline(Static_Control_file, token) ;
                vector<double> parameters({0}) ;
                if (type=="Frictionless")
                {
                    double p1 ;
                    Static_Control_file >> p1 ;
                    getline(Static_Control_file, token) ;
                    vector<double> parameters({p1}) ;
                    Contact_law law ( i, material1, material2, type, length_evolution, parameters ) ;
                    Contact_laws.push_back(law) ;
                }
                else if (type=="Cohesion")
                {
                    double p1, p2, p3, p4 ;
                    Static_Control_file >> p1 >> p2 >> p3 >> p4 ;
                    getline(Static_Control_file, token) ;
                    vector<double> parameters({p1, p2, p3, p4}) ;
                    Contact_law law ( i, material1, material2, type, length_evolution, parameters ) ;
                    Contact_laws.push_back(law) ;
                }
                else if (type=="MohrCoulomb")
                {
                    double p1, p2, p3, p4, p5 ;
                    Static_Control_file >> p1 >> p2 >> p3 >> p4 >> p5 ;
                    getline(Static_Control_file, token) ;
                    vector<double> parameters({p1, p2, p3, p4, p5}) ;
                    Contact_law law ( i, material1, material2, type, length_evolution, parameters ) ;
                    Contact_laws.push_back(law) ;
                }
                else if (type=="DampedMohrCoulomb")
                {
                    double p1, p2, p3, p4, p5, p6 ;
                    Static_Control_file >> p1 >> p2 >> p3 >> p4 >> p5 >> p6 ;
                    getline(Static_Control_file, token) ;
                    vector<double> parameters({p1, p2, p3, p4, p5, p6}) ;
                    Contact_law law ( i, material1, material2, type, length_evolution, parameters ) ;
                    Contact_laws.push_back(law) ;
                    cout << "damp " << p6 << endl ;
                }
                else if (type=="TwoSlopesMohrCoulomb")
                {
                    double p1, p2, p3, p4, p5, p6 ;
                    Static_Control_file >> p1 >> p2 >> p3 >> p4 >> p5 >> p6 ;
                    getline(Static_Control_file, token) ;
                    vector<double> parameters({p1, p2, p3, p4, p5, p6}) ;
                    Contact_law law ( i, material1, material2, type, length_evolution, parameters ) ;
                    Contact_laws.push_back(law) ;
                }
                else if (type=="CZMlinear")
                {
                    double p1, p2, p3, p4, p5, p6 ;
                    Static_Control_file >> p1 >> p2 >> p3 >> p4 >> p5 >> p6 ;
                    getline(Static_Control_file, token) ;
                    vector<double> parameters({p1, p2, p3, p4, p5, p6}) ;
                    Contact_law law ( i, material1, material2, type, length_evolution, parameters ) ;
                    Contact_laws.push_back(law) ;
                }
                else if (type=="CZMfatigue")
                {
                    double p1, p2, p3, p4, p5, p6, p7, p8, p9 ;
                    Static_Control_file >> p1 >> p2 >> p3 >> p4 >> p5 >> p6 >> p7 >> p8 >> p9 ;
                    getline(Static_Control_file, token) ;
                    vector<double> parameters({p1, p2, p3, p4, p5, p6, p7, p8, p9}) ;
                    Contact_law law ( i, material1, material2, type, length_evolution, parameters ) ;
                    Contact_laws.push_back(law) ;
                }
                else if (type=="BondedMohrCoulomb")
                {
                    double p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11 ;
                    Static_Control_file >> p1 >> p2 >> p3 >> p4 >> p5 >> p6 >> p7 >> p8 >> p9 >> p10 >> p11 ;
                    getline(Static_Control_file, token) ;
                    vector<double> parameters({p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11}) ;
                    Contact_law law ( i, material1, material2, type, length_evolution, parameters ) ;
                    Contact_laws.push_back(law) ;
                }
                // OTHER LAWS ?
                getline(Static_Control_file, token) ;
            }
            vector<int> c(Nb_materials) ;
            for (int i=0 ; i<Nb_materials ; i++)
                c[i] = -1 ;
            for (int i=0 ; i<Nb_materials ; i++)
                Contacts_Table.push_back(c) ;//Contacts_Table[i] = c ;
            int nn(0), pp(0) ;
            for (int i=0 ; i<Nb_contact_laws ; i++)
            {
                for (int n(0) ; n<Nb_materials ; n++)
                    if (Materials[n].name == Contact_laws[i].material1 )
                        nn = n ;
                for (int p(0) ; p<Nb_materials ; p++)
                    if (Materials[p].name == Contact_laws[i].material2 )
                        pp = p ;
                Contacts_Table[nn][pp] = i ;
                Contacts_Table[pp][nn] = i ;
            }
            cout << "Contact laws " << Nb_contact_laws << endl ;
            for (int i=0 ; i<Nb_contact_laws ; i++)
            {
                cout << i << ' ' << Contact_laws[i].material1 << ' ' << Contact_laws[i].material2 << ' ' << Contact_laws[i].type << ' ' << Contact_laws[i].length_evolution << ' ' << Contact_laws[i].parameters[0] << endl ;
            }
        }

        if (line.substr(0,6)=="SOLVER")
        {
            Static_Control_file >> Solver ;
            getline(Static_Control_file, token) ;


            if(Solver == "Euler")
            {
                Static_Control_file >> Tini >> Deltat >> Tend ;
                getline(Static_Control_file, token) ;
                Static_Control_file >> Target_error >> Control_parameter >> Accepted_ratio ;
                getline(Static_Control_file, token) ;
                Static_Control_file >> Save_period >> Print_period >> Contact_update_period ;
                getline(Static_Control_file, token) ;
                cout << "Solver " << Solver << ' ' << Tini << ' ' << Deltat << ' ' << Tend << ' ' << Target_error << ' ' << Control_parameter << ' ' << Accepted_ratio << endl ;
                cout << Save_period << ' ' << Print_period << ' ' << Contact_update_period << endl ;
                Inv_Target_error = 1 / Target_error ;
            }
            if(Solver == "Euler_2")
            {
                Static_Control_file >> Tini >> Deltat >> Tend ;
                getline(Static_Control_file, token) ;
                Static_Control_file >> Save_period >> Print_period >> Contact_update_period ;
                getline(Static_Control_file, token) ;
                cout << "Solver " << Solver << ' ' << Tini << ' ' << Deltat << ' ' << Tend << endl ;
                cout << Save_period << ' ' << Print_period << ' ' << Contact_update_period << endl ;
            }
            else if(Solver == "Adaptive_Euler")
            {
                Static_Control_file >> Tini >> Deltat >> Tend ;
                getline(Static_Control_file, token) ;
                Static_Control_file >> Target_error >> Control_parameter >> Accepted_ratio ;
                getline(Static_Control_file, token) ;
                Static_Control_file >> Save_period >> Print_period >> Contact_update_period ;
                getline(Static_Control_file, token) ;
                cout << "Solver " << Solver << ' ' << Tini << ' ' << Deltat << ' ' << Tend << ' ' << Target_error << ' ' << Control_parameter << ' ' << Accepted_ratio << endl ;
                cout << Save_period << ' ' << Print_period << ' ' << Contact_update_period << endl ;
                Inv_Target_error = 1 / Target_error ;
            }
            else if(Solver == "Mass_Scaling")
            {
                Static_Control_file >> Tini >> Deltat >> Tend ;
                getline(Static_Control_file, token) ;
                Static_Control_file >> Target_error >> Control_parameter >> Accepted_ratio ;
                getline(Static_Control_file, token) ;
                Static_Control_file >> Max_mass_scaling >> Control_parameter_mass_scaling >> Error_factor_mass_scaling >> Decrease_factor_mass_scaling ;
                getline(Static_Control_file, token) ;
                Static_Control_file >> Save_period >> Print_period >> Contact_update_period ;
                getline(Static_Control_file, token) ;
                cout << "Solver " << Solver << ' ' << Tini << ' ' << Deltat << ' ' << Tend << ' ' << Target_error << ' ' << Control_parameter << ' ' << Accepted_ratio << endl ;
                cout << "Mass scaling " << Max_mass_scaling << ' ' << Control_parameter_mass_scaling << ' ' << Error_factor_mass_scaling << ' ' << Decrease_factor_mass_scaling << endl ;
                cout << Save_period << ' ' << Print_period << ' ' << Contact_update_period << endl ;
                Inv_Target_error = 1 / Target_error ;
            }
            else if(Solver == "Adaptive_Mass_Scaling")
            {
                Static_Control_file >> Tini >> Deltat >> Tend ;
                getline(Static_Control_file, token) ;
                Static_Control_file >> Target_error >> Control_parameter >> Accepted_ratio ;
                getline(Static_Control_file, token) ;
                Static_Control_file >> Max_mass_scaling >> Control_parameter_mass_scaling >> Error_factor_mass_scaling >> Decrease_factor_mass_scaling ;
                getline(Static_Control_file, token) ;
                Static_Control_file >> Save_period >> Print_period >> Contact_update_period ;
                getline(Static_Control_file, token) ;
                cout << "Solver " << Solver << ' ' << Tini << ' ' << Deltat << ' ' << Tend << ' ' << Target_error << ' ' << Control_parameter << ' ' << Accepted_ratio << endl ;
                cout << "Mass scaling " << Max_mass_scaling << ' ' << Control_parameter_mass_scaling << ' ' << Error_factor_mass_scaling << ' ' << Decrease_factor_mass_scaling << endl ;
                cout << Save_period << ' ' << Print_period << ' ' << Contact_update_period << endl ;
                Inv_Target_error = 1 / Target_error ;
            }




        }

        if (line.substr(0,19)=="PERIODIC_BOUNDARIES")
        {
            Static_Control_file >> Xmin_period >> Xmax_period ;
            getline(Static_Control_file, token) ;
            cout << "Periodic boundaries " << Xmin_period << ' ' << Xmax_period << endl ;
        }

        if (line.substr(0,17)=="PENALTY_PARAMETER")
        {
            Static_Control_file >> Penalty ;
            getline(Static_Control_file, token) ;
            cout << "Penalty " << Penalty << endl ;
        }

        if (line.substr(0,7)=="GRAVITY")
        {
            Static_Control_file >> Xgravity >> Ygravity ;
            getline(Static_Control_file, token) ;
            cout << "Gravity " << Xgravity << ' ' << Ygravity << endl ;
        }

        if (line.substr(0,8)=="PLOTTING")
        {
            Static_Control_file >> Activate_plot ;
            getline(Static_Control_file, token) ;
            Static_Control_file >> Xmin_plot >> Xmax_plot >> Ymin_plot >> Ymax_plot ;
            getline(Static_Control_file, token) ;
        }

        if (line.substr(0,13)=="NUMBER_BODIES")
        {
            Static_Control_file >> Nb_bodies ;
            getline(Static_Control_file, token) ;
            cout << "Nb bodies " << Nb_bodies << endl ;
        }

        //
        if (line.substr(0,10)=="MONITORING")
        {
            Static_Control_file >> Nb_monitored ;
            getline(Static_Control_file, token) ;
            cout << "Nb monitoring " << Nb_monitored << endl ;
            vector<double> New_Read ;
            for (int i(0) ; i < Nb_monitored ; i++)
            {
                Read_spy_line(New_Read, Static_Control_file) ;
                Monitored.push_back(New_Read) ;
                getline(Static_Control_file, token) ;
            }
        }

        if (line.substr(0,12)=="DEACTIVATION")
        {
            Static_Control_file >> Nb_deactivated ;
            getline(Static_Control_file, token) ;
            cout << "Nb possible deactivations " << Nb_deactivated << endl ;
            double threshold ;
            string Keyword ;
            vector<double> New_Read ;
            for (int i(0) ; i < Nb_deactivated ; i++)
            {
                Read_spy_line(New_Read, Static_Control_file) ;
                Static_Control_file >> Keyword >> threshold ;
                if (Keyword == "below")
                    New_Read.push_back(0.) ;
                else if (Keyword == "equal")
                    New_Read.push_back(1.) ;
                else if (Keyword == "above")
                    New_Read.push_back(2.) ;
                New_Read.push_back(threshold) ;
                Deactivated.push_back(New_Read) ;
                getline(Static_Control_file, token) ;
            }
        }

        if (line.substr(0,5)=="SPIES")
        {
            Static_Control_file >> Nb_spies ;
            getline(Static_Control_file, token) ;
            cout << "Nb spies " << Nb_spies << endl ;
            vector<double> New_Read ;
            for (int i(0) ; i < Nb_spies ; i++)
            {
                string filename ;
                int    nb_quantities ;
                double period ;
                vector<vector<double>> quantities ;
                getline(Static_Control_file, token) ;
                Static_Control_file >> filename >> nb_quantities >> period ;
                getline(Static_Control_file, token) ;
                //cout << "spy " << i << " " << filename << " " << nb_quantities << " " << period << endl ;
                for (int j(0) ; j < nb_quantities ; j++)
                {
                    Read_spy_line(New_Read, Static_Control_file) ;
                    //cout << New_Read[0] << ' ' << New_Read[1] << ' ' << New_Read[2] << endl ;
                    quantities.push_back(New_Read) ;
                    getline(Static_Control_file, token) ;
                    //cout << "line " << j << ' ' << New_Read[0] << endl ;
                }
                Spy newspy(filename, nb_quantities, period, quantities) ;
                Spies.push_back(newspy) ;
            }
        }

        if (line.substr(0,7)=="GRAPHIC")
        {
            int j ;

            for (int i(0) ; i < 41 ; i++)
            {
                Static_Control_file >> j >> token ;
                getline(Static_Control_file, token) ;
                To_Plot[i] = j ;
            }
        }

        if (line.substr(0,10)=="DEFORMABLE")
        {
            nb++ ;
            ndef++ ;
            flag_body_rigid = 0 ;
            Static_Control_file >> index_body ;
            getline(Static_Control_file, token) ;
            string mat ;
            Static_Control_file >> mat ;
            getline(Static_Control_file, token) ;
            string periodicity ;
            Static_Control_file >> periodicity ;
            getline(Static_Control_file, token) ;
            double nodal_distance, detection_distance, contact_distance ;
            Static_Control_file >> nodal_distance >> detection_distance >> contact_distance ;
            getline(Static_Control_file, token) ;
            Body body ( index_body, mat, periodicity, "deformable" ) ;
            body.nodal_distance=nodal_distance ;
            body.detection_distance=detection_distance ;
            body.contact_distance=contact_distance ;
            for (int i(0) ; i<Nb_materials ; i++)
            {
                if (Materials[i].name == mat)
                {
                    body.density=Materials[i].parameters[0] ;
                    body.material_index=i ;
                    body.heat_capacity=Materials[i].parameters[5] ;
                    break ;
                }
            }
            Bodies.push_back(body) ;
        }

        if (line.substr(0,5)=="RIGID")
        {
            nb++ ;
            nr++ ;
            flag_body_rigid = 1 ;
            Static_Control_file >> index_body ;
            getline(Static_Control_file, token) ;
            string mat ;
            Static_Control_file >> mat ;
            getline(Static_Control_file, token) ;
            string periodicity ;
            Static_Control_file >> periodicity ;
            getline(Static_Control_file, token) ;
            double nodal_distance, detection_distance, contact_distance ;
            Static_Control_file >> nodal_distance >> detection_distance >> contact_distance ;
            getline(Static_Control_file, token) ;
            Body body( index_body, mat, periodicity, "rigid" ) ;
            body.nodal_distance=nodal_distance ;
            body.detection_distance=detection_distance ;
            body.contact_distance=contact_distance ;
            for (int i(0) ; i<Nb_materials ; i++)
            {
                if (Materials[i].name == mat)
                {
                    body.density=Materials[i].parameters[0] ;
                    body.material_index=i ;
                    body.heat_capacity=Materials[i].parameters[5] ;
                    break ;
                }
            }
            Bodies.push_back(body) ;
        }

        if (line.substr(0,7)=="BORDERS")
        {
            int nb_borders ;
            vector<Border> borders ;
            Static_Control_file >> nb_borders ;
            getline(Static_Control_file, token) ;
            getline(Static_Control_file, token) ;
            string periodicity, interpolant, direction, type1, type2 ;
            double stiffness, damping ;
            int nb_border_nodes ;
            for (int i(0) ; i < nb_borders ; i++)
            {
                // Set list of border nodes //
                Static_Control_file >> periodicity >> interpolant ;
                getline(Static_Control_file, token) ;
                Static_Control_file >> nb_border_nodes ;
                getline(Static_Control_file, token) ;
                Border border ( nb_border_nodes, periodicity, index_body, i) ;
                border.interpolant = interpolant ;

                //
                int number_xdirichlet_values = 0 ;
                vector<double> xdirichlet_instants ;
                vector<double> xdirichlet_values ;
                vector<double> xdirichlet_parameters ;
                int number_xneumann_values = 0 ;
                vector<double> xneumann_instants ;
                vector<double> xneumann_values ;
                Static_Control_file >> direction >> type1 ;
                if (type1=="None")
                {
                    border.x_bc_type="none" ;
                    getline(Static_Control_file, token) ;
                }
                else if (type1=="Dirichlet")
                {
                    Static_Control_file >> type2 ;
                    if (type2=="Soft")
                    {
                        border.x_bc_type="soft" ;
                        Static_Control_file >> stiffness >> damping ;
                        xdirichlet_parameters = {stiffness, damping} ;
                        getline(Static_Control_file, token) ;
                    }
                    else if (type2=="Driven")
                    {
                        border.x_bc_type="driven" ;
                        getline(Static_Control_file, token) ;
                    }
                    //getline(Static_Control_file, token) ;
                    Static_Control_file >> number_xdirichlet_values ;
                    getline(Static_Control_file, token) ;
                    if (number_xdirichlet_values>0)
                    {
                        double xdirichlet_instant ;
                        double xdirichlet_value ;
                        for (int j(0) ; j < number_xdirichlet_values ; j++)
                        {
                            Static_Control_file >> xdirichlet_instant >> xdirichlet_value ;
                            getline(Static_Control_file, token) ;
                            xdirichlet_instants.push_back(xdirichlet_instant) ;
                            xdirichlet_values.push_back(xdirichlet_value) ;
                        }
                    }
                }
                else if (type1=="Neumann")
                {
                    Static_Control_file >> type2 ;
                    if (type2=="Oriented")
                    {
                        border.x_bc_type="oriented" ;
                        getline(Static_Control_file, token) ;
                    }
                    else if (type2=="Following")
                    {
                        border.x_bc_type="following" ;
                        getline(Static_Control_file, token) ;
                    }
                    //getline(Static_Control_file, token) ;
                    Static_Control_file >> number_xneumann_values ;
                    getline(Static_Control_file, token) ;
                    if (number_xneumann_values>0)
                    {
                        double xneumann_instant ;
                        double xneumann_value ;
                        for (int j(0) ; j < number_xneumann_values ; j++)
                        {
                            Static_Control_file >> xneumann_instant >> xneumann_value ;
                            getline(Static_Control_file, token) ;
                            xneumann_instants.push_back(xneumann_instant) ;
                            xneumann_values.push_back(xneumann_value) ;
                        }
                    }
                    border.number_xneumann_values = number_xneumann_values ;
                    border.xneumann_instants = xneumann_instants ;
                    border.xneumann_values = xneumann_values ;
                }
                border.number_xdirichlet_values = number_xdirichlet_values ;
                border.xdirichlet_instants = xdirichlet_instants ;
                border.xdirichlet_values = xdirichlet_values ;
                border.xdirichlet_parameters = xdirichlet_parameters ;
                border.number_xneumann_values = number_xneumann_values ;
                border.xneumann_instants = xneumann_instants ;
                border.xneumann_values = xneumann_values ;

                int number_ydirichlet_values = 0 ;
                vector<double> ydirichlet_instants ;
                vector<double> ydirichlet_values ;
                vector<double> ydirichlet_parameters ;
                int number_yneumann_values = 0 ;
                vector<double> yneumann_instants ;
                vector<double> yneumann_values ;
                Static_Control_file >> direction >> type1 ;
                if (type1=="None")
                {
                    border.y_bc_type="none" ;
                    getline(Static_Control_file, token) ;
                }
                else if (type1=="Dirichlet")
                {
                    Static_Control_file >> type2 ;
                    if (type2=="Soft")
                    {
                        border.y_bc_type="soft" ;
                        Static_Control_file >> stiffness >> damping ;
                        ydirichlet_parameters = {stiffness, damping} ;
                        getline(Static_Control_file, token) ;
                    }
                    else if (type2=="Driven")
                    {
                        border.y_bc_type="driven" ;
                        getline(Static_Control_file, token) ;
                    }
                    //getline(Static_Control_file, token) ;
                    Static_Control_file >> number_ydirichlet_values ;
                    getline(Static_Control_file, token) ;
                    if (number_ydirichlet_values>0)
                    {
                        double ydirichlet_instant ;
                        double ydirichlet_value ;
                        for (int j(0) ; j < number_ydirichlet_values ; j++)
                        {
                            Static_Control_file >> ydirichlet_instant >> ydirichlet_value ;
                            getline(Static_Control_file, token) ;
                            ydirichlet_instants.push_back(ydirichlet_instant) ;
                            ydirichlet_values.push_back(ydirichlet_value) ;
                        }
                    }
                }
                else if (type1=="Neumann")
                {
                    Static_Control_file >> type2 ;
                    if (type2=="Oriented")
                    {
                        border.y_bc_type="oriented" ;
                        getline(Static_Control_file, token) ;
                    }
                    else if (type2=="Following")
                    {
                        border.y_bc_type="following" ;
                        getline(Static_Control_file, token) ;
                    }
                    //getline(Static_Control_file, token) ;
                    Static_Control_file >> number_yneumann_values ;
                    getline(Static_Control_file, token) ;
                    if (number_yneumann_values>0)
                    {
                        double yneumann_instant ;
                        double yneumann_value ;
                        for (int j(0) ; j < number_yneumann_values ; j++)
                        {
                            Static_Control_file >> yneumann_instant >> yneumann_value ;
                            getline(Static_Control_file, token) ;
                            yneumann_instants.push_back(yneumann_instant) ;
                            yneumann_values.push_back(yneumann_value) ;
                        }
                    }
                }
                border.number_ydirichlet_values = number_ydirichlet_values ;
                border.ydirichlet_instants = ydirichlet_instants ;
                border.ydirichlet_values = ydirichlet_values ;
                border.ydirichlet_parameters = ydirichlet_parameters ;
                border.number_yneumann_values = number_yneumann_values ;
                border.yneumann_instants = yneumann_instants ;
                border.yneumann_values = yneumann_values ;

                /*
                // Set history of Dirichlet BC along x //
                getline(Static_Control_file, token) ;
                if (token.substr(0,6)=="Driven") border.x_bc_type="driven" ;
                int number_xdirichlet_values ;
                vector<double> xdirichlet_instants ;
                vector<double> xdirichlet_values ;
                Static_Control_file >> number_xdirichlet_values ;
                getline(Static_Control_file, token) ;
                if (number_xdirichlet_values>0)
                {
                	double xdirichlet_instant ;
                	double xdirichlet_value ;
                	for (int j(0) ; j < number_xdirichlet_values ; j++)
                	{
                		Static_Control_file >> xdirichlet_instant >> xdirichlet_value ;
                		getline(Static_Control_file, token) ;
                		xdirichlet_instants.push_back(xdirichlet_instant) ;
                		xdirichlet_values.push_back(xdirichlet_value) ;
                	}
                }
                border.number_xdirichlet_values = number_xdirichlet_values ;
                border.xdirichlet_instants = xdirichlet_instants ;
                border.xdirichlet_values = xdirichlet_values ;

                // Set history of Dirichlet BC along y //
                getline(Static_Control_file, token) ;
                if (token.substr(0,6)=="Driven") border.y_bc_type="driven" ;
                int number_ydirichlet_values ;
                vector<double> ydirichlet_instants ;
                vector<double> ydirichlet_values ;
                Static_Control_file >> number_ydirichlet_values ;
                getline(Static_Control_file, token) ;
                if (number_ydirichlet_values>0)
                {
                	double ydirichlet_instant ;
                	double ydirichlet_value ;
                	for (int j(0) ; j < number_ydirichlet_values ; j++)
                	{
                		Static_Control_file >> ydirichlet_instant >> ydirichlet_value ;
                		getline(Static_Control_file, token) ;
                		ydirichlet_instants.push_back(ydirichlet_instant) ;
                		ydirichlet_values.push_back(ydirichlet_value) ;
                	}
                }
                border.number_ydirichlet_values = number_ydirichlet_values ;
                border.ydirichlet_instants = ydirichlet_instants ;
                border.ydirichlet_values = ydirichlet_values ;

                // Set history of Neumann BC along x //
                getline(Static_Control_file, token) ;
                int number_xneumann_values ;
                vector<double> xneumann_instants ;
                vector<double> xneumann_values ;
                Static_Control_file >> number_xneumann_values ;
                getline(Static_Control_file, token) ;
                if (number_xneumann_values>0)
                {
                	double xneumann_instant ;
                	double xneumann_value ;
                	for (int j(0) ; j < number_xneumann_values ; j++)
                	{
                		Static_Control_file >> xneumann_instant >> xneumann_value ;
                		getline(Static_Control_file, token) ;
                		xneumann_instants.push_back(xneumann_instant) ;
                		xneumann_values.push_back(xneumann_value) ;
                	}
                }
                border.number_xneumann_values = number_xneumann_values ;
                border.xneumann_instants = xneumann_instants ;
                border.xneumann_values = xneumann_values ;

                // Set history of Neumann BC along y //
                getline(Static_Control_file, token) ;
                int number_yneumann_values ;
                vector<double> yneumann_instants ;
                vector<double> yneumann_values ;
                Static_Control_file >> number_yneumann_values ;
                getline(Static_Control_file, token) ;
                if (number_yneumann_values>0)
                {
                	double yneumann_instant ;
                	double yneumann_value ;
                	for (int j(0) ; j < number_yneumann_values ; j++)
                	{
                		Static_Control_file >> yneumann_instant >> yneumann_value ;
                		getline(Static_Control_file, token) ;
                		yneumann_instants.push_back(yneumann_instant) ;
                		yneumann_values.push_back(yneumann_value) ;
                	}
                }
                border.number_yneumann_values = number_yneumann_values ;
                border.yneumann_instants = yneumann_instants ;
                border.yneumann_values = yneumann_values ;
                */

                // Add to list of borders //
                borders.push_back(border) ;
                getline(Static_Control_file, token) ;
            }
            Bodies[index_body].borders=borders ;
            Bodies[index_body].stored_borders=borders ;
            Bodies[index_body].nb_borders=nb_borders ;
            for (int j=0 ; j<nb_borders ; j++)
                Bodies[index_body].nb_border_nodes += borders[j].number_border_nodes ;
            for (int i(0) ; i < nb_borders ; i++)
            {
                cout << "Border " << i << endl ;
                cout << Bodies[index_body].borders[i].number_xdirichlet_values << " X dirichlet instants" << endl ;
                for (int j(0) ; j<Bodies[index_body].borders[i].number_xdirichlet_values ; j++)
                {
                    cout << Bodies[index_body].borders[i].xdirichlet_instants[j] << ' ' << Bodies[index_body].borders[i].xdirichlet_values[j] << endl ;
                }
                cout << Bodies[index_body].borders[i].number_ydirichlet_values << " Y dirichlet instants" << endl ;
                for (int j(0) ; j<Bodies[index_body].borders[i].number_ydirichlet_values ; j++)
                {
                    cout << Bodies[index_body].borders[i].ydirichlet_instants[j] << ' ' << Bodies[index_body].borders[i].ydirichlet_values[j] << endl ;
                }
                cout << Bodies[index_body].borders[i].number_xneumann_values << " X neumann instants" << endl ;
                for (int j(0) ; j<Bodies[index_body].borders[i].number_xneumann_values ; j++)
                {
                    cout << Bodies[index_body].borders[i].xneumann_instants[j] << ' ' << Bodies[index_body].borders[i].xneumann_values[j] << endl ;
                }
                cout << Bodies[index_body].borders[i].number_yneumann_values << " Y neumann instants" << endl ;
                for (int j(0) ; j<Bodies[index_body].borders[i].number_yneumann_values ; j++)
                {
                    cout << Bodies[index_body].borders[i].yneumann_instants[j] << ' ' << Bodies[index_body].borders[i].yneumann_values[j] << endl ;
                }
            }
        }

        if (line.substr(0,4)=="ALID")
        {
            int number_xalid_values, number_yalid_values ;
            double xmin, xmax, ymin, ymax, range, alphain, alphaout, betain, betaout, exponent ;
            Static_Control_file >> number_xalid_values >> number_yalid_values >> xmin >> xmax >> ymin >> ymax >> range >> alphain >> alphaout >> betain >> betaout >> exponent ;
            cout << number_xalid_values << ' ' << number_yalid_values << ' ' << xmin << ' ' << xmax << ' ' << ymin << ' ' << ymax << ' ' << range << endl ;
            getline(Static_Control_file, token) ;
            getline(Static_Control_file, token) ;
            Bodies[index_body].alid_xmin = xmin ;
            Bodies[index_body].alid_xmax = xmax ;
            Bodies[index_body].alid_ymin = ymin ;
            Bodies[index_body].alid_ymax = ymax ;
            Bodies[index_body].alid_range = range ;
            Bodies[index_body].alid_alphain = alphain ;
            Bodies[index_body].alid_alphaout = alphaout ;
            Bodies[index_body].alid_betain = betain ;
            Bodies[index_body].alid_betaout = betaout ;
            Bodies[index_body].alid_exponent = exponent ;
            if ((number_xalid_values==0) & (number_yalid_values==0))
            {
                Bodies[index_body].flag_alid = 0 ;
            }
            else if ((number_xalid_values>0) & (number_yalid_values==0))
            {
                vector<double> xalid_instants ;
                double xalid_instant ;
                vector<double> xalid_values ;
                double xalid_value ;
                for (int j(0) ; j < number_xalid_values ; j++)
                {
                    Static_Control_file >> xalid_instant >> xalid_value ;
                    getline(Static_Control_file, token) ;
                    xalid_instants.push_back(xalid_instant) ;
                    xalid_values.push_back(xalid_value) ;
                }
                Bodies[index_body].flag_alid = 1 ;
                Bodies[index_body].xalid_instants = xalid_instants ;
                Bodies[index_body].xalid_values = xalid_values ;
            }
            else if ((number_xalid_values==0) & (number_yalid_values>0))
            {
                vector<double> yalid_instants ;
                double yalid_instant ;
                vector<double> yalid_values ;
                double yalid_value ;
                for (int j(0) ; j < number_yalid_values ; j++)
                {
                    Static_Control_file >> yalid_instant >> yalid_value ;
                    getline(Static_Control_file, token) ;
                    yalid_instants.push_back(yalid_instant) ;
                    yalid_values.push_back(yalid_value) ;
                }
                Bodies[index_body].flag_alid = 1 ;
                Bodies[index_body].yalid_instants = yalid_instants ;
                Bodies[index_body].yalid_values = yalid_values ;
            }
            else if ((number_xalid_values>0) & (number_yalid_values>0))
            {
                vector<double> xalid_instants ;
                double xalid_instant ;
                vector<double> xalid_values ;
                double xalid_value ;
                vector<double> yalid_instants ;
                double yalid_instant ;
                vector<double> yalid_values ;
                double yalid_value ;
                for (int j(0) ; j < number_xalid_values ; j++)
                {
                    Static_Control_file >> xalid_instant >> xalid_value ;
                    getline(Static_Control_file, token) ;
                    xalid_instants.push_back(xalid_instant) ;
                    xalid_values.push_back(xalid_value) ;
                    cout << xalid_instant << ' ' << xalid_value << endl ;
                }
                getline(Static_Control_file, token) ;
                for (int j(0) ; j < number_yalid_values ; j++)
                {
                    Static_Control_file >> yalid_instant >> yalid_value ;
                    getline(Static_Control_file, token) ;
                    yalid_instants.push_back(yalid_instant) ;
                    yalid_values.push_back(yalid_value) ;
                    cout << yalid_instant << ' ' << yalid_value << endl ;
                }
                Bodies[index_body].flag_alid = 1 ;
                Bodies[index_body].xalid_instants = xalid_instants ;
                Bodies[index_body].xalid_values = xalid_values ;
                Bodies[index_body].yalid_instants = yalid_instants ;
                Bodies[index_body].yalid_values = yalid_values ;
            }
        }

    }
    Static_Control_file.close () ;

    cout << endl ;
    cout << "Loading static data" << endl ;
    ifstream Static_Data_file ("STATIC_DATA.asc") ;
    while(getline(Static_Data_file, line))
    {
        if (line.substr(0,10)=="DEFORMABLE")
        {
            flag_body_rigid = 0 ;
            Static_Data_file >> index_body ;
            getline(Static_Data_file, token) ;
            cout << "Deformable " << index_body << endl ;
        }

        if (line.substr(0,5)=="RIGID")
        {
            flag_body_rigid = 1 ;
            nd += 3 ;
            Static_Data_file >> index_body ;
            getline(Static_Data_file, token) ;
            double x, y, m, i, invm, invi ;
            Static_Data_file >> x >> y ;
            getline(Static_Data_file, token) ;
            Static_Data_file >> m >> i >> invm >> invi ;
            getline(Static_Data_file, token) ;
            Bodies[index_body].x_initial = x ;
            Bodies[index_body].y_initial = y ;
            Bodies[index_body].r_initial = 0. ;
            Bodies[index_body].mass = m ;
            Bodies[index_body].inertia = i ;
            Bodies[index_body].inverse_mass = invm ;
            Bodies[index_body].inverse_inertia = invi ;
            cout << "Rigid " << index_body << endl ;
        }

        if (line.substr(0,5)=="NODES")
        {
            int nb_nodes ;
            vector<Node> nodes ;
            Static_Data_file >> nb_nodes ;
            nn += nb_nodes ;
            getline(Static_Data_file, token) ;
            getline(Static_Data_file, token) ;
            if (flag_body_rigid==0)
            {
                nd += 2 * nb_nodes ;
                double j, ka, kb, temp, xini, yini, d, mx, my, imx, imy ;
                for (int i(0) ; i < nb_nodes ; i++)
                {
                    Static_Data_file >> j
                                     >> xini
                                     >> yini
                                     >> d
                                     >> mx
                                     >> my
                                     >> imx
                                     >> imy ;
                    getline(Static_Data_file, token) ;
                    Node node ( xini, yini, d, mx, my, imx, imy ) ;
                    //
                    node.alid_alpha = 0. ;
                    node.alid_beta = 0. ;
                    if (Bodies[index_body].flag_alid==1)
                    {
                        if (yini < Bodies[index_body].alid_ymin + Bodies[index_body].alid_range)
                        {
                            ka = (Bodies[index_body].alid_alphaout - Bodies[index_body].alid_alphain) / pow ( -Bodies[index_body].alid_range, Bodies[index_body].alid_exponent) ;
                            kb = (Bodies[index_body].alid_betaout - Bodies[index_body].alid_betain) / pow ( -Bodies[index_body].alid_range, Bodies[index_body].alid_exponent) ;
                            node.alid_alpha = Bodies[index_body].alid_alphain + ka * pow ( yini - ( Bodies[index_body].alid_ymin + Bodies[index_body].alid_range ), Bodies[index_body].alid_exponent ) ;
                            node.alid_beta = Bodies[index_body].alid_betain + kb * pow ( yini - ( Bodies[index_body].alid_ymin + Bodies[index_body].alid_range ), Bodies[index_body].alid_exponent ) ;
                        }
                        else if (yini > Bodies[index_body].alid_ymax - Bodies[index_body].alid_range)
                        {
                            ka = (Bodies[index_body].alid_alphaout - Bodies[index_body].alid_alphain) / pow ( Bodies[index_body].alid_range, Bodies[index_body].alid_exponent) ;
                            kb = (Bodies[index_body].alid_betaout - Bodies[index_body].alid_betain) / pow ( Bodies[index_body].alid_range, Bodies[index_body].alid_exponent) ;
                            node.alid_alpha = Bodies[index_body].alid_alphain + ka * pow ( yini - ( Bodies[index_body].alid_ymax - Bodies[index_body].alid_range ), Bodies[index_body].alid_exponent ) ;
                            node.alid_beta = Bodies[index_body].alid_betain + kb * pow ( yini - ( Bodies[index_body].alid_ymax - Bodies[index_body].alid_range ), Bodies[index_body].alid_exponent ) ;
                        }
                        if (xini < Bodies[index_body].alid_xmin + Bodies[index_body].alid_range)
                        {
                            ka = (Bodies[index_body].alid_alphaout - Bodies[index_body].alid_alphain) / pow( -Bodies[index_body].alid_range, Bodies[index_body].alid_exponent) ;
                            kb = (Bodies[index_body].alid_betaout - Bodies[index_body].alid_betain) / pow( -Bodies[index_body].alid_range, Bodies[index_body].alid_exponent) ;
                            temp = Bodies[index_body].alid_alphain + ka * pow ( xini - ( Bodies[index_body].alid_xmin + Bodies[index_body].alid_range ), Bodies[index_body].alid_exponent ) ;
                            if (temp > node.alid_alpha)
                                node.alid_alpha = temp ;
                            temp = Bodies[index_body].alid_betain + kb * pow ( xini - ( Bodies[index_body].alid_xmin + Bodies[index_body].alid_range ), Bodies[index_body].alid_exponent ) ;
                            if (temp > node.alid_beta)
                                node.alid_beta = temp ;
                        }
                        else if (xini > Bodies[index_body].alid_xmax - Bodies[index_body].alid_range)
                        {
                            ka = (Bodies[index_body].alid_alphaout - Bodies[index_body].alid_alphain) / pow( Bodies[index_body].alid_range, Bodies[index_body].alid_exponent) ;
                            kb = (Bodies[index_body].alid_betaout - Bodies[index_body].alid_betain) / pow( Bodies[index_body].alid_range, Bodies[index_body].alid_exponent) ;
                            temp = Bodies[index_body].alid_alphain + ka * pow ( xini - ( Bodies[index_body].alid_xmax - Bodies[index_body].alid_range ), Bodies[index_body].alid_exponent ) ;
                            if (temp > node.alid_alpha)
                                node.alid_alpha = temp ;
                            temp = Bodies[index_body].alid_betain + kb * pow ( xini - ( Bodies[index_body].alid_xmax - Bodies[index_body].alid_range ), Bodies[index_body].alid_exponent ) ;
                            if (temp > node.alid_beta)
                                node.alid_beta = temp ;
                        }
                        //cout << index_body << ' ' << node.alid_parameter << endl ;
                    }
                    //
                    /*
                    if (Bodies[index_body].flag_alid==1)
                    {
                        if (Bodies[index_body].alid_ymin<yini & yini<Bodies[index_body].alid_ymax)
                        {
                            double k = (Bodies[index_body].alid_betamax-Bodies[index_body].alid_betamin)*pow((Bodies[index_body].alid_ymax-Bodies[index_body].alid_ymin),-Bodies[index_body].alid_exponent);
                            if (Bodies[index_body].alid_betamax>Bodies[index_body].alid_betamin)    node.alid_parameter = Bodies[index_body].alid_betamin + k * pow( ( yini - Bodies[index_body].alid_ymin ) , Bodies[index_body].alid_exponent ) ;
                            else                                                                    node.alid_parameter = Bodies[index_body].alid_betamax - k * pow( ( Bodies[index_body].alid_ymax - yini ) , Bodies[index_body].alid_exponent ) ;
                        }
                        else node.alid_parameter = 0. ;
                        //cout << index_body << ' ' << node.alid_parameter << endl ;
                    }
                    */
                    nodes.push_back(node) ;
                }
                Bodies[index_body].nodes=nodes ;
                Bodies[index_body].stored_nodes=nodes ;
                Bodies[index_body].nb_nodes=nb_nodes ;
            }
            else if (flag_body_rigid==1)
            {
                double j, xini, yini ;
                for (int i(0) ; i < nb_nodes ; i++)
                {
                    Static_Data_file >> j
                                     >> xini
                                     >> yini ;
                    getline(Static_Data_file, token) ;
                    Node node ( xini, yini, 0, 0, 0, 0, 0 ) ;
                    node.inverse_distance_estimator = 1. / Bodies[index_body].nodal_distance ;
                    node.x_internal_force = 0. ;
                    node.y_internal_force = 0. ;
                    nodes.push_back(node) ;
                }
                Bodies[index_body].nodes=nodes ;
                Bodies[index_body].nb_nodes=nb_nodes ;
            }
        }

        if (line.substr(0,7)=="BORDERS")
        {
            int nb_borders ;
            int border_before ;
            int border_after ;
            Static_Data_file >> nb_borders ;
            getline(Static_Data_file, token) ;
            getline(Static_Data_file, token) ;
            int nb_border_nodes ;
            for (int i(0) ; i < nb_borders ; i++)
            {
                // Set list of border nodes //
                getline(Static_Data_file, token) ;
                Static_Data_file >> nb_border_nodes >> border_before >> border_after ;
                getline(Static_Data_file, token) ;
                vector<int> list_nodes ;
                for (int j(0) ; j < nb_border_nodes ; j++)
                {
                    int n ;
                    Static_Data_file >> n ;
                    getline(Static_Data_file, token) ;
                    list_nodes.push_back(n) ;
                }
                Bodies[index_body].borders[i].border_nodes = list_nodes ;
                Bodies[index_body].borders[i].border_before = border_before ;
                Bodies[index_body].borders[i].border_after = border_after ;
                getline(Static_Data_file, token) ;
            }
            Bodies[index_body].Initialize_segments(Xmin_period, Xmax_period) ;
        }

        if (line.substr(0,11)=="INTEGRATION")
        {
            string integration ;
            vector<Gauss> gpoints ;
            vector<vector<int>> region_gpoints ;
            int nb_gauss, nb_regions, nb_per_region ;
            getline(Static_Data_file, integration) ;
            if (integration.substr(0,5)=="Gauss")
            {
                Static_Data_file >> nb_gauss >> nb_regions ;
                ng += nb_gauss ;
                getline(Static_Data_file, token) ;
                getline(Static_Data_file, token) ;
                for (int i(0) ; i < nb_gauss ; i++)
                {
                    getline(Static_Data_file, token) ;
                    double xgauss, ygauss, weight, jacobian ;
                    Static_Data_file >> xgauss >> ygauss >> weight >> jacobian ;
                    getline(Static_Data_file, token) ;
                    int nb_influencing ;
                    Static_Data_file >> nb_influencing ;
                    getline(Static_Data_file, token) ;
                    Gauss gauss ( xgauss, ygauss, weight, jacobian, nb_influencing ) ;
                    vector<int> influencing_nodes ;
                    vector<double> shape_functions ;
                    vector<double> shape_xderiv ;
                    vector<double> shape_yderiv ;
                    for (int j(0) ; j < nb_influencing ; j++)
                    {
                        int influencing_node, n1, n2 ;
                        double shape_function, xderiv, yderiv ;
                        Static_Data_file >> influencing_node >> n1 >> n2 >> shape_function >> xderiv >> yderiv ;
                        getline(Static_Data_file, token) ;
                        influencing_nodes.push_back(influencing_node) ;
                        shape_functions.push_back(shape_function) ;
                        shape_xderiv.push_back(xderiv) ;
                        shape_yderiv.push_back(yderiv) ;
                    }
                    gauss.influencing_nodes = influencing_nodes ;
                    gauss.shape_functions = shape_functions ;
                    gauss.shape_xderiv = shape_xderiv ;
                    gauss.shape_yderiv = shape_yderiv ;
                    gpoints.push_back(gauss) ;
                    getline(Static_Data_file, token) ;
                }
            }
            else if (integration.substr(0,5)=="Nodal")
            {
                Static_Data_file >> nb_gauss >> nb_regions ;
                ng += nb_gauss ;
                getline(Static_Data_file, token) ;
                getline(Static_Data_file, token) ;
                for (int i(0) ; i < nb_gauss ; i++)
                {
                    getline(Static_Data_file, token) ;
                    double xgauss, ygauss, weight ;
                    Static_Data_file >> xgauss >> ygauss >> weight ;
                    getline(Static_Data_file, token) ;
                    int nb_influencing ;
                    Static_Data_file >> nb_influencing ;
                    getline(Static_Data_file, token) ;
                    Gauss gauss ( xgauss, ygauss, weight, 1., nb_influencing ) ;
                    vector<int> influencing_nodes ;
                    vector<double> shape_xderiv ;
                    vector<double> shape_yderiv ;
                    for (int j(0) ; j < nb_influencing ; j++)
                    {
                        int influencing_node, n1, n2 ;
                        double xderiv, yderiv ;
                        Static_Data_file >> influencing_node >> n1 >> n2 >> xderiv >> yderiv ;
                        getline(Static_Data_file, token) ;
                        influencing_nodes.push_back(influencing_node) ;
                        shape_xderiv.push_back(xderiv) ;
                        shape_yderiv.push_back(yderiv) ;
                    }
                    gauss.influencing_nodes = influencing_nodes ;
                    gauss.shape_xderiv = shape_xderiv ;
                    gauss.shape_yderiv = shape_yderiv ;
                    gpoints.push_back(gauss) ;
                    getline(Static_Data_file, token) ;
                }
            }
            nb_per_region = ceil( (double)nb_gauss / (double)nb_regions ) ;
            for (int i=0 ; i<nb_regions ; i++)
            {
                vector<int> per_region ;
                for (int j=0 ; j<nb_per_region ; j++)
                {
                    if (i*nb_per_region+j<nb_gauss)
                        per_region.push_back(i*nb_per_region+j) ;
                }
                region_gpoints.push_back(per_region) ;
                Nb_regions++ ;
                Regions.push_back({index_body, i}) ;
            }
            Bodies[index_body].gpoints=gpoints ;
            Bodies[index_body].nb_gauss=nb_gauss ;
            Bodies[index_body].nb_regions=nb_regions ;
            Bodies[index_body].region_gpoints=region_gpoints ;

            vector<double> v(Bodies[index_body].nb_regions) ;
            for (int i=0 ; i<nb_regions ; i++)
                v[i]=0.;
            Bodies[index_body].Elastic_energy = v ;
            for (int i(0) ; i< Bodies[index_body].nb_nodes ; i++)
            {
                Bodies[index_body].nodes[i].x_regions_internal_forces = v ;
                Bodies[index_body].nodes[i].y_regions_internal_forces = v ;
            }
        }

        if (line.substr(0,13)=="TRIANGULATION")
        {
            int nb_cells ;
            vector<vector<int>> cells ;
            Static_Data_file >> nb_cells ;
            getline(Static_Data_file, token) ;
            for (int i(0) ; i < nb_cells ; i++)
            {
                int n0, n1, n2, n3 ;
                Static_Data_file >> n0 >> n1 >> n2 >> n3 ;
                getline(Static_Data_file, token) ;
                vector<int> cell({ n1, n2, n3 }) ;
                cells.push_back(cell) ;
            }
            getline(Static_Data_file, token) ;
            Bodies[index_body].triangulation=cells ;
            Bodies[index_body].nb_cells=nb_cells ;
        }

        if (line.substr(0,11)=="DOF_TO_DISP")
        {
            int nb_nodes = Bodies[index_body].nb_nodes ;
            for (int i(0) ; i < nb_nodes ; i++)
            {
                getline(Static_Data_file, token) ;
                getline(Static_Data_file, token) ;
                int nb_influencing ;
                Static_Data_file >> nb_influencing ;
                getline(Static_Data_file, token) ;
                vector<int> influencing_nodes ;
                vector<double> shape_functions ;
                vector<double> shape_xderiv ;
                vector<double> shape_yderiv ;
                for (int j(0) ; j < nb_influencing ; j++)
                {
                    int nn ;
                    double ss ;
                    double sx ;
                    double sy ;
                    Static_Data_file >> nn >> ss >> sx >> sy ;
                    getline(Static_Data_file, token) ;
                    influencing_nodes.push_back(nn) ;
                    shape_functions.push_back(ss) ;
                    shape_xderiv.push_back(sx) ;
                    shape_yderiv.push_back(sy) ;
                }
                Bodies[index_body].nodes[i].number_influencing_nodes = nb_influencing ;
                Bodies[index_body].nodes[i].inverse_distance_estimator = 1. / ( Bodies[index_body].nodes[i].domain_size * sqrt ( 3.1415927/4. / nb_influencing ) ) ;
                Bodies[index_body].nodes[i].influencing_nodes = influencing_nodes ;
                Bodies[index_body].nodes[i].shape_functions = shape_functions ;
                Bodies[index_body].nodes[i].shape_xderiv = shape_xderiv ;
                Bodies[index_body].nodes[i].shape_yderiv = shape_yderiv ;
            }
        }

        if (line.substr(0,8)=="MATRICES")
        {
            int nb_matrix ;
            vector<vector<int>> matrix_coordinates ;
            vector<double> stiffness ;
            vector<double> damping ;
            Static_Data_file >> nb_matrix ;
            getline(Static_Data_file, token) ;
            getline(Static_Data_file, token) ;
            for (int i(0) ; i < nb_matrix ; i++)
            {
                int n1, n2 ;
                double s, d ;
                Static_Data_file >> n1 >> n2 >> s >> d ;
                getline(Static_Data_file, token) ;
                vector<int> nn({ n1, n2 }) ;
                matrix_coordinates.push_back(nn) ;
                stiffness.push_back(s) ;
                damping.push_back(d) ;
            }
            Bodies[index_body].matrix_coordinates=matrix_coordinates ;
            Bodies[index_body].stiffness=stiffness ;
            Bodies[index_body].damping=damping ;
            Bodies[index_body].nb_matrix=nb_matrix ;
        }
    }
    Static_Data_file.close () ;

    cout << "    " << nb << " bodies" << endl ;
    cout << "    " << nr << " rigid" << endl ;
    cout << "    " << ndef << " deformable" << endl ;
    cout << "    " << Nb_regions << " regions" << endl ;
    cout << "    " << nn << " nodes" << endl ;
    cout << "    " << nd << " degrees of freedom" << endl ;
    cout << "    " << ng << " Gauss points" << endl ;
}




//********************************************//
//** LOAD DYNAMIC ****************************//
//********************************************//

void Load_dynamic(
    int Number_save,
    int& Number_print,
    double& Time,
    int& Number_iteration,
    double& Deltat,
    string& Solver,
    double& Next_save,
    double& Next_print,
    double& Next_contact_update,
    double& Xmin_period,
    double& Xmax_period,
    vector<Body>& Bodies,
    vector<int>& flags )
{
    cout << endl ;
    cout << "Loading dynamic data of file " << Number_save << endl ;
    int index_body = 0 ;
    int flag_body_rigid(0) ;
    //flags = { 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 } ;
    string filename ;
    stringstream sfilename ;
    sfilename << Number_save ;
    if      (Number_save<10)
        filename="DYNAMIC_0000"+sfilename.str()+".asc" ;
    else if (Number_save<100)
        filename="DYNAMIC_000"+sfilename.str()+".asc" ;
    else if (Number_save<1000)
        filename="DYNAMIC_00"+sfilename.str()+".asc" ;
    else if (Number_save<10000)
        filename="DYNAMIC_0"+sfilename.str()+".asc" ;
    else if (Number_save<100000)
        filename="DYNAMIC_"+sfilename.str()+".asc" ;
    ifstream Dynamic_file (filename) ;
    string line ;
    string token ;
    string status ;
    while(getline(Dynamic_file, line))
    {
        if( flags[7] == 0 )
            cout << line << endl ;
        if (line.substr(0,13)=="KILL_VELOCITY")
            flags[0] = 1 ;
        if (line.substr(0,14)=="MONITOR_ENERGY")
            flags[1] = 1 ;
        if (line.substr(0,18)=="MONITOR_BOUNDARIES")
            flags[2] = 1 ;
        if (line.substr(0,14)=="INITIALIZE_CZM")
            flags[3] = 1 ;
        if (line.substr(0,18)=="UPDATE_MASS_MATRIX")
            flags[4] = 1 ;
        if (line.substr(0,23)=="UPDATE_STIFFNESS_MATRIX")
            flags[5] = 1 ;
        if (line.substr(0,21)=="UPDATE_DAMPING_MATRIX")
            flags[6] = 1 ;
        if (line.substr(0,6) =="NO_LOG")
            flags[7] = 1 ;
        if (line.substr(0,13)=="NO_MONITORING")
            flags[8] = 1 ;
        if (line.substr(0,15)=="NO_SELF_CONTACT")
            flags[9] = 1 ;
        if (line.substr(0,10)=="RESET_WORK")
            flags[10] = 1 ;


        if (line.substr(0,6)=="SOLVER")
        {
            Dynamic_file >> Solver ;
            getline(Dynamic_file, token) ;
            Dynamic_file >> Number_save >> Number_iteration >> Time >> Number_print ;
            getline(Dynamic_file, token) ;
            Dynamic_file >> Deltat ;
            getline(Dynamic_file, token) ;
            Dynamic_file >> Next_save >> Next_print >> Next_contact_update ;
            getline(Dynamic_file, token) ;
            cout << "Solver " << Solver << ' ' << Number_save << ' ' << Number_iteration << ' ' << Time << ' ' << Number_print << endl ;
            cout << Deltat << ' ' << Next_save << ' ' << Next_print << ' ' << Next_contact_update << endl ;
        }

        if (line.substr(0,10)=="DEFORMABLE")
        {
            Dynamic_file >> status ;
            getline(Dynamic_file, token) ;
            Dynamic_file >> index_body ;
            getline(Dynamic_file, token) ;
            if( flags[7] == 0 )
                cout << "Deformable " << index_body << endl ;
            Bodies[index_body].status = status ;
            flag_body_rigid = 0 ;
        }

        if (line.substr(0,5)=="RIGID")
        {
            Dynamic_file >> status ;
            getline(Dynamic_file, token) ;
            Dynamic_file >> index_body ;
            getline(Dynamic_file, token) ;
            if( flags[7] == 0 )
                cout << "Rigid " << index_body << endl ;
            Bodies[index_body].status = status ;
            flag_body_rigid = 1 ;
        }

        if (line.substr(0,10)=="KINEMATICS")
        {
            if (flag_body_rigid==0)
            {
                int nb_nodes ;
                Dynamic_file >> nb_nodes ;
                getline(Dynamic_file, token) ;
                int j ;
                double xcurrent, ycurrent ;
                double xdisp, ydisp, xpdisp, ypdisp ;
                double xvel, yvel, xpvel, ypvel ;
                double xacc, yacc, xpacc, ypacc ;
                for (int i(0) ; i < nb_nodes ; i++)
                {
                    Dynamic_file >> j >> xcurrent >> ycurrent ;
                    Dynamic_file >> xdisp >> ydisp ;
                    Dynamic_file >> xvel >> yvel ;
                    Dynamic_file >> xacc >> yacc ;
                    Dynamic_file >> xpdisp >> ypdisp ;
                    Dynamic_file >> xpvel >> ypvel ;
                    Dynamic_file >> xpacc >> ypacc ;
                    getline(Dynamic_file, token) ;
                    Bodies[index_body].nodes[i].x_displacement = xdisp ;
                    Bodies[index_body].nodes[i].y_displacement = ydisp ;
                    Bodies[index_body].nodes[i].x_displacement_parameter = xpdisp ;
                    Bodies[index_body].nodes[i].y_displacement_parameter = ypdisp ;
                    if ( flags[0] == 0 )
                    {
                        Bodies[index_body].nodes[i].x_velocity = xvel ;
                        Bodies[index_body].nodes[i].y_velocity = yvel ;
                        Bodies[index_body].nodes[i].x_velocity_parameter = xpvel ;
                        Bodies[index_body].nodes[i].y_velocity_parameter = ypvel ;
                    }
                    else
                    {
                        Bodies[index_body].nodes[i].x_velocity = 0. ;
                        Bodies[index_body].nodes[i].y_velocity = 0. ;
                        Bodies[index_body].nodes[i].x_velocity_parameter = 0. ;
                        Bodies[index_body].nodes[i].y_velocity_parameter = 0. ;
                    }
                    Bodies[index_body].nodes[i].x_acceleration = xacc ;
                    Bodies[index_body].nodes[i].y_acceleration = yacc ;
                    Bodies[index_body].nodes[i].x_acceleration_parameter = xpacc ;
                    Bodies[index_body].nodes[i].y_acceleration_parameter = ypacc ;
                    Bodies[index_body].nodes[i].Update_current_positions() ;
                }
            }
            else if (flag_body_rigid==1)
            {
                double xcurrent, ycurrent, rcurrent ;
                double xdisp, ydisp, rdisp ;
                double xvel, yvel, rvel ;
                Dynamic_file >> xcurrent >> ycurrent >> rcurrent ;
                getline(Dynamic_file, token) ;
                Dynamic_file >> xdisp >> ydisp >> rdisp ;
                getline(Dynamic_file, token) ;
                Dynamic_file >> xvel >> yvel >> rvel ;
                getline(Dynamic_file, token) ;
                Bodies[index_body].x_current = xcurrent ;
                Bodies[index_body].y_current = ycurrent ;
                Bodies[index_body].r_current = rcurrent ;
                Bodies[index_body].x_displacement = xdisp ;
                Bodies[index_body].y_displacement = ydisp ;
                Bodies[index_body].r_displacement = rdisp ;
                if ( flags[0] == 0 )
                {
                    Bodies[index_body].x_velocity = xvel ;
                    Bodies[index_body].y_velocity = yvel ;
                    Bodies[index_body].r_velocity = rvel ;
                }
                else
                {
                    Bodies[index_body].x_velocity = 0. ;
                    Bodies[index_body].y_velocity = 0. ;
                    Bodies[index_body].r_velocity = 0. ;
                }
                Bodies[index_body].Update_current_positions() ;
            }
            //cout << "Kinematics" << endl ;
        }

        if (line.substr(0,6)=="FORCES")
        {
            if (flag_body_rigid==0)
            {
                int j ;
                double x_internal_force, y_internal_force ;
                double x_contact_force, y_contact_force ;
                double x_self_contact_force, y_self_contact_force ;
                double x_body_force, y_body_force ;
                double x_dirichlet_force, y_dirichlet_force ;
                double x_neumann_force, y_neumann_force ;
                double x_damping_force, y_damping_force ;
                for (int i(0) ; i < Bodies[index_body].nb_nodes ; i++)
                {
                    Dynamic_file >> j >> x_internal_force >> y_internal_force ;
                    Dynamic_file >> x_contact_force >> y_contact_force ;
                    Dynamic_file >> x_self_contact_force >> y_self_contact_force ;
                    Dynamic_file >> x_body_force >> y_body_force ;
                    Dynamic_file >> x_dirichlet_force >> y_dirichlet_force ;
                    Dynamic_file >> x_neumann_force >> y_neumann_force ;
                    Dynamic_file >> x_damping_force >> y_damping_force ;
                    getline(Dynamic_file, token) ;
                    Bodies[index_body].nodes[i].x_internal_force = x_internal_force ;
                    Bodies[index_body].nodes[i].y_internal_force = y_internal_force ;
                    Bodies[index_body].nodes[i].x_contact_force = x_contact_force ;
                    Bodies[index_body].nodes[i].y_contact_force = y_contact_force ;
                    Bodies[index_body].nodes[i].x_self_contact_force = x_self_contact_force ;
                    Bodies[index_body].nodes[i].y_self_contact_force = y_self_contact_force ;
                    Bodies[index_body].nodes[i].x_body_force = x_body_force ;
                    Bodies[index_body].nodes[i].y_body_force = y_body_force ;
                    Bodies[index_body].nodes[i].x_dirichlet_force = x_dirichlet_force ;
                    Bodies[index_body].nodes[i].y_dirichlet_force = y_dirichlet_force ;
                    Bodies[index_body].nodes[i].x_neumann_force = x_neumann_force ;
                    Bodies[index_body].nodes[i].y_neumann_force = y_neumann_force ;
                    Bodies[index_body].nodes[i].x_damping_force = x_damping_force ;
                    Bodies[index_body].nodes[i].y_damping_force = y_damping_force ;
                    Bodies[index_body].nodes[i].Sum_up_forces() ;
                }
            }
            else if (flag_body_rigid==1)
            {
                double x_contact_force, y_contact_force, r_contact_force ;
                double x_body_force, y_body_force, r_body_force ;
                double x_dirichlet_force, y_dirichlet_force, r_dirichlet_force ;
                double x_neumann_force, y_neumann_force, r_neumann_force ;
                double x_damping_force, y_damping_force, r_damping_force ;
                Dynamic_file >> x_contact_force >> y_contact_force >> r_contact_force ;
                getline(Dynamic_file, token) ;
                Dynamic_file >> x_body_force >> y_body_force >> r_body_force ;
                getline(Dynamic_file, token) ;
                Dynamic_file >> x_dirichlet_force >> y_dirichlet_force >> r_dirichlet_force ;
                getline(Dynamic_file, token) ;
                Dynamic_file >> x_neumann_force >> y_neumann_force >> r_neumann_force ;
                getline(Dynamic_file, token) ;
                Dynamic_file >> x_damping_force >> y_damping_force >> r_damping_force ;
                getline(Dynamic_file, token) ;
                getline(Dynamic_file, token) ;
                Bodies[index_body].x_contact_force = x_contact_force ;
                Bodies[index_body].y_contact_force = y_contact_force ;
                Bodies[index_body].r_contact_force = r_contact_force ;
                Bodies[index_body].x_body_force = x_body_force ;
                Bodies[index_body].y_body_force = y_body_force ;
                Bodies[index_body].r_body_force = r_body_force ;
                Bodies[index_body].x_dirichlet_force = x_dirichlet_force ;
                Bodies[index_body].y_dirichlet_force = y_dirichlet_force ;
                Bodies[index_body].r_dirichlet_force = r_dirichlet_force ;
                Bodies[index_body].x_neumann_force = x_neumann_force ;
                Bodies[index_body].y_neumann_force = y_neumann_force ;
                Bodies[index_body].r_neumann_force = r_neumann_force ;
                Bodies[index_body].x_damping_force = x_damping_force ;
                Bodies[index_body].y_damping_force = y_damping_force ;
                Bodies[index_body].r_damping_force = r_damping_force ;
                Bodies[index_body].Sum_up_forces() ;
            }
            //cout << "Forces" << endl ;
        }

        if (line.substr(0,5)=="WORKS")
        {
            double w0, w1, w2, w3, w4, w5, w6 ;
            Dynamic_file >> w0 >> w1 >> w2 >> w3 >> w4 >> w5 >> w6 ;
            getline(Dynamic_file, token) ;
            if ( flags[10] == 0 )
            {
                Bodies[index_body].internal_work = w0 ;
                Bodies[index_body].contact_work = w1 ;
                Bodies[index_body].body_work = w2 ;
                Bodies[index_body].dirichlet_work = w3 ;
                Bodies[index_body].neumann_work = w4 ;
                Bodies[index_body].damping_work = w5 ;
                Bodies[index_body].alid_work = w6 ;
            }
            else
            {
                Bodies[index_body].internal_work = 0. ;
                Bodies[index_body].contact_work = 0. ;
                Bodies[index_body].body_work = 0. ;
                Bodies[index_body].dirichlet_work = 0. ;
                Bodies[index_body].neumann_work = 0. ;
                Bodies[index_body].damping_work = 0. ;
                Bodies[index_body].alid_work = 0. ;
            }
            //cout << "Neighbours" << endl ;
        }

        if (line.substr(0,10)=="NEIGHBOURS")
        {
            int nb_neighbours ;
            Dynamic_file >> nb_neighbours ;
            getline(Dynamic_file, token) ;
            vector<vector<int>> neighbours ;
            int n1, n2, n3 ;
            for (int i(0) ; i < nb_neighbours ; i++)
            {
                Dynamic_file >> n1 >> n2 >> n3 ;
                getline(Dynamic_file, token) ;
                neighbours.push_back(vector<int>( {n1, n2, n3} )) ;
            }
            Bodies[index_body].nb_neighbours=nb_neighbours ;
            Bodies[index_body].neighbours=neighbours ;
            //cout << "Neighbours" << endl ;
        }

        if (line.substr(0,7)=="BORDERS")
        {
            for (int i(0) ; i < Bodies[index_body].nb_borders ; i++)
            {
                getline(Dynamic_file, token) ;
                getline(Dynamic_file, token) ;
                /*
                vector<double> length ;
                vector<int>    node0 ;
                vector<int>    shift0 ;
                vector<int>    node1 ;
                vector<int>    shift1 ;
                vector<int>    node2 ;
                vector<int>    shift2 ;
                vector<int>    node3 ;
                vector<int>    shift3 ;
                */
                int nb_border_nodes ;
                Dynamic_file >> nb_border_nodes ;
                getline(Dynamic_file, token) ;
                //int n0 , n1 , n2 , n3 ;
                //double l , s0 , s1 , s2 , s3 ;
                for (int j(0) ; j < nb_border_nodes ; j++)
                {
                    //Dynamic_file >> l >> n0 >> s0 >> n1 >> s1 >> n2 >> s2 >> n3 >> s3 ;
                    getline(Dynamic_file, token) ;
                    /*
                    length.push_back(l) ;
                    node0.push_back(n0) ;
                    shift0.push_back(s0) ;
                    node1.push_back(n1) ;
                    shift1.push_back(s1) ;
                    node2.push_back(n2) ;
                    shift2.push_back(s2) ;
                    node3.push_back(n3) ;
                    shift3.push_back(s3) ;
                    */
                }
                //Bodies[index_body].borders[i].length = length ;
                //Bodies[index_body].Initialize_segments(Xmin_period , Xmax_period) ;
            }
            //cout << "Borders" << endl ;
        }

        if (line.substr(0,8)=="CONTACTS")
        {
            int nb_contact_elements ;
            Dynamic_file >> nb_contact_elements ;
            getline(Dynamic_file, token) ;
            vector<Contact_element> contact_elements(nb_contact_elements) ;
            for (int i(0) ; i < nb_contact_elements ; i++)
            {
                vector<double> internal ;
                int bos, ns, bns, bm, bom, sm, sem, n0, s0, n1, s1, n2, s2, n3, s3, ni ;
                double gn, gt, xs, xn, yn, xt, yt, l, fx, fy, N0, N1, N2, N3, vi, m ;
                Contact_element new_element ;
                Dynamic_file >> bos >> bns >> ns >> bm >> bom >> sm >> sem >> gn >> gt >> xs ;
                Dynamic_file >> xn >> yn >> xt >> yt ;
                Dynamic_file >> l >> fx >> fy ;
                Dynamic_file >> n0 >> N0 >> s0 ;
                Dynamic_file >> n1 >> N1 >> s1 ;
                Dynamic_file >> n2 >> N2 >> s2 ;
                Dynamic_file >> n3 >> N3 >> s3 ;
                Dynamic_file >> m >> ni ;
                if (ni==0)
                    internal = {0} ;
                else if (ni>0)
                {
                    for (int j(0) ; j<ni ; j++)
                    {
                        Dynamic_file >> vi ;
                        internal.push_back(vi) ;
                        //internal[j] = vi ;
                    }
                }
                new_element.borderS = bos ;
                new_element.border_nodeS = bns ;
                new_element.nodeS = ns ;
                new_element.bodyM = bm ;
                new_element.borderM = bom ;
                new_element.shiftM = sm ;
                new_element.segmentM = sem ;
                new_element.gapn = gn ;
                new_element.gapt = gt ;
                new_element.xsi = xs ;
                new_element.xnorm = xn ;
                new_element.ynorm = yn ;
                new_element.xtan = xt ;
                new_element.ytan = yt ;
                new_element.length = l ;
                new_element.fx = fx ;
                new_element.fy = fy ;
                new_element.nodeM0 = n0 ;
                new_element.shapeM0 = N0 ;
                new_element.shiftM0 = s0 ;
                new_element.nodeM0 = n1 ;
                new_element.shapeM0 = N1 ;
                new_element.shiftM0 = s1 ;
                new_element.nodeM0 = n2 ;
                new_element.shapeM0 = N2 ;
                new_element.shiftM0 = s2 ;
                new_element.nodeM0 = n3 ;
                new_element.shapeM0 = N3 ;
                new_element.shiftM0 = s3 ;
                new_element.effective_mass = m ;
                new_element.nb_internal = ni ;
                new_element.internal = internal ;
                contact_elements[i] = new_element ;
            }
            Bodies[index_body].contact_elements = contact_elements ;
            Bodies[index_body].nb_contact_elements = nb_contact_elements ;
            //cout << "Contacts" << endl ;
        }

        if (line.substr(0,17)=="CONTACT_PRESSURES")
        {
            for (int i(0) ; i < Bodies[index_body].nb_borders ; i++)
            {
                getline(Dynamic_file, token) ;
                getline(Dynamic_file, token) ;
                vector<double> f1 ;
                vector<double> f2 ;
                int nb_border_nodes ;
                Dynamic_file >> nb_border_nodes;
                getline(Dynamic_file, token) ;
                double ff1, ff2 ;
                for (int j(0) ; j < nb_border_nodes ; j++)
                {
                    Dynamic_file >> ff1 >> ff2 ;
                    getline(Dynamic_file, token) ;
                    f1.push_back(ff1) ;
                    f2.push_back(ff2) ;
                }
                Bodies[index_body].borders[i].x_contact_pressure = f1 ;
                Bodies[index_body].borders[i].y_contact_pressure = f2 ;
            }
            //cout << "Contact pressures" << endl ;
        }

        if (line.substr(0,12)=="BC_PRESSURES")
        {
            for (int i(0) ; i < Bodies[index_body].nb_borders ; i++)
            {
                getline(Dynamic_file, token) ;
                getline(Dynamic_file, token) ;
                vector<double> f1 ;
                vector<double> f2 ;
                int nb_border_nodes ;
                Dynamic_file >> nb_border_nodes ;
                getline(Dynamic_file, token) ;
                double ff1, ff2 ;
                for (int j(0) ; j < nb_border_nodes ; j++)
                {
                    Dynamic_file >> ff1 >> ff2 ;
                    getline(Dynamic_file, token) ;
                    f1.push_back(ff1) ;
                    f2.push_back(ff2) ;
                }
                Bodies[index_body].borders[i].x_bc_pressure = f1 ;
                Bodies[index_body].borders[i].y_bc_pressure = f2 ;
            }
            //cout << "BC pressures" << endl ;
        }
    }
    Dynamic_file.close () ;
}



//********************************************//
//** WRITE DYNAMIC ***************************//
//********************************************//

void Write_dynamic(
    string& Simulation_name,
    int Number_save,
    int Number_print,
    double& Time,
    int& Number_iteration,
    double& Deltat,
    string& Solver,
    double& Next_save,
    double& Next_print,
    double& Next_contact_update,
    int Nb_bodies,
    vector<Body>& Bodies,
    vector<int>& flags )
{
    if (flags[7]==0)
        cout << "Writing dynamic data of file " << Number_save << endl ;
    string filename ;

    stringstream sfilename ;
    sfilename << Number_save;
    if      (Number_save<10)
        filename="DYNAMIC_0000"+sfilename.str()+".asc" ;
    else if (Number_save<100)
        filename="DYNAMIC_000"+sfilename.str()+".asc" ;
    else if (Number_save<1000)
        filename="DYNAMIC_00"+sfilename.str()+".asc" ;
    else if (Number_save<10000)
        filename="DYNAMIC_0"+sfilename.str()+".asc" ;
    else if (Number_save<100000)
        filename="DYNAMIC_"+sfilename.str()+".asc" ;
    ofstream Dynamic_file (filename) ;
    Dynamic_file << endl ;
    Dynamic_file << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl ;
    Dynamic_file << "%%%%%%%          GENERAL DATA          %%%%%%%" << endl ;
    Dynamic_file << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl ;
    Dynamic_file << endl ;

    Dynamic_file << "SIMULATION_NAME" << endl ;
    Dynamic_file << Simulation_name << endl ;
    Dynamic_file << endl ;

    Dynamic_file << "SOLVER" << endl ;
    Dynamic_file << Solver << endl ;
    Dynamic_file << setprecision (16) << Number_save << ' ' << Number_iteration << ' ' << Time << ' ' << Number_print << endl ;
    Dynamic_file << setprecision (10) << Deltat << endl ;
    Dynamic_file << Next_save << ' ' << Next_print << ' ' << Next_contact_update << endl ;
    Dynamic_file << endl ;

    if (flags[1]==1)
        Dynamic_file << "MONITOR_ENERGY" << endl ;
    if (flags[2]==1)
        Dynamic_file << "MONITOR_BOUNDARIES" << endl ;
    if (flags[4]==1)
        Dynamic_file << "UPDATE_MASS_MATRIX" << endl ;
    if (flags[5]==1)
        Dynamic_file << "UPDATE_STIFFNESS_MATRIX" << endl ;
    if (flags[6]==1)
        Dynamic_file << "UPDATE_DAMPING_MATRIX" << endl ;
    if (flags[7]==1)
        Dynamic_file << "NO_LOG" << endl ;
    if (flags[8]==1)
        Dynamic_file << "NO_MONITORING" << endl ;
    if (flags[9]==1)
        Dynamic_file << "NO_SELF_CONTACT" << endl ;

    Dynamic_file << endl ;

    for (int i(0) ; i<Nb_bodies ; i++)
    {
        Dynamic_file << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl ;
        Dynamic_file << "%%%%%%%             BODY " << i << "             %%%%%%%" << endl ;
        Dynamic_file << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl ;
        Dynamic_file << endl ;

        if (Bodies[i].type=="deformable")
        {
            Dynamic_file << "DEFORMABLE" << endl ;
            Dynamic_file << Bodies[i].status << endl ;
            Dynamic_file << i << endl ;
            Dynamic_file << endl ;

            Dynamic_file << "KINEMATICS" << endl ;
            Dynamic_file << Bodies[i].nb_nodes << endl ;
            for (int j(0) ; j<Bodies[i].nb_nodes ; j++)
            {
                Dynamic_file << j << ' '
                             << Bodies[i].nodes[j].x_current << ' '
                             << Bodies[i].nodes[j].y_current << ' '
                             << Bodies[i].nodes[j].x_displacement << ' '
                             << Bodies[i].nodes[j].y_displacement << ' '
                             << Bodies[i].nodes[j].x_velocity << ' '
                             << Bodies[i].nodes[j].y_velocity << ' '
                             << Bodies[i].nodes[j].x_acceleration << ' '
                             << Bodies[i].nodes[j].y_acceleration << ' '
                             << Bodies[i].nodes[j].x_displacement_parameter << ' '
                             << Bodies[i].nodes[j].y_displacement_parameter << ' '
                             << Bodies[i].nodes[j].x_velocity_parameter << ' '
                             << Bodies[i].nodes[j].y_velocity_parameter << ' '
                             << Bodies[i].nodes[j].x_acceleration_parameter << ' '
                             << Bodies[i].nodes[j].y_acceleration_parameter << endl ;
            }
            Dynamic_file << endl ;

            Dynamic_file << "FORCES" << endl ;
            for (int j(0) ; j<Bodies[i].nb_nodes ; j++)
            {
                Dynamic_file << j << ' '
                             << Bodies[i].nodes[j].x_internal_force << ' '
                             << Bodies[i].nodes[j].y_internal_force << ' '
                             << Bodies[i].nodes[j].x_contact_force << ' '
                             << Bodies[i].nodes[j].y_contact_force << ' '
                             << Bodies[i].nodes[j].x_self_contact_force << ' '
                             << Bodies[i].nodes[j].y_self_contact_force << ' '
                             << Bodies[i].nodes[j].x_body_force << ' '
                             << Bodies[i].nodes[j].y_body_force << ' '
                             << Bodies[i].nodes[j].x_dirichlet_force << ' '
                             << Bodies[i].nodes[j].y_dirichlet_force << ' '
                             << Bodies[i].nodes[j].x_neumann_force << ' '
                             << Bodies[i].nodes[j].y_neumann_force << ' '
                             << Bodies[i].nodes[j].x_damping_force << ' '
                             << Bodies[i].nodes[j].y_damping_force << ' '
                             << Bodies[i].nodes[j].x_force << ' '
                             << Bodies[i].nodes[j].y_force << endl ;
            }
            Dynamic_file << endl ;
        }
        else if (Bodies[i].type=="rigid")
        {
            Dynamic_file << "RIGID" << endl ;
            Dynamic_file << Bodies[i].status << endl ;
            Dynamic_file << i << endl ;
            Dynamic_file << endl ;

            Dynamic_file << "KINEMATICS" << endl ;
            Dynamic_file << Bodies[i].x_current << ' ' << Bodies[i].y_current << ' ' << Bodies[i].r_current << endl ;
            Dynamic_file << Bodies[i].x_displacement << ' ' << Bodies[i].y_displacement << ' ' << Bodies[i].r_displacement << endl ;
            Dynamic_file << Bodies[i].x_velocity << ' ' << Bodies[i].y_velocity << ' ' << Bodies[i].r_velocity << endl ;
            Dynamic_file << endl ;

            Dynamic_file << "FORCES" << endl ;
            Dynamic_file << Bodies[i].x_contact_force << ' ' << Bodies[i].y_contact_force << ' ' << Bodies[i].r_contact_force << endl ;
            Dynamic_file << Bodies[i].x_body_force << ' ' << Bodies[i].y_body_force << ' ' << Bodies[i].r_body_force << endl ;
            Dynamic_file << Bodies[i].x_dirichlet_force << ' ' << Bodies[i].y_dirichlet_force << ' ' << Bodies[i].r_dirichlet_force << endl ;
            Dynamic_file << Bodies[i].x_neumann_force << ' ' << Bodies[i].y_neumann_force << ' ' << Bodies[i].r_neumann_force << endl ;
            Dynamic_file << Bodies[i].x_damping_force << ' ' << Bodies[i].y_damping_force << ' ' << Bodies[i].r_damping_force << endl ;
            Dynamic_file << Bodies[i].x_force << ' ' << Bodies[i].y_force << ' ' << Bodies[i].r_force << endl ;
            Dynamic_file << endl ;
        }

        Dynamic_file << "WORKS" << endl ;
        Dynamic_file << Bodies[i].internal_work << ' '
                     << Bodies[i].contact_work << ' '
                     << Bodies[i].body_work << ' '
                     << Bodies[i].dirichlet_work << ' '
                     << Bodies[i].neumann_work << ' '
                     << Bodies[i].damping_work << ' '
                     << Bodies[i].alid_work << endl ;
        Dynamic_file << endl ;

        Dynamic_file << "NEIGHBOURS" << endl ;
        Dynamic_file << Bodies[i].nb_neighbours << endl ;
        for (int j(0) ; j<Bodies[i].nb_neighbours ; j++)
        {
            Dynamic_file << Bodies[i].neighbours[j][0] << ' '
                         << Bodies[i].neighbours[j][1] << ' '
                         << Bodies[i].neighbours[j][2] << endl ;
        }
        Dynamic_file << endl ;

        Dynamic_file << "BORDERS" << endl ;
        Dynamic_file << endl ;
        for (int j(0) ; j<Bodies[i].nb_borders ; j++)
        {
            Dynamic_file << j << endl ;
            Dynamic_file << Bodies[i].borders[j].number_border_nodes << endl ;
            vector<double> length=Bodies[i].borders[j].length ;
            vector<int>    node0=Bodies[i].borders[j].node0 ;
            vector<int>    shift0=Bodies[i].borders[j].shift0 ;
            vector<int>    node1=Bodies[i].borders[j].node1 ;
            vector<int>    shift1=Bodies[i].borders[j].shift1 ;
            vector<int>    node2=Bodies[i].borders[j].node2 ;
            vector<int>    shift2=Bodies[i].borders[j].shift2 ;
            vector<int>    node3=Bodies[i].borders[j].node3 ;
            vector<int>    shift3=Bodies[i].borders[j].shift3 ;
            for (int k(0) ; k<Bodies[i].borders[j].number_border_nodes ; k++)
            {
                Dynamic_file << length[k] << ' ' ;
                Dynamic_file << node0[k] << ' ' << shift0[k] << ' ' << node1[k] << ' ' << shift1[k] << ' ' ;
                Dynamic_file << node2[k] << ' ' << shift2[k] << ' ' << node3[k] << ' ' << shift3[k] << endl ;
            }
            Dynamic_file << endl ;
        }

        Dynamic_file << "CONTACTS" << endl ;
        Dynamic_file << Bodies[i].nb_contact_elements << endl ;
        for (int k(0) ; k<Bodies[i].nb_contact_elements ; k++)
        {
            int borderS=Bodies[i].contact_elements[k].borderS ;
            int border_nodeS=Bodies[i].contact_elements[k].border_nodeS ;
            int nodeS=Bodies[i].contact_elements[k].nodeS ;
            int bodyM=Bodies[i].contact_elements[k].bodyM ;
            int borderM=Bodies[i].contact_elements[k].borderM ;
            int shiftM=Bodies[i].contact_elements[k].shiftM ;
            int segmentM=Bodies[i].contact_elements[k].segmentM ;
            double gapn=Bodies[i].contact_elements[k].gapn ;
            double gapt=Bodies[i].contact_elements[k].gapt ;
            double xsi=Bodies[i].contact_elements[k].xsi ;
            double xnorm=Bodies[i].contact_elements[k].xnorm ;
            double ynorm=Bodies[i].contact_elements[k].ynorm ;
            double xtan=Bodies[i].contact_elements[k].xtan ;
            double ytan=Bodies[i].contact_elements[k].ytan ;
            double length=Bodies[i].contact_elements[k].length ;
            double fx=Bodies[i].contact_elements[k].fx ;
            double fy=Bodies[i].contact_elements[k].fy ;
            int nodeM0=Bodies[i].contact_elements[k].nodeM0 ;
            double shapeM0=Bodies[i].contact_elements[k].shapeM0 ;
            int shiftM0=Bodies[i].contact_elements[k].shiftM0 ;
            int nodeM1=Bodies[i].contact_elements[k].nodeM1 ;
            double shapeM1=Bodies[i].contact_elements[k].shapeM1 ;
            int shiftM1=Bodies[i].contact_elements[k].shiftM1 ;
            int nodeM2=Bodies[i].contact_elements[k].nodeM2 ;
            double shapeM2=Bodies[i].contact_elements[k].shapeM2 ;
            int shiftM2=Bodies[i].contact_elements[k].shiftM2 ;
            int nodeM3=Bodies[i].contact_elements[k].nodeM3 ;
            double shapeM3=Bodies[i].contact_elements[k].shapeM3 ;
            int shiftM3=Bodies[i].contact_elements[k].shiftM3 ;
            double m=Bodies[i].contact_elements[k].effective_mass ;
            int ni = Bodies[i].contact_elements[k].nb_internal ;
            vector<double> internal = Bodies[i].contact_elements[k].internal ;
            Dynamic_file << borderS << ' ' << border_nodeS << ' ' << nodeS << ' ' ;
            Dynamic_file << bodyM << ' ' << borderM << ' ' << shiftM << ' ' << segmentM << ' ' ;
            Dynamic_file << gapn << ' ' << gapt << ' ' << xsi << ' ' ;
            Dynamic_file << xnorm << ' ' << ynorm << ' ' << xtan << ' ' << ytan << ' ' ;
            Dynamic_file << length << ' ' << fx << ' ' << fy << ' ' ;
            Dynamic_file << nodeM0 << ' ' << shapeM0 << ' ' << shiftM0 << ' ' ;
            Dynamic_file << nodeM1 << ' ' << shapeM1 << ' ' << shiftM1 << ' ' ;
            Dynamic_file << nodeM2 << ' ' << shapeM2 << ' ' << shiftM2 << ' ' ;
            Dynamic_file << nodeM3 << ' ' << shapeM3 << ' ' << shiftM3 << ' ' ;
            Dynamic_file << m << ' ' << ni ;
            if (ni>0)
            {
                for (int kk(0) ; kk<ni ; kk++)
                {
                    Dynamic_file << ' ' << internal[kk] ;
                }
            }
            Dynamic_file << endl ;
        }
        Dynamic_file << endl ;

        Dynamic_file << "CONTACT_PRESSURES" << endl ;
        Dynamic_file << endl ;
        for (int j(0) ; j<Bodies[i].nb_borders ; j++)
        {
            Dynamic_file << j << endl ;
            Dynamic_file << Bodies[i].borders[j].number_border_nodes << endl ;
            vector<double> x_contact_pressure=Bodies[i].borders[j].x_contact_pressure ;
            vector<double> y_contact_pressure=Bodies[i].borders[j].y_contact_pressure ;
            for (int k(0) ; k<Bodies[i].borders[j].number_border_nodes ; k++)
            {
                Dynamic_file << x_contact_pressure[k] << ' ' << y_contact_pressure[k] << endl ;
            }
            Dynamic_file << endl ;
        }
        Dynamic_file << endl ;

        Dynamic_file << "BC_PRESSURES" << endl ;
        Dynamic_file << endl ;
        for (int j(0) ; j<Bodies[i].nb_borders ; j++)
        {
            Dynamic_file << j << endl ;
            Dynamic_file << Bodies[i].borders[j].number_border_nodes << endl ;
            vector<double> x_bc_pressure=Bodies[i].borders[j].x_bc_pressure ;
            vector<double> y_bc_pressure=Bodies[i].borders[j].y_bc_pressure ;
            for (int k(0) ; k<Bodies[i].borders[j].number_border_nodes ; k++)
            {
                Dynamic_file << x_bc_pressure[k] << ' ' << y_bc_pressure[k] << endl ;
            }
            Dynamic_file << endl ;
        }
        Dynamic_file << endl ;
    }

    Dynamic_file.close() ;
}



//********************************************//
//** WRITE MONITORING ************************//
//********************************************//

void Write_monitoring(vector<vector<double>>& monitoring,
                      vector<int>& flags)
{
    if (flags[7]==0)
        cout << "Writing monitored data " << endl ;
    string filename="MONITORING.asc" ;
    vector<double> m ;
    int sizei, sizej ;
    ofstream Monitoring_file(filename, ios::app) ;
    sizei = (int) monitoring.size() ;
    Monitoring_file << setprecision (10) ;
    for (int i(0) ; i<sizei ; i++)
    {
        m = monitoring[i] ;
        sizej = (int) m.size() ;
        for (int j(0) ; j<sizej ; j++)
        {
            Monitoring_file << m[j] << ' ' ;
        }
        Monitoring_file << endl ;
    }
    monitoring.clear() ;
}



//********************************************//
//** WRITE SPYING ****************************//
//********************************************//

void Write_spying(vector<vector<vector<double>>>& spying,
                  int& Nb_spies,
                  vector<Spy>& Spies,
                  vector<int>& flags)
{
    if (flags[7]==0)
        cout << "Writing spied data " << endl ;
    for (int n(0) ; n < Nb_spies ; n++)
    {
        string filename = Spies[n].filename + ".asc" ;
        vector<double> m ;
        int sizei, sizej ;
        ofstream Monitoring_file(filename, ios::app) ;
        sizei = (int) spying[n].size() ;
        Monitoring_file << setprecision (10) ;
        for (int i(0) ; i<sizei ; i++)
        {
            m = spying[n][i] ;
            sizej = (int) m.size() ;
            for (int j(0) ; j<sizej ; j++)
            {
                Monitoring_file << m[j] << ' ' ;
            }
            Monitoring_file << endl ;
        }
        spying[n].clear() ;
    }
}

#endif
