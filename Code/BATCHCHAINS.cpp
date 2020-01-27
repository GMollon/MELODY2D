#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>
using namespace std ;
#include "Tools.h"
#include "Material.h"
#include "Contact_Law.h"
#include "Node.h"
#include "Gauss.h"
#include "Border.h"
#include "Contact_Element.h"
#include "Body.h"
#include "Spy.h"
#include "IO.h"
#include "Solver.h"
#include "Proximity.h"
#include "Graphic.h"
#include "Monitoring.h"

main(int argc, char **argv)
{
    // CREATE VARIABLES //
    string Simulation_name ;
    int Nb_materials ;
    vector<Material> Materials ;
    int Nb_contact_laws ;
    vector<Contact_law> Contact_laws ;
    string Solver ;
    double Tini, Deltat, Tend, Time ;
    double Target_error, Control_parameter, Accepted_ratio ;
    double Save_period, Print_period, Contact_update_period ;
    double Next_save, Next_print, Next_contact_update ;
    int Number_save, Number_print, Number_iteration ;
    double Xmin_period, Xmax_period, Penalty ;
    double Xgravity, Ygravity ;
    int Activate_plot ;
    double Xmin_plot, Xmax_plot, Ymin_plot, Ymax_plot ;
    int Nb_bodies ;
    vector<Body> Bodies ;
    int Nb_monitored ;
    vector<vector<double>> Monitored ;
    int Nb_deactivated ;
    vector<vector<double>> Deactivated ;
    int Nb_spies ;
    vector<Spy> Spies ;
    int Nb_regions = 0 ;
    vector<vector<int>> Regions ;
    vector<int> flags(11) ;
    vector<int> To_Plot(40) ;
    vector<vector<int>> Contacts_Table ;

    // LOAD STATIC DATA //
    Load_static( Simulation_name, Nb_materials, Materials,
                 Nb_contact_laws, Contact_laws, Contacts_Table,
                 Solver, Tini,	Deltat, Tend,
                 Target_error, Control_parameter, Accepted_ratio,
                 Save_period, Print_period, Contact_update_period,
                 Xmin_period, Xmax_period, Penalty, Xgravity, Ygravity,
                 Activate_plot,	Xmin_plot,	Xmax_plot,	Ymin_plot,	Ymax_plot,
                 Nb_monitored, Monitored, Nb_deactivated, Deactivated, Nb_spies, Spies,
                 Nb_regions, Regions, Nb_bodies, Bodies, To_Plot ) ;

    double Typical_pressure = atof(argv[1]) ;
    double Size_ratio = atof(argv[2]) ;
    int istart = atoi(argv[3]) ;
    int interval = atoi(argv[4]) ;
    int iend = atoi(argv[5]) ;

    for (int i(istart) ; i<=iend ; i+=interval)
    {
        Number_save = i ;
        Load_dynamic( Number_save, Number_print,
                      Time, Number_iteration, Deltat, Solver,
                      Next_save,	Next_print, Next_contact_update,
                      Xmin_period, Xmax_period,
                      Bodies, flags ) ;
        Number_print = i ;
        Write_chains(Nb_bodies, Bodies, Number_iteration, Number_save, Number_print, Time, Xmin_period, Xmax_period, Typical_pressure, Size_ratio) ;
    }
}
