#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>
#include <omp.h>
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
    cout << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl ;
    cout << "%%                                                       %%" << endl ;
    cout << "%%                       MELODY 2D                       %%" << endl ;
    cout << "%%            Multibody ELement-free Open code           %%" << endl ;
    cout << "%%              for DYnamic simulation in 2D             %%" << endl ;
    cout << "%%                                                       %%" << endl ;
    cout << "%%                     Main Program                      %%" << endl ;
    cout << "%%        Version 3.81 ; 12th of September 2018          %%" << endl ;
    cout << "%%                                                       %%" << endl ;
    cout << "%%                Author: Guilhem Mollon                 %%" << endl ;
    cout << "%%                                                       %%" << endl ;
    cout << "%%            Developed at LaMCoS, INSA Lyon,            %%" << endl ;
    cout << "%%                     Lyon, France                      %%" << endl ;
    cout << "%%                      Year 2015                        %%" << endl ;
    cout << "%%                                                       %%" << endl ;
    cout << "%%            Please read enclosed .txt file             %%" << endl ;
    cout << "%%                                                       %%" << endl ;
    cout << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl ;

    // CREATE VARIABLES //
    string Simulation_name ;
    int Nb_materials ;
    vector<Material> Materials ;
    int Nb_contact_laws ;
    vector<Contact_law> Contact_laws ;
    string Solver ;
    double Tini, Deltat, Tend, Time ;
    double Target_error, Inv_Target_error, Control_parameter, Accepted_ratio ;
    double Max_mass_scaling, Control_parameter_mass_scaling, Error_factor_mass_scaling, Decrease_factor_mass_scaling ;
    double Save_period, Print_period, Contact_update_period ;
    double Next_save, Next_print, Next_contact_update ;
    int Number_save(0), Number_print, Number_iteration ;
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
    vector<vector<double>> monitoring ;
    int neval ;
    double max_error, mean_error ;
    double total_mass, max_mass ;
    int Nb_regions = 0 ;
    vector<vector<int>> Regions ;
    vector<int> flags(9) ;
    int flag_failure = 0 ;
    vector<int> To_Plot(34) ;
    vector<vector<int>> Contacts_Table ;

    // LOAD STATIC DATA //
    Load_static( Simulation_name, Nb_materials, Materials,
                 Nb_contact_laws, Contact_laws, Contacts_Table,
                 Solver, Tini,	Deltat, Tend,
                 Target_error, Inv_Target_error, Control_parameter, Accepted_ratio,
                 Max_mass_scaling, Control_parameter_mass_scaling, Error_factor_mass_scaling, Decrease_factor_mass_scaling,
                 Save_period, Print_period, Contact_update_period,
                 Xmin_period, Xmax_period, Penalty, Xgravity, Ygravity,
                 Activate_plot,	Xmin_plot,	Xmax_plot,	Ymin_plot,	Ymax_plot,
                 Nb_monitored, Monitored, Nb_deactivated, Deactivated, Nb_spies, Spies,
                 Nb_regions, Regions, Nb_bodies, Bodies, To_Plot ) ;

    Number_save = atoi(argv[1]) ;

    // LOAD DYNAMIC DATA //
    Load_dynamic( Number_save, Number_print,
                  Time, Number_iteration, Deltat, Solver,
                  Next_save,	Next_print, Next_contact_update,
                  Xmin_period, Xmax_period,
                  Bodies, flags ) ;

    // INITIALIZATION //
    if (flags[7]==0)
        cout << "Initializing" << endl ;
    for (int i=0 ; i<Nb_bodies ; i++)
        Bodies[i].Update_borders(Xmin_period, Xmax_period) ;
    Update_proximity(Nb_bodies, Bodies, Xmin_period, Xmax_period) ;
    if (flags[3]==1)
        Initialize_CZM( Nb_bodies, Bodies, Nb_contact_laws, Contact_laws, flags, Xmin_period, Xmax_period ) ;
    cout << "Updating Material Properties" << endl ;
    #pragma omp parallel
    {
        #pragma omp for schedule(dynamic)
        for (int i=0 ; i<Nb_bodies ; i++)
            Bodies[i].Update_material( Nb_materials, Materials, flags ) ;
    }
    for (int i=0 ; i<Nb_spies ; i++)
        Spies[i].next_time += Time ;
    vector<vector<vector<double>>> spying(Nb_spies) ;

    // PRINT INITIAL STATE //
    if ( Number_save == 0 )
        Write_graphic(Nb_bodies, Bodies, Number_iteration, Number_save, Number_print, Time, Xmin_period, Xmax_period, Nb_materials, Materials, To_Plot) ;

    // MAIN PROGRAM //
    if (flags[7]==0)
        cout << "Running main program" << endl ;

    while ( Time < Tend )
    {
        Number_iteration++ ;                                                         // SEQUENTIAL //

        //********************************************//
        //** SOLVER **********************************//
        //********************************************//
        if ( Solver == "Euler" )
        {
            Euler_step(Nb_bodies, Bodies, Nb_materials,
                       Materials, Nb_contact_laws, Contact_laws, Contacts_Table,
                       Tend, Xmin_period, Xmax_period, Penalty,
                       Xgravity, Ygravity,
                       Time, Deltat, Number_iteration) ;
        }
        else if ( Solver == "Adaptive_Euler" || Solver == "Mass_Scaling" || Solver == "Adaptive_Mass_Scaling" || Solver == "Euler_2")
        {
            Solver_step(Nb_bodies, Bodies, Nb_regions, Regions,
                                Nb_materials, Materials, Nb_contact_laws, Contact_laws, Contacts_Table,
                                Tend, Xmin_period, Xmax_period, Penalty,
                                Xgravity, Ygravity,
                                Time, Deltat, Number_iteration,
                                Next_contact_update, Contact_update_period,
                                Target_error, Inv_Target_error, Control_parameter, Accepted_ratio,
                                Max_mass_scaling, Control_parameter_mass_scaling, Error_factor_mass_scaling, Decrease_factor_mass_scaling,
                                neval, max_error, mean_error, flags, total_mass, max_mass, Solver) ;
            //cout << "   SOLVER : " << Solver << " " << endl ;
        }


        //********************************************//
        //** MONITORING ******************************//
        //********************************************//
        if (flags[8]==0)
            Monitoring(Nb_bodies, Bodies,
                       Time, Deltat, Number_iteration,
                       neval, max_error, mean_error, flags,
                       monitoring, Nb_monitored, Monitored) ;                           // SEQUENTIAL //

        //********************************************//
        //** SPYING **********************************//
        //********************************************//
        Spying(Nb_bodies, Bodies,
               Time, Deltat, Number_iteration,
               spying, Nb_spies, Spies) ;                                               // SEQUENTIAL //
        //cout << "45" << endl ;

        //********************************************//
        //** CONTACT UPDATE **************************//
        //********************************************//
        if ( Time > Next_contact_update-1.e-6*Deltat )                                  // SEQUENTIAL //
        {
            Next_contact_update = Next_contact_update + Contact_update_period ;
            Update_proximity(Nb_bodies, Bodies, Xmin_period, Xmax_period) ;
            Update_status(Nb_bodies, Bodies, Nb_deactivated, Deactivated) ;             // SEQUENTIAL //
        }
        //cout << "46" << endl ;

        //********************************************//
        //** WRITE GRAPHIC ***************************//
        //********************************************//
        if ( Time > Next_print-1.e-6*Deltat )                                           // SEQUENTIAL //
        {
            Number_print = Number_print + 1 ;                                           // SEQUENTIAL //
            Next_print = Next_print + Print_period ;                                    // SEQUENTIAL //
            Write_graphic(Nb_bodies, Bodies, Number_iteration, Number_save, Number_print, Time, Xmin_period, Xmax_period, Nb_materials, Materials, To_Plot) ;
        }
        //cout << "47" << endl ;

        //********************************************//
        //** SAVE ************************************//
        //********************************************//
        if ( Time > Next_save-1.e-6*Deltat )                                            // SEQUENTIAL //
        {
            Number_save = Number_save + 1 ;                                             // SEQUENTIAL //
            Next_save = Next_save + Save_period ;                                       // SEQUENTIAL //
            Write_dynamic( Simulation_name,  Number_save, Number_print,
                           Time, Number_iteration, Deltat, Solver,
                           Next_save,	Next_print, Next_contact_update,
                           Nb_bodies, Bodies, flags ) ;                                 // SEQUENTIAL //
            if (flags[8]==0)
                Write_monitoring( monitoring, flags ) ;                                 // SEQUENTIAL //
            Write_spying( spying, Nb_spies, Spies, flags ) ;                            // SEQUENTIAL //
        }
        //cout << "48" << endl ;


        for (int i=0 ; i<Nb_bodies ; i++) // Check failure
        {
            if (std::isnan(Bodies[i].total_error))
            {
                flag_failure = 1 ;
                break ;
            }
        }
        //cout << "49" << endl ;
        if (flag_failure == 1)
            break ;
    }
}
