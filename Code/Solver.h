#include "Proximity.h"
#ifndef DEF_SOLVER
#define DEF_SOLVER


//********************************************//
//** EULER STEP ******************************//
//********************************************//

void Euler_step(
    int Nb_bodies,
    vector<Body>& Bodies,
    int Nb_regions,
    vector<vector<int>> Regions,
    int Nb_materials,
    vector<Material> Materials,
    int Nb_contact_laws,
    vector<Contact_law> Contact_laws,
    vector<vector<int>>& Contacts_Table,
    double Tend,
    double Xmin_period,
    double Xmax_period,
    double Penalty,
    double Xgravity,
    double Ygravity,
    double& Time,
    double& Deltat,
    int& Number_iteration,
    vector<int>& flags)
{
    if (flags[7]==0)                                                                                    // SEQUENTIAL //
    {
        cout << endl ;                                                                                  // SEQUENTIAL //
        cout << Number_iteration << " ; t " << Time ;                                                   // SEQUENTIAL //
    }
    #pragma omp parallel
    {
        #pragma omp for schedule(dynamic)
        for (int i=0 ; i<Nb_bodies ; i++)
        {
            if (Bodies[i].status == "inactive")
                continue ;
            Bodies[i].Initialize_contact_forces() ;
            //Bodies[i].Update_bc(Time) ;
            Bodies[i].Update_borders(Xmin_period, Xmax_period) ;
            Bodies[i].Update_contacts(Bodies, Xmin_period, Xmax_period) ;
        }
        #pragma omp for schedule(dynamic)
        for (int i=0 ; i<Nb_bodies ; i++)
        {
            if (Bodies[i].status == "inactive")
                continue ;
            Bodies[i].Update_bc_forces() ;
            Bodies[i].Update_body_forces(Xgravity, Ygravity) ;
            Bodies[i].Update_damping_forces() ;
            Bodies[i].Update_alid_forces(Time) ;
            Bodies[i].Update_contact_forces(Deltat, Bodies, Nb_contact_laws, Contact_laws, Contacts_Table, Xmin_period, Xmax_period) ;
        }
    }

    for (int i=0 ; i<Nb_bodies ; i++)
        Bodies[i].Send_contact_forces(Bodies) ;

    #pragma omp parallel
    {
        #pragma omp for schedule(dynamic)
        for (int i=0 ; i<Nb_regions ; i++)
        {
            if (Bodies[Regions[i][0]].status == "inactive")
                continue ;
            Bodies[Regions[i][0]].Update_internal_forces(Regions[i][1], Nb_materials, Materials) ;
        }
        #pragma omp for schedule(dynamic)
        for (int i=0 ; i<Nb_bodies ; i++)
        {
            if (Bodies[i].status == "inactive")
                continue ;
            Bodies[i].Sum_up_forces() ;
            Bodies[i].Apply_Newton() ;
            Bodies[i].Update_bc(Time+Deltat) ;
            Bodies[i].Apply_Euler(Deltat) ;
            Bodies[i].Update_kinematics() ;
            Bodies[i].Update_current_positions() ;
        }
    }
    Time += Deltat ;                                                                            // SEQUENTIAL //
    Update_contact_pressures( Nb_bodies, Bodies ) ;                                             // SEQUENTIAL //
}



//********************************************//
//** ADAPTIVE EULER STEP *********************//
//********************************************//

void Adaptive_Euler_step(
    int Nb_bodies,
    vector<Body>& Bodies,
    int Nb_regions,
    vector<vector<int>> Regions,
    int Nb_materials,
    vector<Material> Materials,
    int Nb_contact_laws,
    vector<Contact_law> Contact_laws,
    vector<vector<int>>& Contacts_Table,
    double Tend,
    double Xmin_period,
    double Xmax_period,
    double Penalty,
    double Xgravity,
    double Ygravity,
    double& Time,
    double& Deltat,
    int& Number_iteration,
    double& Next_contact_update,
    double& Contact_update_period,
    double Target_error,
    double Control_parameter,
    double Accepted_ratio,
    int& neval,
    double& max_error,
    double& mean_error,
    vector<int>& flags)
{
    if (flags[7]==0)                                                                                    // SEQUENTIAL //
    {
        cout << endl ;                                                                                  // SEQUENTIAL //
        cout << Number_iteration << " ; t " << Time ;                                                   // SEQUENTIAL //
    }
    #pragma omp parallel
    {
        #pragma omp for schedule(dynamic)
        for (int i=0 ; i<Nb_bodies ; i++)
            Bodies[i].Store() ;
    }
    int bodyformax, nodeformax ;                                                                        // SEQUENTIAL //
    int flag_success = 0 ;                                                                              // SEQUENTIAL //
    neval = 0 ;                                                                                         // SEQUENTIAL //
    while (flag_success == 0)
    {
        neval++ ;                                                                                       // SEQUENTIAL //
        #pragma omp parallel
        {
            #pragma omp for schedule(dynamic)
            for (int i=0 ; i<Nb_bodies ; i++)
            {
                if (Bodies[i].status == "inactive")
                    continue ;
                Bodies[i].Initialize_contact_forces() ;
                //Bodies[i].Update_bc(Time) ;
                Bodies[i].Update_borders(Xmin_period, Xmax_period) ;
                Bodies[i].Update_contacts(Bodies, Xmin_period, Xmax_period) ;
            }
            #pragma omp for schedule(dynamic)
            for (int i=0 ; i<Nb_bodies ; i++)
            {
                if (Bodies[i].status == "inactive")
                    continue ;
                Bodies[i].Update_bc_forces() ;
                Bodies[i].Update_body_forces(Xgravity, Ygravity) ;
                Bodies[i].Update_damping_forces() ;
                Bodies[i].Update_alid_forces(Time) ;
                Bodies[i].Update_contact_forces(Deltat*0.5, Bodies, Nb_contact_laws, Contact_laws, Contacts_Table, Xmin_period, Xmax_period) ;
            }
        }

        for (int i=0 ; i<Nb_bodies ; i++)
            Bodies[i].Send_contact_forces(Bodies) ;

        #pragma omp parallel
        {
            #pragma omp for schedule(dynamic)
            for (int i=0 ; i<Nb_regions ; i++)
            {
                if (Bodies[Regions[i][0]].status == "inactive")
                    continue ;
                Bodies[Regions[i][0]].Update_internal_forces(Regions[i][1], Nb_materials, Materials) ;
            }
            #pragma omp for schedule(dynamic)
            for (int i=0 ; i<Nb_bodies ; i++)
            {
                if (Bodies[i].status == "inactive")
                    continue ;
                Bodies[i].Sum_up_forces() ;
                Bodies[i].Apply_Newton() ;
                Bodies[i].Apply_Euler_temporary(Deltat) ;
                Bodies[i].Update_kinematics_temporary() ;
                Bodies[i].Apply_Euler(Deltat*0.5) ;
                Bodies[i].Update_kinematics() ;
                Bodies[i].Update_current_positions() ;
            }
            #pragma omp for schedule(dynamic)
            for (int i=0 ; i<Nb_bodies ; i++)
            {
                if (Bodies[i].status == "inactive")
                    continue ;
                Bodies[i].Initialize_contact_forces() ;
                Bodies[i].Update_bc(Time+Deltat*0.5) ;
                Bodies[i].Update_borders(Xmin_period, Xmax_period) ;
                Bodies[i].Update_contacts(Bodies, Xmin_period, Xmax_period) ;
            }
            #pragma omp for schedule(dynamic)
            for (int i=0 ; i<Nb_bodies ; i++)
            {
                if (Bodies[i].status == "inactive")
                    continue ;
                Bodies[i].Update_bc_forces() ;
                Bodies[i].Update_body_forces(Xgravity, Ygravity) ;
                Bodies[i].Update_damping_forces() ;
                Bodies[i].Update_alid_forces(Time+Deltat*0.5) ;
                Bodies[i].Update_contact_forces(Deltat*0.5, Bodies, Nb_contact_laws, Contact_laws, Contacts_Table, Xmin_period, Xmax_period) ;
            }
        }

        for (int i=0 ; i<Nb_bodies ; i++)
            Bodies[i].Send_contact_forces(Bodies) ;

        #pragma omp parallel
        {
            #pragma omp for schedule(dynamic)
            for (int i=0 ; i<Nb_regions ; i++)
            {
                if (Bodies[Regions[i][0]].status == "inactive")
                    continue ;
                Bodies[Regions[i][0]].Update_internal_forces(Regions[i][1], Nb_materials, Materials) ;
            }
            #pragma omp for schedule(dynamic)
            for (int i=0 ; i<Nb_bodies ; i++)
            {
                if (Bodies[i].status == "inactive")
                    continue ;
                Bodies[i].Sum_up_forces() ;
                Bodies[i].Apply_Newton() ;
                Bodies[i].Update_bc(Time+Deltat) ;
                Bodies[i].Apply_Euler(Deltat*0.5) ;
                Bodies[i].Update_kinematics() ;
                Bodies[i].Update_current_positions() ;
                Bodies[i].Compute_error() ;
            }
        }
        mean_error = 0. ;                                                                               // SEQUENTIAL //
        max_error = 0. ;                                                                                // SEQUENTIAL //
        int nnodes = 0 ;                                                                                // SEQUENTIAL //
        for (int i=0 ; i<Nb_bodies ; i++)                                                               // SEQUENTIAL //
        {
            if (Bodies[i].status == "inactive")
                continue ;                                              // SEQUENTIAL //
            mean_error += Bodies[i].total_error ;                                                       // SEQUENTIAL //
            nnodes += Bodies[i].nb_nodes ;                                                              // SEQUENTIAL //
            if (Bodies[i].max_error > max_error)                                                        // SEQUENTIAL //
            {
                max_error = Bodies[i].max_error ;                                                       // SEQUENTIAL //
                bodyformax = i ;                                                                        // SEQUENTIAL //
                nodeformax = Bodies[i].node_for_max_error ;                                             // SEQUENTIAL //
            }
        }
        mean_error /= nnodes ;                                                                          // SEQUENTIAL //

        if (flags[7]==0)
            cout << " ; dt " << Deltat << " ; err " << max_error << " -" << bodyformax << "-" << nodeformax << "- (" << mean_error << ")" << endl ;
        double error = max_error ;                                                                      // SEQUENTIAL //
        if (error > Target_error * Accepted_ratio)                                                      // SEQUENTIAL //
        {
            flag_success = 0 ;                                                                          // SEQUENTIAL //
            #pragma omp parallel
            {
                #pragma omp for schedule(dynamic)
                for (int i=0 ; i<Nb_bodies ; i++)
                    Bodies[i].Restore() ;
            }
            double correction1 = error / Target_error ;                                                 // SEQUENTIAL //
            double correction2 = 1. / Accepted_ratio ;                                                  // SEQUENTIAL //
            if (correction1>correction2)
                Deltat = Deltat / pow(correction1, Control_parameter) ;     // SEQUENTIAL //
            else
                Deltat = Deltat / pow(correction2, Control_parameter) ;     // SEQUENTIAL //
            if (Deltat>Contact_update_period)
                Deltat = Contact_update_period ;                        // SEQUENTIAL //
            if (flags[7]==0)
                cout << "                 Failed" ;                                        // SEQUENTIAL //
        }
        else                                                                                            // SEQUENTIAL //
        {
            flag_success = 1 ;                                                                          // SEQUENTIAL //
            Time += Deltat ;                                                                            // SEQUENTIAL //
            double correction1 = error / Target_error ;                                                 // SEQUENTIAL //
            double correction2 = 1. / Accepted_ratio ;                                                  // SEQUENTIAL //
            if (correction1>correction2)
                Deltat = Deltat / pow(correction1, Control_parameter) ;     // SEQUENTIAL //
            else
                Deltat = Deltat / pow(correction2, Control_parameter) ;     // SEQUENTIAL //
            if (Deltat>Contact_update_period)
                Deltat = Contact_update_period ;                        // SEQUENTIAL //
            Update_contact_pressures( Nb_bodies, Bodies ) ;                                             // SEQUENTIAL //
        }
    }
}

#endif
