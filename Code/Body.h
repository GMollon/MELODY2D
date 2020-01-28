#ifndef DEF_BODY
#define DEF_BODY
#include "Model.h"

// Header of the class

class Body
{
public :

    // Static attributes
    int index ;
    string material_name ;
    int material_index ;
    string periodicity ;
    string type ;
    string integration ;
    double nodal_distance ;
    double detection_distance ;
    double contact_distance ;
    int nb_nodes ;
    int nb_borders ;
    int nb_border_nodes ;
    int nb_gauss ;
    int nb_cells ;
    int nb_matrix ;
    int nb_regions ;

    // Dynamic attributes
    int nb_neighbours ;
    int nb_contact_elements ;
    vector<Node> nodes ;
    vector<Node> stored_nodes ;
    vector<Border> borders ;
    vector<Border> stored_borders ;
    vector<Contact_element> contact_elements ;
    vector<Contact_element> stored_contact_elements ;
    vector<Gauss> gpoints ;
    vector<vector<int>> triangulation ;
    vector<vector<int>> matrix_coordinates ;
    vector<double> stiffness ;
    vector<double> damping ;
    vector<vector<int>> neighbours ;
    vector<vector<int>> region_gpoints ;
    vector<double> Elastic_energy ;
    string status ;

    // Rigid attributes
    double mass ;
    double inertia ;
    double inverse_mass ;
    double inverse_inertia ;
    double x_initial ;
    double y_initial ;
    double r_initial ;
    double x_current ;
    double y_current ;
    double r_current ;
    double stored_x_current ;
    double stored_y_current ;
    double stored_r_current ;
    double x_displacement ;
    double y_displacement ;
    double r_displacement ;
    double x_displacement_temporary ;
    double y_displacement_temporary ;
    double r_displacement_temporary ;
    double stored_x_displacement ;
    double stored_y_displacement ;
    double stored_r_displacement ;
    double x_velocity ;
    double y_velocity ;
    double r_velocity ;
    double x_velocity_temporary ;
    double y_velocity_temporary ;
    double r_velocity_temporary ;
    double stored_x_velocity ;
    double stored_y_velocity ;
    double stored_r_velocity ;
    double x_acceleration ;
    double y_acceleration ;
    double r_acceleration ;
    double stored_x_acceleration ;
    double stored_y_acceleration ;
    double stored_r_acceleration ;
    double x_contact_force ;
    double y_contact_force ;
    double r_contact_force ;
    //vector<double> x_master_contact_forces ;
    //vector<double> y_master_contact_forces ;
    //vector<double> r_master_contact_forces ;
    double x_body_force ;
    double y_body_force ;
    double r_body_force ;
    double x_dirichlet_force ;
    double y_dirichlet_force ;
    double r_dirichlet_force ;
    double x_neumann_force ;
    double y_neumann_force ;
    double r_neumann_force ;
    double x_damping_force ;
    double y_damping_force ;
    double r_damping_force ;
    double x_force ;
    double y_force ;
    double r_force ;
    double unbalanced ;
    double total_error ;
    double max_error ;
    int node_for_max_error ;
    double damage ;
    vector<double> drivendof ;
    double internal_work ;
    double contact_work ;
    double body_work ;
    double dirichlet_work ;
    double neumann_work ;
    double damping_work ;
    double alid_work ;
    double kinetic_energy ;
    double deformation_energy ;
    double stored_internal_work ;
    double stored_contact_work ;
    double stored_body_work ;
    double stored_dirichlet_work ;
    double stored_neumann_work ;
    double stored_damping_work ;
    double stored_alid_work ;
    int flag_alid ;
    double alid_xmin ;
    double alid_xmax ;
    double alid_ymin ;
    double alid_ymax ;
    double alid_range ;
    double alid_alphain ;
    double alid_alphaout ;
    double alid_betain ;
    double alid_betaout ;
    double alid_exponent ;
    vector<double> xalid_instants ;
    vector<double> xalid_values ;
    vector<double> yalid_instants ;
    vector<double> yalid_values ;
    double density ;
    vector<vector<double>> contact_forces_to_send ;
    vector<vector<int>> contact_forces_to_send_to ;

    // Mass scaling attributes
    double delta_factor_mass_scaling ;
    double stored_delta_factor_mass_scaling ;
    double factor_mass_scaling ;
    double max_factor_mass_scaling ;
    double mass_mass_scaling;

    //Thermal attributes
    double heat_capacity ;
    double temperature ;

    // Constructor and Destructor
    Body(int i, string m, string p, string t) ;
    ~Body() ;

    // Accessors

    // Modifiers

    // Methods
    void Update_bc(double t) ;
    void Initialize_segments(double xmin, double xmax) ;
    void Update_borders(double xmin, double xmax) ;
    void Update_bc_forces() ;
    void Update_body_forces(double xg, double yg) ;
    void Update_internal_forces(int region, int Nb_materials, vector<Material>& Materials) ;
    void Update_alid_forces(double time) ;
    void Update_damping_forces() ;
    void Initialize_contact_forces() ;
    void Update_contacts(vector<Body>& b, double xmin, double xmax) ;
    void Update_contact_forces(double dt, vector<Body>& b, int Nb_contact_laws, vector<Contact_law>& Contact_laws, vector<vector<int>>& Contacts_Table, double xmin, double xmax) ;
    void Send_contact_forces(vector<Body>& b) ;
    void Sum_up_forces() ;
    void Apply_Newton() ;
    void Apply_Euler(double dt) ;
    void Apply_Euler_temporary(double dt) ;
    void Update_kinematics() ;
    void Update_kinematics_temporary() ;
    void Update_current_positions() ;
    void Compute_nodal_stresses(int Nb_materials, vector<Material>& Materials) ;
    void Compute_error() ;
    void Compute_mass_scaling(double Target_error, double Inv_Target_error, double Control_parameter_mass_scaling, double Max_mass_scaling, double Error_factor_mass_scaling, double Accepted_ratio, double Decrease_factor_mass_scaling) ;
    void Store() ;
    void Restore() ;
    void Update_damage() ;
    void Update_material(int Nb_materials, vector<Material>& Materials, vector<int> flags) ;
} ;



// Content of the class

// Constructor and Destructor

Body::Body (int i, string m, string p, string t)
{
    index = i ;
    material_name = m ;
    periodicity = p ;
    type = t ;
    drivendof = {0.,0.,0.,0.,0.,0.,0.,0.} ;
    damage = 0. ;
    nb_border_nodes = 0 ;
    Elastic_energy = {0.} ;
    nb_regions = 0 ;

    nodal_distance = 0. ;
    detection_distance = 0. ;
    contact_distance = 0. ;
    mass = 0. ;
    inertia = 0. ;
    inverse_mass = 0. ;
    inverse_inertia = 0. ;
    delta_factor_mass_scaling = 0. ;
    factor_mass_scaling = 1. ;
    x_initial = 0. ;
    y_initial = 0. ;
    r_initial = 0. ;
    x_current = 0. ;
    y_current = 0. ;
    r_current = 0. ;
    stored_x_current = 0. ;
    stored_y_current = 0. ;
    stored_r_current = 0. ;
    x_displacement = 0. ;
    y_displacement = 0. ;
    r_displacement = 0. ;
    x_displacement_temporary = 0. ;
    y_displacement_temporary = 0. ;
    r_displacement_temporary = 0. ;
    stored_x_displacement = 0. ;
    stored_y_displacement = 0. ;
    stored_r_displacement = 0. ;
    x_velocity = 0. ;
    y_velocity = 0. ;
    r_velocity = 0. ;
    x_velocity_temporary = 0. ;
    y_velocity_temporary = 0. ;
    r_velocity_temporary = 0. ;
    stored_x_velocity = 0. ;
    stored_y_velocity = 0. ;
    stored_r_velocity = 0. ;
    x_acceleration = 0. ;
    y_acceleration = 0. ;
    r_acceleration = 0. ;
    stored_x_acceleration = 0. ;
    stored_y_acceleration = 0. ;
    stored_r_acceleration = 0. ;
    x_contact_force = 0. ;
    y_contact_force = 0. ;
    r_contact_force = 0. ;
//    x_master_contact_forces = {0.} ;
//    y_master_contact_forces = {0.} ;
//    r_master_contact_forces = {0.} ;
    x_body_force = 0. ;
    y_body_force = 0. ;
    r_body_force = 0. ;
    x_dirichlet_force = 0. ;
    y_dirichlet_force = 0. ;
    r_dirichlet_force = 0. ;
    x_neumann_force = 0. ;
    y_neumann_force = 0. ;
    r_neumann_force = 0. ;
    x_damping_force = 0. ;
    y_damping_force = 0. ;
    r_damping_force = 0. ;
    x_force = 0. ;
    y_force = 0. ;
    r_force = 0. ;
    unbalanced = 0. ;
    total_error = 0. ;
    max_error = 0. ;
    node_for_max_error = 0 ;
    damage = 0. ;
    internal_work = 0. ;
    contact_work = 0. ;
    body_work = 0. ;
    dirichlet_work = 0. ;
    neumann_work = 0. ;
    damping_work = 0. ;
    alid_work = 0. ;
    kinetic_energy = 0. ;
    deformation_energy = 0. ;
    flag_alid = 0. ;
    alid_xmin = 0. ;
    alid_xmax = 0. ;
    alid_ymin = 0. ;
    alid_ymax = 0. ;
    alid_range = 0. ;
    alid_alphain = 0. ;
    alid_alphaout = 0. ;
    alid_betain = 0. ;
    alid_betaout = 0. ;
    alid_exponent = 0. ;
    stored_internal_work = 0. ;
    stored_contact_work = 0. ;
    stored_body_work = 0. ;
    stored_dirichlet_work = 0. ;
    stored_neumann_work = 0. ;
    stored_damping_work = 0. ;
    stored_alid_work = 0. ;
    stored_x_current = 0. ;
    stored_y_current = 0. ;
    stored_r_current = 0. ;
    stored_x_displacement = 0. ;
    stored_y_displacement = 0. ;
    stored_r_displacement = 0. ;
    stored_x_velocity = 0. ;
    stored_y_velocity = 0. ;
    stored_r_velocity = 0. ;
    density = 0. ;
}

Body::~Body()
{
}

// Accessors

// Modifiers

// Methods



//********************************************//
//** UPDATE BC *******************************//
//********************************************//

void Body::Update_bc(double Time)
{
    for (int n(0) ; n<nb_borders ; n++)
    {
        borders[n].Update_bc(Time) ;
        if (borders[n].x_bc_type == "driven")
        {
            double t0, t1, p0, p1 ;
            int flag_exit = 0 ;
            int nb = 0 ;
            drivendof[0] = 1. ;
            drivendof[1] = 0. ;//drivendof[2] = 0. ;drivendof[3] = 0. ;
            while ( flag_exit == 0 )
            {
                t0 = borders[n].xdirichlet_instants[nb] ;
                p0 = borders[n].xdirichlet_values[nb] ;
                t1 = borders[n].xdirichlet_instants[nb+1] ;
                p1 = borders[n].xdirichlet_values[nb+1] ;
                if ( Time >= t1 )
                {
                    drivendof[1] += ( t1 - t0 ) * ( p1 + p0 ) * 0.5 ;
                }
                else
                {
                    drivendof[1] += ( Time - t0 ) * ( p0 + ( p1 - p0 ) / ( t1 - t0 ) * ( ( Time - t0 ) * 0.5 ) ) ;
                    drivendof[2] = ( Time - t0 ) * ( p1 - p0 ) / ( t1 - t0 ) ;
                    drivendof[3] = ( p1 - p0 ) / ( t1 - t0 ) ;
                    flag_exit = 1 ;
                }
                nb++ ;
            }
        }
        if (borders[n].y_bc_type == "driven")
        {
            double t0, t1, p0, p1 ;
            int flag_exit = 0 ;
            int nb = 0 ;
            drivendof[4] = 1. ;
            drivendof[5] = 0. ;//drivendof[6] = 0. ;drivendof[7] = 0. ;
            while ( flag_exit == 0 )
            {
                t0 = borders[n].ydirichlet_instants[nb] ;
                p0 = borders[n].ydirichlet_values[nb] ;
                t1 = borders[n].ydirichlet_instants[nb+1] ;
                p1 = borders[n].ydirichlet_values[nb+1] ;
                if ( Time >= t1 )
                {
                    drivendof[5] += ( t1 - t0 ) * ( p1 + p0 ) * 0.5 ;
                }
                else
                {
                    drivendof[5] += ( Time - t0 ) * ( p0 + ( p1 - p0 ) / ( t1 - t0 ) * ( ( Time - t0 ) * 0.5 ) ) ;
                    drivendof[6]= ( Time - t0 ) * ( p1 - p0 ) / ( t1 - t0 ) ;
                    drivendof[7] = ( p1 - p0 ) / ( t1 - t0 ) ;
                    flag_exit = 1 ;
                }
                nb++ ;
            }
        }
    }
}



//********************************************//
//** INITIALIZE SEGMENTS *********************//
//********************************************//

void Body::Initialize_segments(double xmin, double xmax)
{
    for (int i(0) ; i<nb_borders ; i++)
    {
        int number_border_nodes = borders[i].number_border_nodes ;
        string periodicity = borders[i].periodicity ;
        vector<double> length(number_border_nodes, 0) ;
        vector<double> initial_length(number_border_nodes, 0) ;
        vector<int> node0(number_border_nodes, 0) ;
        vector<int> shift0(number_border_nodes, 0) ;
        vector<int> node1(number_border_nodes, 0) ;
        vector<int> shift1(number_border_nodes, 0) ;
        vector<int> node2(number_border_nodes, 0) ;
        vector<int> shift2(number_border_nodes, 0) ;
        vector<int> node3(number_border_nodes, 0) ;
        vector<int> shift3(number_border_nodes, 0) ;

        // First node of the border
        int border_before = borders[i].border_before ;
        int border_after = borders[i].border_after ;
        node0[0] = borders[border_before].border_nodes[borders[border_before].number_border_nodes-2] ;
        node1[0] = borders[i].border_nodes[0] ;
        node2[0] = borders[i].border_nodes[1] ;
        node3[0] = borders[i].border_nodes[2] ;
        if (periodicity == "Periodic" && i==0)
            shift0[0] = -1 ;
        else if (periodicity == "Periodic")
            shift0[0] = 1 ;
        /*
        if (periodicity == "Closed")
        {
        	node0[0] = borders[i].border_nodes[number_border_nodes-1] ;
        }
        else if (periodicity == "Simple")
        {
        	if (i==0)
        	{
        		int nbn = borders[nb_borders-1].number_border_nodes-2 ;
        		node0[0] = borders[nb_borders-1].border_nodes[nbn] ;
        	}
        	else
        	{
        		int nbn = borders[i-1].number_border_nodes-2 ;
        		node0[0] = borders[i-1].border_nodes[nbn] ;
        	}
        }
        else if (periodicity == "Periodic")
        {
        	node0[0] = borders[i].border_nodes[number_border_nodes-1] ;
        	if (i==0)
        	{
        		shift0[0] = -1 ;
        	}
        	else
        	{
        		shift0[0] = 1 ;
        	}
        }
        */

        // Current nodes
        for (int j(1) ; j<number_border_nodes-2 ; j++)
        {
            node0[j] = borders[i].border_nodes[j-1] ;
            node1[j] = borders[i].border_nodes[j] ;
            node2[j] = borders[i].border_nodes[j+1] ;
            node3[j] = borders[i].border_nodes[j+2] ;
        }

        // Second to last node of the border
        node0[number_border_nodes-2] = borders[i].border_nodes[number_border_nodes-3] ;
        node1[number_border_nodes-2] = borders[i].border_nodes[number_border_nodes-2] ;
        node2[number_border_nodes-2] = borders[i].border_nodes[number_border_nodes-1] ;
        node3[number_border_nodes-2] = borders[border_after].border_nodes[1] ;
        if (periodicity == "Periodic" && i==0)
        {
            shift2[number_border_nodes-2] = 1 ;
            shift3[number_border_nodes-2] = 1 ;
        }
        else if (periodicity == "Periodic")
        {
            shift2[number_border_nodes-2] = -1 ;
            shift3[number_border_nodes-2] = -1 ;
        }
        /*
        node0[number_border_nodes-2] = borders[i].border_nodes[number_border_nodes-3] ;
        node1[number_border_nodes-2] = borders[i].border_nodes[number_border_nodes-2] ;
        node2[number_border_nodes-2] = borders[i].border_nodes[number_border_nodes-1] ;
        if (periodicity == "Closed")
        {
        	node3[number_border_nodes-2] = borders[i].border_nodes[0] ;
        }
        else if (periodicity == "Simple")
        {
        	if (i==nb_borders-1)
        	{
        		node3[number_border_nodes-2] = borders[0].border_nodes[1] ;
        	}
        	else
        	{
        		node3[number_border_nodes-2] = borders[i+1].border_nodes[1] ;
        	}
        }
        else if (periodicity == "Periodic")
        {
        	node3[number_border_nodes-2] = borders[i].border_nodes[0] ;
        	if (i==0)
        	{
        		shift3[number_border_nodes-2] = 1 ;
        	}
        	else
        	{
        		shift3[number_border_nodes-2] = -1 ;
        	}
        }
        */

        // Last node of the border
        node0[number_border_nodes-1] = borders[i].border_nodes[number_border_nodes-2] ;
        node1[number_border_nodes-1] = borders[i].border_nodes[number_border_nodes-1] ;
        node2[number_border_nodes-1] = borders[border_after].border_nodes[1] ;
        node3[number_border_nodes-1] = borders[border_after].border_nodes[2] ;
        if (periodicity == "Periodic")
        {
            if (i==0)
            {
                shift1[number_border_nodes-1] = 1 ;
                shift2[number_border_nodes-1] = 1 ;
                shift3[number_border_nodes-1] = 1 ;
            }
            else
            {
                shift1[number_border_nodes-1] = -1 ;
                shift2[number_border_nodes-1] = -1 ;
                shift3[number_border_nodes-1] = -1 ;
            }
        }
        /*
        node0[number_border_nodes-1] = borders[i].border_nodes[number_border_nodes-2] ;
        node1[number_border_nodes-1] = borders[i].border_nodes[number_border_nodes-1] ;
        if (periodicity == "Closed")
        {
        	node2[number_border_nodes-1] = borders[i].border_nodes[0] ;
        	node3[number_border_nodes-1] = borders[i].border_nodes[1] ;
        }
        else if (periodicity == "Simple")
        {
        	if (i==nb_borders-1)
        	{
        		node2[number_border_nodes-1] = borders[0].border_nodes[1] ;
        		node3[number_border_nodes-1] = borders[0].border_nodes[2] ;
        	}
        	else
        	{
        		node2[number_border_nodes-1] = borders[i+1].border_nodes[1] ;
        		node3[number_border_nodes-1] = borders[i+1].border_nodes[2] ;
        	}
        }
        else if (periodicity == "Periodic")
        {
        	node2[number_border_nodes-1] = borders[i].border_nodes[0] ;
        	node3[number_border_nodes-1] = borders[i].border_nodes[1] ;
        	if (i==0)
        	{
        		shift2[number_border_nodes-1] = 1 ;
        		shift3[number_border_nodes-1] = 1 ;
        	}
        	else
        	{
        		shift2[number_border_nodes-1] = -1 ;
        		shift3[number_border_nodes-1] = -1 ;
        	}
        }
        */

        // Lengths
        double period = xmax - xmin ;
        double x0, y0, x1, y1, x2, y2 ;
        for (int j(0) ; j<number_border_nodes ; j++)
        {
            x0 = nodes[node0[j]].x_initial + shift0[j] * period ;
            y0 = nodes[node0[j]].y_initial ;
            x1 = nodes[node1[j]].x_initial + shift1[j] * period ;
            y1 = nodes[node1[j]].y_initial ;
            x2 = nodes[node2[j]].x_initial + shift2[j] * period ;
            y2 = nodes[node2[j]].y_initial ;
            length[j] = 0.5 * sqrt((x1-x0)*(x1-x0)+(y1-y0)*(y1-y0)) + 0.5 * sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1)) ;
            initial_length[j] = length[j] ;
        }

        // Update of the border
        borders[i].length = length ;
        borders[i].initial_length = initial_length ;
        borders[i].node0 = node0 ;
        borders[i].shift0 = shift0 ;
        borders[i].node1 = node1 ;
        borders[i].shift1 = shift1 ;
        borders[i].node2 = node2 ;
        borders[i].shift2 = shift2 ;
        borders[i].node3 = node3 ;
        borders[i].shift3 = shift3 ;
    }
}




//********************************************//
//** UPDATE BORDERS **************************//
//********************************************//

void Body::Update_borders(double xmin, double xmax)
{
    if (type=="deformable")
    {
        double x0, y0, x1, y1, x2, y2 ;
        double period = xmax - xmin ;
        for (int i(0) ; i<nb_borders ; i++)
        {
            string border_periodicity = borders[i].periodicity ;
            int number_border_nodes = borders[i].number_border_nodes ;
            vector<double> length(number_border_nodes) ;
            vector<int> node0=borders[i].node0 ;
            vector<int> shift0=borders[i].shift0 ;
            vector<int> node1=borders[i].node1 ;
            vector<int> shift1=borders[i].shift1 ;
            vector<int> node2=borders[i].node2 ;
            vector<int> shift2=borders[i].shift2 ;
            vector<int> node3=borders[i].node3 ;
            vector<int> shift3=borders[i].shift3 ;
            for (int j(0) ; j<number_border_nodes ; j++)
            {
                x0 = nodes[node0[j]].x_current + period * shift0[j] ;
                y0 = nodes[node0[j]].y_current ;
                x1 = nodes[node1[j]].x_current + period * shift1[j] ;
                y1 = nodes[node1[j]].y_current ;
                x2 = nodes[node2[j]].x_current + period * shift2[j] ;
                y2 = nodes[node2[j]].y_current ;
                if (j==0)
                {
                    if (border_periodicity == "Simple")
                        length[j] = sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1))*0.5 ;
                    else
                        length[j] = (sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1))+sqrt((x1-x0)*(x1-x0)+(y1-y0)*(y1-y0)))*0.5 ;
                }
                else if (j==number_border_nodes-1)
                {
                    if (border_periodicity == "Simple")
                        length[j] = sqrt((x1-x0)*(x1-x0)+(y1-y0)*(y1-y0))*0.5 ;
                    else
                        length[j] = (sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1))+sqrt((x1-x0)*(x1-x0)+(y1-y0)*(y1-y0)))*0.5 ;
                }
                else
                    length[j] = (sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1))+sqrt((x1-x0)*(x1-x0)+(y1-y0)*(y1-y0)))*0.5 ;
            }
            borders[i].length = length ;
        }
    }
}



//********************************************//
//** UPDATE BC FORCES ************************//
//********************************************//

void Body::Update_bc_forces()
{
    double penalty, damping ;
    if (type=="deformable")
    {
        string x_bc_type, y_bc_type ;
        double x_bc_value, y_bc_value, x_bc_velocity, y_bc_velocity, fxn, fyn, fxd, fyd ;
        double  tx, ty, tpx, tpy, x0, y0, x2, y2, xp0, yp0, xp2, yp2, invl, invlp ;
        double length ;
        int border_node, number_influencing_nodes ;
        vector<int> border_nodes, influencing_nodes ;
        vector<double> shape_functions ;
        for (int n(0) ; n<nb_nodes ; n++)
        {
            nodes[n].x_neumann_force = 0. ;
            nodes[n].y_neumann_force = 0. ;
            nodes[n].x_dirichlet_force = 0. ;
            nodes[n].y_dirichlet_force = 0. ;
        }
        //if ( == "Simple")   n = nb_borders ;
        //else                n = nb_borders - 1 ;
        for (int i(0) ; i<nb_borders ; i++)
        {
            x_bc_type = borders[i].x_bc_type ;
            y_bc_type = borders[i].y_bc_type ;
            vector<double> px( borders[i].number_border_nodes ) ;
            vector<double> py( borders[i].number_border_nodes ) ;
            if ((x_bc_type == "none") && (y_bc_type == "none"))
            {
                borders[i].x_bc_pressure = px ;
                borders[i].y_bc_pressure = py ;
                continue ;
            }
            x_bc_value = borders[i].x_bc_value ;
            y_bc_value = borders[i].y_bc_value ;
            x_bc_velocity = borders[i].x_bc_velocity ;
            y_bc_velocity = borders[i].y_bc_velocity ;
            border_nodes = borders[i].border_nodes ;
            for (int n(0) ; n<borders[i].number_border_nodes - 1 ; n++)
            {
                border_node = border_nodes[n] ;
                length = borders[i].length[n] ;
                number_influencing_nodes = nodes[border_node].number_influencing_nodes ;
                influencing_nodes = nodes[border_node].influencing_nodes ;
                shape_functions = nodes[border_node].shape_functions ;
                //if (x_bc_type == "oriented" || x_bc_type == "following") //NB : both cases should be treated differently
                //{
                //    px[n] = x_bc_value ;
                //}
                if (x_bc_type == "oriented")
                {
                    px[n] = x_bc_value ;
                }
                else if (x_bc_type == "following")
                {
                    x0 = nodes[borders[i].node0[n]].x_initial ;
                    y0 = nodes[borders[i].node0[n]].y_initial ;
                    x2 = nodes[borders[i].node2[n]].x_initial ;
                    y2 = nodes[borders[i].node2[n]].y_initial ;
                    xp0 = nodes[borders[i].node0[n]].x_current ;
                    yp0 = nodes[borders[i].node0[n]].y_current ;
                    xp2 = nodes[borders[i].node2[n]].x_current ;
                    yp2 = nodes[borders[i].node2[n]].y_current ;
                    invl = pow( (x2 - x0) * (x2 - x0) + (y2 - y0) * (y2 - y0), -0.5 ) ;
                    invlp = pow( (xp2 - xp0) * (xp2 - xp0) + (yp2 - yp0) * (yp2 - yp0), -0.5 ) ;
                    tx = ( x2 - x0 ) * invl ;
                    ty = ( y2 - y0 ) * invl ;
                    tpx = ( xp2 - xp0 ) * invlp ;
                    tpy = ( yp2 - yp0 ) * invlp ;
                    //cout << tx << ' ' << ty << ' ' << tpx << ' ' << tpy << endl ;
                    px[n] += x_bc_value * ( tx * tpx + ty * tpy) ;
                    py[n] += x_bc_value * ( tx * tpy - ty * tpx) ;
                }
                else if (x_bc_type == "soft")
                {
                    penalty = borders[i].xdirichlet_parameters[0] ;
                    damping = borders[i].xdirichlet_parameters[1] ;
                    px[n] = penalty * ( x_bc_value-nodes[border_node].x_displacement ) - damping * ( x_bc_velocity-nodes[border_node].x_velocity ) ; //NB : sign for damping ?
                }
                //if (y_bc_type == "oriented" || y_bc_type == "following")
                //{
                //    py[n] = y_bc_value ;
                //}
                if (y_bc_type == "oriented")
                {
                    py[n] = y_bc_value ;
                }
                else if (y_bc_type == "following")
                {
                    x0 = nodes[borders[i].node0[n]].x_initial ;
                    y0 = nodes[borders[i].node0[n]].y_initial ;
                    x2 = nodes[borders[i].node2[n]].x_initial ;
                    y2 = nodes[borders[i].node2[n]].y_initial ;
                    xp0 = nodes[borders[i].node0[n]].x_current ;
                    yp0 = nodes[borders[i].node0[n]].y_current ;
                    xp2 = nodes[borders[i].node2[n]].x_current ;
                    yp2 = nodes[borders[i].node2[n]].y_current ;
                    invl = pow( (x2 - x0) * (x2 - x0) + (y2 - y0) * (y2 - y0), -0.5 ) ;
                    invlp = pow( (xp2 - xp0) * (xp2 - xp0) + (yp2 - yp0) * (yp2 - yp0), -0.5 ) ;
                    tx = ( x2 - x0 ) * invl ;
                    ty = ( y2 - y0 ) * invl ;
                    tpx = ( xp2 - xp0 ) * invlp ;
                    tpy = ( yp2 - yp0 ) * invlp ;
                    //cout << tx << ' ' << ty << ' ' << tpx << ' ' << tpy << endl ;
                    px[n] += y_bc_value * ( ty * tpx - tx * tpy) ;
                    py[n] += y_bc_value * ( ty * tpy + tx * tpx) ;
                }
                else if (y_bc_type == "soft")
                {
                    penalty = borders[i].ydirichlet_parameters[0] ;
                    damping = borders[i].ydirichlet_parameters[1] ;
                    py[n] = penalty * ( y_bc_value-nodes[border_node].y_displacement ) - damping * ( y_bc_velocity-nodes[border_node].y_velocity ) ; //NB : sign for damping ?
                }
                for (int j(0) ; j<number_influencing_nodes ; j++)
                {
                    fxn = 0. ;
                    fyn = 0. ;
                    fxd = 0. ;
                    fyd = 0. ;
                    if (x_bc_type == "oriented")
                        fxn = px[n] * length * shape_functions[j] ;
                    else if (x_bc_type == "following")
                    {
                        fxn = px[n] * length * shape_functions[j] ;
                        fyn = py[n] * length * shape_functions[j] ;
                    }
                    else if (x_bc_type == "soft")
                        fxd = px[n] * length * shape_functions[j] ;
                    if (y_bc_type == "oriented")
                        fyn = py[n] * length * shape_functions[j] ;
                    else if (y_bc_type == "following")
                    {
                        fxn = px[n] * length * shape_functions[j] ;
                        fyn = py[n] * length * shape_functions[j] ;
                    }
                    else if (y_bc_type == "soft")
                        fyd = py[n] * length * shape_functions[j] ;
                    nodes[influencing_nodes[j]].x_neumann_force += fxn ;
                    nodes[influencing_nodes[j]].y_neumann_force += fyn ;
                    nodes[influencing_nodes[j]].x_dirichlet_force += fxd ;
                    nodes[influencing_nodes[j]].y_dirichlet_force += fyd ;
                }
            }
            borders[i].x_bc_pressure = px ;
            borders[i].y_bc_pressure = py ;
        }
    }
    else if (type=="rigid")
    {
        string x_bc_type, y_bc_type ;
        double x_bc_value, y_bc_value, x_bc_velocity, y_bc_velocity, fxn, fyn, fxd, fyd, dx, dy ;
        vector<int> border_nodes ;
        double length ;
        int border_node ;
        x_neumann_force = 0. ;
        y_neumann_force = 0. ;
        r_neumann_force = 0. ;
        x_dirichlet_force = 0. ;
        y_dirichlet_force = 0. ;
        r_dirichlet_force = 0. ;
        for (int i(0) ; i<nb_borders ; i++)
        {
            x_bc_type = borders[i].x_bc_type ;
            y_bc_type = borders[i].y_bc_type ;
            vector<double> px( borders[i].number_border_nodes ) ;
            vector<double> py( borders[i].number_border_nodes ) ;
            if ((x_bc_type == "none") && (y_bc_type == "none"))
            {
                borders[i].x_bc_pressure = px ;
                borders[i].y_bc_pressure = py ;
                continue ;
            }
            x_bc_value = borders[i].x_bc_value ;
            y_bc_value = borders[i].y_bc_value ;
            x_bc_velocity = borders[i].x_bc_velocity ;
            y_bc_velocity = borders[i].y_bc_velocity ;
            border_nodes = borders[i].border_nodes ;
            for (int n(0) ; n<borders[i].number_border_nodes ; n++)
            {
                fxn = 0. ;
                fyn = 0. ;
                fxd = 0. ;
                fyd = 0. ;
                border_node = border_nodes[n] ;
                length = borders[i].length[n] ;
                if (x_bc_type == "oriented" || x_bc_type == "following")
                {
                    px[n] = x_bc_value ;
                    fxn = px[n] * length ;
                }
                else if (x_bc_type == "soft")
                {
                    penalty = borders[i].xdirichlet_parameters[0] ;
                    damping = borders[i].xdirichlet_parameters[1] ;
                    px[n] = penalty * ( x_bc_value-nodes[border_node].x_displacement ) + damping * ( x_bc_velocity-nodes[border_node].x_velocity ) ; //NB : sign for damping ?
                    fxd = px[n] * length ;
                }
                else if (x_bc_type == "driven")
                {
                    px[n] = 0. ;
                    fxd = 0. ;
                }
                if (y_bc_type == "oriented" || y_bc_type == "following")
                {
                    py[n] = y_bc_value ;
                    fyn = py[n] * length ;
                }
                else if (y_bc_type == "soft")
                {
                    penalty = borders[i].ydirichlet_parameters[0] ;
                    damping = borders[i].ydirichlet_parameters[1] ;
                    py[n] = penalty * ( y_bc_value-nodes[border_node].y_displacement ) + damping * ( y_bc_velocity-nodes[border_node].y_velocity ) ; //NB : sign for damping ?
                    fyd = py[n] * length ;
                }
                else if (y_bc_type == "driven")
                {
                    py[n] = 0. ;
                    fyd = 0. ;
                }
                dx = nodes[border_node].x_current - x_current ;
                dy = nodes[border_node].y_current - y_current ;
                x_neumann_force += fxn ;
                y_neumann_force += fyn ;
                r_neumann_force += -fxn * dy + fyn * dx ;
                x_dirichlet_force += fxd ;
                y_dirichlet_force += fyd ;
                r_dirichlet_force += -fxd * dy + fyd * dx ;
            }
            borders[i].x_bc_pressure = px ;
            borders[i].y_bc_pressure = py ;
        }
    }
}



//********************************************//
//** UPDATE BODY FORCES **********************//
//********************************************//

void Body::Update_body_forces(double Xgravity, double Ygravity)
{
    if (status=="nogravity")
    {
        Xgravity = 0. ;
        Ygravity = 0. ;
    }
    if (type=="deformable")
    {
        for (int i(0) ; i<nb_nodes ; i++)
        {
            nodes[i].x_body_force = nodes[i].x_mass * Xgravity ;
            nodes[i].y_body_force = nodes[i].y_mass * Ygravity ;
        }
    }
    else if (type=="rigid")
    {
        x_body_force = mass * Xgravity ;
        y_body_force = mass * Ygravity ;
    }
}



//********************************************//
//** UPDATE INTERNAL FORCES ******************//
//********************************************//

void Body::Update_internal_forces( int region, int Nb_materials, vector<Material>& Materials )
{
    string material_type ;
    vector<double> parameters ;
    /*
    for (int i(0) ; i<Nb_materials ; i++)
    {
    	if (Materials[i].name == material_name)
    	{
    		material_type=Materials[i].type ;
    		parameters=Materials[i].parameters ;
    		break ;
    	}
    }
    */
    material_type=Materials[material_index].type ;
    parameters=Materials[material_index].parameters ;
    //vector<vector<double>> S={{0,0,0},{0,0,0},{0,0,0}};
    //vector<vector<double>> E={{0,0,0},{0,0,0},{0,0,0}};
    //vector<vector<double>> Eref={{0,0,0},{0,0,0},{0,0,0}};
    //vector<vector<double>> Fref={{1,0,0},{0,1,0},{0,0,1}};
    //vector<vector<double>> F ;
    double E11, E12, E22 ;
    double S11, S12, S22, S33 ;
    double F11, F12, F21, F22 ;
    //
    //double a , b , c , d ;
    int n ;
    //
    Elastic_energy[region] = 0. ;
    int i ;
    double energy, c0(0), J, Mu, Kappa, Lambda ;
    for (int ii(0) ; ii<(int)region_gpoints[region].size() ; ii++)
    {
        i = region_gpoints[region][ii] ;
        if (material_type == "ElasticLinear")
        {
            //E = Eref ;
            E11 = 0. ;
            E12 = 0. ;
            E22 = 0. ;
            for (int j(0) ; j<gpoints[i].number_influencing_nodes ; j++)
            {
                //E[0][0] += gpoints[i].shape_xderiv[j] * nodes[gpoints[i].influencing_nodes[j]].x_displacement_parameter ;
                //E[0][1] += gpoints[i].shape_xderiv[j] * nodes[gpoints[i].influencing_nodes[j]].y_displacement_parameter + gpoints[i].shape_yderiv[j] * nodes[gpoints[i].influencing_nodes[j]].x_displacement_parameter ;
                //E[1][1] += gpoints[i].shape_yderiv[j] * nodes[gpoints[i].influencing_nodes[j]].y_displacement_parameter ;
                E11 += gpoints[i].shape_xderiv[j] * nodes[gpoints[i].influencing_nodes[j]].x_displacement_parameter ;
                E12 += gpoints[i].shape_xderiv[j] * nodes[gpoints[i].influencing_nodes[j]].y_displacement_parameter + gpoints[i].shape_yderiv[j] * nodes[gpoints[i].influencing_nodes[j]].x_displacement_parameter ;
                E22 += gpoints[i].shape_yderiv[j] * nodes[gpoints[i].influencing_nodes[j]].y_displacement_parameter ;
            }
            //E[0][1] = 0.5 * E[0][1] ;
            //E[1][0] = E[0][1] ;
            E12 = 0.5 * E12 ;
            Mu = parameters[3] ;
            Lambda = parameters[4] ;
            Apply_Elastic_Linear( S11, S12, S22, S33, E11, E12, E22, Mu, Lambda, J, energy ) ;
            c0 = gpoints[i].weight * gpoints[i].jacobian ;
            for (int j(0) ; j<gpoints[i].number_influencing_nodes ; j++)
            {
                nodes[gpoints[i].influencing_nodes[j]].x_regions_internal_forces[region] -= c0 * ( gpoints[i].shape_xderiv[j] * S11 + gpoints[i].shape_yderiv[j] * S12 ) ;
                nodes[gpoints[i].influencing_nodes[j]].y_regions_internal_forces[region] -= c0 * ( gpoints[i].shape_yderiv[j] * S22 + gpoints[i].shape_xderiv[j] * S12 ) ;
            }
        }
        else if (material_type == "NeoHookean")
        {
            F11 = 1. ;
            F12 = 0. ;
            F21 = 0. ;
            F22 = 1. ;
            //
            n = gpoints[i].number_influencing_nodes ;
            for (int j(0) ; j<n ; j++)
            {
                F11 += gpoints[i].shape_xderiv[j] * nodes[gpoints[i].influencing_nodes[j]].x_displacement_parameter ;
                F12 += gpoints[i].shape_yderiv[j] * nodes[gpoints[i].influencing_nodes[j]].x_displacement_parameter ;
                F21 += gpoints[i].shape_xderiv[j] * nodes[gpoints[i].influencing_nodes[j]].y_displacement_parameter ;
                F22 += gpoints[i].shape_yderiv[j] * nodes[gpoints[i].influencing_nodes[j]].y_displacement_parameter ;
            }
            /*
            for (int j(0) ; j<gpoints[i].number_influencing_nodes ; j++)
            {
                F11 += gpoints[i].shape_xderiv[j] * nodes[gpoints[i].influencing_nodes[j]].x_displacement_parameter ;
                F12 += gpoints[i].shape_yderiv[j] * nodes[gpoints[i].influencing_nodes[j]].x_displacement_parameter ;
                F21 += gpoints[i].shape_xderiv[j] * nodes[gpoints[i].influencing_nodes[j]].y_displacement_parameter ;
                F22 += gpoints[i].shape_yderiv[j] * nodes[gpoints[i].influencing_nodes[j]].y_displacement_parameter ;
            }
            */
            Mu = parameters[3] ;
            Kappa = parameters[4] ;
            Apply_NeoHookean( S11, S12, S22, S33, F11, F12, F21, F22, Mu, Kappa, J, energy ) ;
            c0 = gpoints[i].weight * gpoints[i].jacobian ;
            //
            for (int j(0) ; j<n ; j++)
            {
                nodes[gpoints[i].influencing_nodes[j]].x_regions_internal_forces[region] -= c0 * ( gpoints[i].shape_xderiv[j] * ( F11 * S11 + F12 * S12 ) +
                        gpoints[i].shape_yderiv[j] * ( F12 * S22 + F11 * S12 ) ) ;
                nodes[gpoints[i].influencing_nodes[j]].y_regions_internal_forces[region] -= c0 * ( gpoints[i].shape_xderiv[j] * ( F21 * S11 + F22 * S12 ) +
                        gpoints[i].shape_yderiv[j] * ( F22 * S22 + F21 * S12 ) ) ;
            }
            /*
            for (int j(0) ; j<gpoints[i].number_influencing_nodes ; j++)
            {
                nodes[gpoints[i].influencing_nodes[j]].x_regions_internal_forces[region] -= c0 * ( gpoints[i].shape_xderiv[j] * ( F11 * S11 + F12 * S12 ) +
                                                                                                   gpoints[i].shape_yderiv[j] * ( F12 * S22 + F11 * S12 ) ) ;
                nodes[gpoints[i].influencing_nodes[j]].y_regions_internal_forces[region] -= c0 * ( gpoints[i].shape_xderiv[j] * ( F21 * S11 + F22 * S12 ) +
                                                                                                   gpoints[i].shape_yderiv[j] * ( F22 * S22 + F21 * S12 ) ) ;
            }
            */
        }
        // OTHER MATERIALS ??
        Elastic_energy[region] += energy * c0 ;
    }
}



//********************************************//
//** UPDATE DAMPING FORCES *******************//
//********************************************//

void Body::Update_damping_forces()
{
    if (type=="deformable")
    {
        int dof0, dof1 ;
        int node0, node1 ;
        int dir0, dir1 ;
        //double d , dx , dy ;
        double d ;
        //vector<int> test10={1,0} ;
        //vector<int> test01={0,1} ;
        //
        for (int i(0) ; i<nb_nodes ; i++)
        {
            nodes[i].x_damping_force = 0. ;
            nodes[i].y_damping_force = 0. ;
        }
        for (int i(0) ; i<nb_matrix ; i++)
        {
            dof0 = matrix_coordinates[i][0] ;
            dof1 = matrix_coordinates[i][1] ;
            d = damping[i] ;
            node0 = dof0 / 2 ;
            dir0 = dof0 - node0 * 2 ;
            node1 = dof1 / 2 ;
            dir1 = dof1 - node1 * 2 ;
            /*
            if ((dir0 == 0) && (dir1 == 0))
            {
                dx = -nodes[node1].x_velocity_parameter * d ;
                nodes[node0].x_damping_force += dx ;
                if (dof0 != dof1)
                {
                    dx = -nodes[node0].x_velocity_parameter * d ;
                    nodes[node1].x_damping_force += dx ;
                }
            }
            else if ((dir0 == 0) && (dir1 == 1))
            {
                dx = -nodes[node1].y_velocity_parameter * d ;
                nodes[node0].x_damping_force += dx ;
                dy = -nodes[node0].x_velocity_parameter * d ;
                nodes[node1].y_damping_force += dy ;
            }
            else if ((dir0 == 1) && (dir1 == 0))
            {
                dy = -nodes[node1].x_velocity_parameter * d ;
                nodes[node0].y_damping_force += dy ;
                dx = -nodes[node0].y_velocity_parameter * d ;
                nodes[node1].x_damping_force += dx ;
            }
            else if ((dir0 == 1) && (dir1 == 1))
            {
                dy = -nodes[node1].y_velocity_parameter * d ;
                nodes[node0].y_damping_force += dy ;
                if (dof0 != dof1)
                {
                    dy = -nodes[node0].y_velocity_parameter * d ;
                    nodes[node1].y_damping_force += dy ;
                }
            }
            */
            /*
            nodes[node0].x_damping_force -= ( test10[dir0] * test10[dir1] * nodes[node1].x_velocity_parameter + test10[dir0] * test01[dir1] * nodes[node1].y_velocity_parameter ) * d ;
            nodes[node0].y_damping_force -= ( test01[dir0] * test10[dir1] * nodes[node1].x_velocity_parameter + test01[dir0] * test01[dir1] * nodes[node1].y_velocity_parameter ) * d ;
            nodes[node1].x_damping_force -= test01[dir0] * test10[dir1] * nodes[node0].y_velocity_parameter * d ;
            nodes[node1].y_damping_force -= test10[dir0] * test01[dir1] * nodes[node0].x_velocity_parameter * d ;
            if (dof0 != dof1)
            {
                nodes[node1].x_damping_force -= test10[dir0] * test10[dir1] * nodes[node0].x_velocity_parameter * d ;
                nodes[node1].y_damping_force -= test01[dir0] * test01[dir1] * nodes[node0].y_velocity_parameter * d ;
            }
            */
            if ((dir0 == 0) && (dir1 == 0))
            {
                nodes[node0].x_damping_force -= nodes[node1].x_velocity_parameter * d ;
                if (dof0 != dof1)
                    nodes[node1].x_damping_force -= nodes[node0].x_velocity_parameter * d ;
            }
            else if ((dir0 == 0) && (dir1 == 1))
            {
                nodes[node0].x_damping_force -= nodes[node1].y_velocity_parameter * d ;
                nodes[node1].y_damping_force -= nodes[node0].x_velocity_parameter * d ;
            }
            else if ((dir0 == 1) && (dir1 == 0))
            {
                nodes[node0].y_damping_force -= nodes[node1].x_velocity_parameter * d ;
                nodes[node1].x_damping_force -= nodes[node0].y_velocity_parameter * d ;
            }
            else if ((dir0 == 1) && (dir1 == 1))
            {
                nodes[node0].y_damping_force -= nodes[node1].y_velocity_parameter * d ;
                if (dof0 != dof1)
                    nodes[node1].y_damping_force -= nodes[node0].y_velocity_parameter * d ;
            }
        }
    }
}




//********************************************//
//** UPDATE ALID FORCES **********************//
//********************************************//

void Body::Update_alid_forces(double Time)
{
    if (type=="rigid")
        return ;
    else if (flag_alid==0)
        return ;
    double beta, pvx(0), pvy(0) ;
    int flag_alidx, flag_alidy ;
    int nb ;
    double t0, t1, p0, p1 ;
    int flag_exit ;
    if ((int)xalid_instants.size()==0)
        flag_alidx = 0 ;
    else
    {
        flag_exit = 0 ;
        flag_alidx = 1 ;
        nb = 0 ;
        while ( flag_exit == 0 )
        {
            t0 = xalid_instants[nb] ;
            p0 = xalid_values[nb] ;
            t1 = xalid_instants[nb+1] ;
            p1 = xalid_values[nb+1] ;
            if ( (t0 <= Time) & (Time <= t1) )
            {
                pvx = p0 + ( p1 - p0 ) / ( t1 - t0 ) * ( Time - t0 ) ;
                flag_exit = 1 ;
            }
            nb++ ;
        }
    }
    if ((int)yalid_instants.size()==0)
        flag_alidy = 0 ;
    else
    {
        flag_exit = 0 ;
        flag_alidy = 1 ;
        nb = 0 ;
        while ( flag_exit == 0 )
        {
            t0 = yalid_instants[nb] ;
            p0 = yalid_values[nb] ;
            t1 = yalid_instants[nb+1] ;
            p1 = yalid_values[nb+1] ;
            if ( (t0 <= Time) & (Time <= t1) )
            {
                pvy = p0 + ( p1 - p0 ) / ( t1 - t0 ) * ( Time - t0 ) ;
                flag_exit = 1 ;
            }
            nb++ ;
        }
    }
    for (int i(0) ; i<nb_nodes ; i++)
    {
        beta = nodes[i].alid_beta ;
        if (flag_alidx==0)
            nodes[i].x_alid_force = 0. ;
        else
            nodes[i].x_alid_force = - beta * nodes[i].x_mass * ( nodes[i].x_velocity - pvx ) ;
        if (flag_alidy==0)
            nodes[i].y_alid_force = 0. ;
        else
            nodes[i].y_alid_force = - beta * nodes[i].y_mass * ( nodes[i].y_velocity - pvy ) ;
    }

    int dof0, dof1 ;
    int node0, node1 ;
    int dir0, dir1 ;
    double d, dx, dy ;
    for (int i(0) ; i<nb_matrix ; i++)
    {
        dof0 = matrix_coordinates[i][0] ;
        dof1 = matrix_coordinates[i][1] ;
        node0 = dof0 / 2 ;
        dir0 = dof0 - node0 * 2 ;
        node1 = dof1 / 2 ;
        dir1 = dof1 - node1 * 2 ;
        d = stiffness[i] / (1. / nodes[node0].alid_alpha + 1. / nodes[node1].alid_alpha) ;
        //d = stiffness[i] * 0.5 * (nodes[node0].alid_alpha + nodes[node1].alid_alpha) ;
        if ((dir0 == 0) && (dir1 == 0))
        {
            if (flag_alidx==1)
            {
                dx = -(nodes[node1].x_velocity_parameter - pvx) * d ;
                nodes[node0].x_alid_force += dx ;
            }
            if ((dof0 != dof1) && (flag_alidx==1))
            {
                dx = -(nodes[node0].x_velocity_parameter - pvx) * d ;
                nodes[node1].x_alid_force += dx ;
            }
        }
        else if ((dir0 == 0) && (dir1 == 1))
        {
            if (flag_alidx==1)
            {
                dx = -(nodes[node1].y_velocity_parameter - pvy) * d ;
                nodes[node0].x_alid_force += dx ;
            }
            if (flag_alidy==1)
            {
                dy = -(nodes[node0].x_velocity_parameter - pvx) * d ;
                nodes[node1].y_alid_force += dy ;
            }
        }
        else if ((dir0 == 1) && (dir1 == 0))
        {
            if (flag_alidy==1)
            {
                dy = -(nodes[node1].x_velocity_parameter - pvx) * d ;
                nodes[node0].y_alid_force += dy ;
            }
            if (flag_alidx==1)
            {
                dx = -(nodes[node0].y_velocity_parameter - pvy) * d ;
                nodes[node1].x_alid_force += dx ;
            }
        }
        else if ((dir0 == 1) && (dir1 == 1))
        {
            if (flag_alidy==1)
            {
                dy = -(nodes[node1].y_velocity_parameter - pvy) * d ;
                nodes[node0].y_alid_force += dy ;
            }
            if ((dof0 != dof1) && (flag_alidy==1))
            {
                dy = -(nodes[node0].y_velocity_parameter - pvy) * d ;
                nodes[node1].y_alid_force += dy ;
            }
        }
    }
}




//********************************************//
//** INITIALIZE CONTACT FORCES ***************//
//********************************************//

void Body::Initialize_contact_forces()
{
    if (type=="deformable")
    {
        for (int i(0) ; i<nb_nodes ; i++)
        {
            nodes[i].x_contact_force = 0. ;
            nodes[i].y_contact_force = 0. ;
        }
        contact_forces_to_send = {} ;
        contact_forces_to_send_to = {} ;
    }
    else if (type=="rigid")
    {
        x_contact_force = 0. ;
        y_contact_force = 0. ;
        r_contact_force = 0. ;
        contact_forces_to_send = {} ;
        contact_forces_to_send_to = {} ;
    }
    //for (int i(0) ; i<nb_borders ; i++)
    //{
    //	for (int j(0) ; j<borders[i].number_border_nodes ; j++)
    //	{
    //		borders[i].x_contact_pressure[j] = 0. ;
    //		borders[i].y_contact_pressure[j] = 0. ;
    //	}
    //}
}



//********************************************//
//** UPDATE CONTACTS *************************//
//********************************************//

void Body::Update_contacts(vector<Body>& Bodies, double xmin, double xmax)
{
    double period = xmax - xmin ;
    vector<int> border_nodesS ;
    int bodyM, borderM, shiftM, segmentM ;
    double gapn, gapt, xsi, xnorm, ynorm, xtan, ytan, xclosest, yclosest, xsi_previous ;
    double vgapn, vgapt, vgapx, vgapy ;
    int    nodeS, nodeM0, shiftM0, nodeM1, shiftM1, nodeM2, shiftM2, nodeM3, shiftM3 ;
    double shapeM0, shapeM1, shapeM2, shapeM3 ;
    string periodicity, interpolantM ;
    double xs, ys, x0, y0, x1, y1, x2, y2, x3, y3 ;
    int flag_detect ;
    for (int icontact(0) ; icontact<nb_contact_elements ; icontact++)
    {
        nodeS = contact_elements[icontact].nodeS ;
        bodyM = contact_elements[icontact].bodyM ;
        borderM = contact_elements[icontact].borderM ;
        interpolantM = Bodies[bodyM].borders[borderM].interpolant ;
        shiftM = contact_elements[icontact].shiftM ;
        segmentM = contact_elements[icontact].segmentM ;
        periodicity = Bodies[bodyM].borders[borderM].periodicity ;
        gapn = contact_elements[icontact].gapn ;
        gapt = contact_elements[icontact].gapt ;
        xsi = contact_elements[icontact].xsi ;
        xnorm = contact_elements[icontact].xnorm ;
        ynorm = contact_elements[icontact].ynorm ;
        xtan = contact_elements[icontact].xtan ;
        ytan = contact_elements[icontact].ytan ;
        nodeM0 = Bodies[bodyM].borders[borderM].node0[segmentM] ;
        shiftM0 = Bodies[bodyM].borders[borderM].shift0[segmentM] ;
        nodeM1 = Bodies[bodyM].borders[borderM].node1[segmentM] ;
        shiftM1 = Bodies[bodyM].borders[borderM].shift1[segmentM] ;
        nodeM2 = Bodies[bodyM].borders[borderM].node2[segmentM] ;
        shiftM2 = Bodies[bodyM].borders[borderM].shift2[segmentM] ;
        nodeM3 = Bodies[bodyM].borders[borderM].node3[segmentM] ;
        shiftM3 = Bodies[bodyM].borders[borderM].shift3[segmentM] ;
        xsi_previous = xsi ;
        xs = nodes[nodeS].x_current ;
        ys = nodes[nodeS].y_current ;
        x0 = Bodies[bodyM].nodes[nodeM0].x_current + period * ( shiftM + shiftM0 ) ;
        y0 = Bodies[bodyM].nodes[nodeM0].y_current ;
        x1 = Bodies[bodyM].nodes[nodeM1].x_current + period * ( shiftM + shiftM1 ) ;
        y1 = Bodies[bodyM].nodes[nodeM1].y_current ;
        x2 = Bodies[bodyM].nodes[nodeM2].x_current + period * ( shiftM + shiftM2 ) ;
        y2 = Bodies[bodyM].nodes[nodeM2].y_current ;
        x3 = Bodies[bodyM].nodes[nodeM3].x_current + period * ( shiftM + shiftM3 ) ;
        y3 = Bodies[bodyM].nodes[nodeM3].y_current ;
        if (index==bodyM && (nodeM0==nodeS || nodeM3==nodeS))
        {
            Closest_2_segment_self(xs, ys, x1, y1, x2, y2,
                                   1.e-16, gapn, xsi,
                                   xnorm, ynorm, xtan, ytan, xclosest, yclosest,
                                   shapeM0, shapeM1, shapeM2, shapeM3, flag_detect ) ;
        }
        else
        {
            Closest_2_segment(interpolantM, xs, ys, x0, y0, x1, y1, x2, y2, x3, y3,
                              1.e-16, gapn, xsi,
                              xnorm, ynorm, xtan, ytan, xclosest, yclosest,
                              shapeM0, shapeM1, shapeM2, shapeM3, flag_detect ) ;
        }
        if (flag_detect == -2 || flag_detect == -3)
        {
            gapt = gapt + 0.5 * ( -1. - xsi_previous ) * sqrt( (x2-x1)*(x2-x1) + (y2-y1)*(y2-y1) ) ;
            //
            if (periodicity=="Simple")
            {
                if (segmentM==0)
                {
                    borderM = Bodies[bodyM].borders[borderM].border_before ;
                    segmentM = Bodies[bodyM].borders[borderM].number_border_nodes - 2 ;
                }
                else
                    segmentM-- ;
            }
            else if (periodicity=="Periodic")
            {
                if (segmentM==0)
                {
                    if (borderM==0)
                        shiftM-- ;
                    else if (borderM==1)
                        shiftM++ ;
                    segmentM = Bodies[bodyM].borders[borderM].number_border_nodes - 2 ;
                }
                else
                    segmentM-- ;
            }
            else
            {
                if (segmentM==0)
                    segmentM = Bodies[bodyM].borders[borderM].number_border_nodes - 2 ;
                else
                    segmentM-- ;
            }
            /*
            if (periodicity=="Simple")
            {
                if (segmentM==0)
                {
            		if (borderM==0) borderM = Bodies[bodyM].nb_borders - 1 ;
            		else borderM-- ;
                    segmentM = Bodies[bodyM].borders[borderM].number_border_nodes - 1 ;
                }
                else segmentM-- ;
            }
            else if (periodicity=="Periodic")
            {
                if (segmentM==0)
                {
                    if (borderM==0) shiftM-- ;
                    else if (borderM==1) shiftM++ ;
                    segmentM = Bodies[bodyM].borders[borderM].number_border_nodes - 1 ;
                }
                else segmentM-- ;
            }
            else
            {
                if (segmentM==0) segmentM = Bodies[bodyM].borders[borderM].number_border_nodes - 1 ;
                else segmentM-- ;
            }
            */
            interpolantM = Bodies[bodyM].borders[borderM].interpolant ;
            nodeM0 = Bodies[bodyM].borders[borderM].node0[segmentM] ;
            shiftM0 = Bodies[bodyM].borders[borderM].shift0[segmentM] ;
            nodeM1 = Bodies[bodyM].borders[borderM].node1[segmentM] ;
            shiftM1 = Bodies[bodyM].borders[borderM].shift1[segmentM] ;
            nodeM2 = Bodies[bodyM].borders[borderM].node2[segmentM] ;
            shiftM2 = Bodies[bodyM].borders[borderM].shift2[segmentM] ;
            nodeM3 = Bodies[bodyM].borders[borderM].node3[segmentM] ;
            shiftM3 = Bodies[bodyM].borders[borderM].shift3[segmentM] ;
            x0 = Bodies[bodyM].nodes[nodeM0].x_current + period * ( shiftM + shiftM0 ) ;
            y0 = Bodies[bodyM].nodes[nodeM0].y_current ;
            x1 = Bodies[bodyM].nodes[nodeM1].x_current + period * ( shiftM + shiftM1 ) ;
            y1 = Bodies[bodyM].nodes[nodeM1].y_current ;
            x2 = Bodies[bodyM].nodes[nodeM2].x_current + period * ( shiftM + shiftM2 ) ;
            y2 = Bodies[bodyM].nodes[nodeM2].y_current ;
            x3 = Bodies[bodyM].nodes[nodeM3].x_current + period * ( shiftM + shiftM3 ) ;
            y3 = Bodies[bodyM].nodes[nodeM3].y_current ;
            xsi = 1. ;
            if (index==bodyM && (nodeM0==nodeS || nodeM3==nodeS))
            {
                Closest_2_segment_self(xs, ys, x1, y1, x2, y2,
                                       1.e-16, gapn, xsi,
                                       xnorm, ynorm, xtan, ytan, xclosest, yclosest,
                                       shapeM0, shapeM1, shapeM2, shapeM3, flag_detect ) ;
            }
            else
            {
                Closest_2_segment(interpolantM, xs, ys, x0, y0, x1, y1, x2, y2, x3, y3,
                                  1.e-16, gapn, xsi,
                                  xnorm, ynorm, xtan, ytan, xclosest, yclosest,
                                  shapeM0, shapeM1, shapeM2, shapeM3, flag_detect ) ;
            }
            gapt = gapt + 0.5 * ( -1. + xsi ) * sqrt( (x2-x1)*(x2-x1) + (y2-y1)*(y2-y1) ) ;
            contact_elements[icontact].bodyM = bodyM ;
            contact_elements[icontact].borderM = borderM ;
            contact_elements[icontact].shiftM = shiftM ;
            contact_elements[icontact].segmentM = segmentM ;
        }
        else if (flag_detect == 2 || flag_detect == 3)
        {
            gapt = gapt + 0.5 * ( -1. + xsi_previous ) * sqrt( (x2-x1)*(x2-x1) + (y2-y1)*(y2-y1) ) ;
            //
            if (periodicity=="Simple")
            {
                if (segmentM==Bodies[bodyM].borders[borderM].number_border_nodes - 2)
                {
                    borderM = Bodies[bodyM].borders[borderM].border_after ;
                    segmentM = 0 ;
                }
                else
                    segmentM++ ;
            }
            else if (periodicity=="Periodic")
            {
                if (segmentM==Bodies[bodyM].borders[borderM].number_border_nodes - 2)
                {
                    if (borderM==0)
                        shiftM++ ;
                    else if (borderM==1)
                        shiftM-- ;
                    segmentM = 0 ;
                }
                else
                    segmentM++ ;
            }
            else
            {
                if (segmentM==Bodies[bodyM].borders[borderM].number_border_nodes - 2)
                    segmentM = 0 ;
                else
                    segmentM++ ;
            }
            /*
            if (periodicity=="Simple")
            {
                if (segmentM==Bodies[bodyM].borders[borderM].number_border_nodes - 1)
                {
                    if (borderM==Bodies[bodyM].nb_borders - 1 ) borderM = 0 ;
                    else borderM++ ;
                    segmentM = 0 ;
                }
                else segmentM++ ;
            }
            else if (periodicity=="Periodic")
            {
                if (segmentM==Bodies[bodyM].borders[borderM].number_border_nodes - 1)
                {
                    if (borderM==0) shiftM++ ;
                    else if (borderM==1) shiftM-- ;
                    segmentM = 0 ;
                }
                else segmentM++ ;
            }
            else
            {
                if (segmentM==Bodies[bodyM].borders[borderM].number_border_nodes - 1) segmentM = 0 ;
                else segmentM++ ;
            }
            */
            interpolantM = Bodies[bodyM].borders[borderM].interpolant ;
            nodeM0 = Bodies[bodyM].borders[borderM].node0[segmentM] ;
            shiftM0 = Bodies[bodyM].borders[borderM].shift0[segmentM] ;
            nodeM1 = Bodies[bodyM].borders[borderM].node1[segmentM] ;
            shiftM1 = Bodies[bodyM].borders[borderM].shift1[segmentM] ;
            nodeM2 = Bodies[bodyM].borders[borderM].node2[segmentM] ;
            shiftM2 = Bodies[bodyM].borders[borderM].shift2[segmentM] ;
            nodeM3 = Bodies[bodyM].borders[borderM].node3[segmentM] ;
            shiftM3 = Bodies[bodyM].borders[borderM].shift3[segmentM] ;
            x0 = Bodies[bodyM].nodes[nodeM0].x_current + period * ( shiftM + shiftM0 ) ;
            y0 = Bodies[bodyM].nodes[nodeM0].y_current ;
            x1 = Bodies[bodyM].nodes[nodeM1].x_current + period * ( shiftM + shiftM1 ) ;
            y1 = Bodies[bodyM].nodes[nodeM1].y_current ;
            x2 = Bodies[bodyM].nodes[nodeM2].x_current + period * ( shiftM + shiftM2 ) ;
            y2 = Bodies[bodyM].nodes[nodeM2].y_current ;
            x3 = Bodies[bodyM].nodes[nodeM3].x_current + period * ( shiftM + shiftM3 ) ;
            y3 = Bodies[bodyM].nodes[nodeM3].y_current ;
            xsi = -1. ;
            if (index==bodyM && (nodeM0==nodeS || nodeM3==nodeS))
            {
                Closest_2_segment_self(xs, ys, x1, y1, x2, y2,
                                       1.e-16, gapn, xsi,
                                       xnorm, ynorm, xtan, ytan, xclosest, yclosest,
                                       shapeM0, shapeM1, shapeM2, shapeM3, flag_detect ) ;
            }
            else
            {
                Closest_2_segment(interpolantM, xs, ys, x0, y0, x1, y1, x2, y2, x3, y3,
                                  1.e-16, gapn, xsi,
                                  xnorm, ynorm, xtan, ytan, xclosest, yclosest,
                                  shapeM0, shapeM1, shapeM2, shapeM3, flag_detect ) ;
            }
            gapt = gapt + 0.5 * ( 1. + xsi ) * sqrt( (x2-x1)*(x2-x1) + (y2-y1)*(y2-y1) ) ;
            contact_elements[icontact].bodyM = bodyM ;
            contact_elements[icontact].borderM = borderM ;
            contact_elements[icontact].shiftM = shiftM ;
            contact_elements[icontact].segmentM = segmentM ;
        }
        else
            gapt = gapt + 0.5 * ( xsi - xsi_previous ) * sqrt( (x2-x1)*(x2-x1) + (y2-y1)*(y2-y1) ) ;
        vgapx = nodes[nodeS].x_velocity - ( shapeM0 * Bodies[bodyM].nodes[nodeM0].x_velocity +
                                            shapeM1 * Bodies[bodyM].nodes[nodeM1].x_velocity +
                                            shapeM2 * Bodies[bodyM].nodes[nodeM2].x_velocity +
                                            shapeM3 * Bodies[bodyM].nodes[nodeM3].x_velocity ) ;
        vgapy = nodes[nodeS].y_velocity - ( shapeM0 * Bodies[bodyM].nodes[nodeM0].y_velocity +
                                            shapeM1 * Bodies[bodyM].nodes[nodeM1].y_velocity +
                                            shapeM2 * Bodies[bodyM].nodes[nodeM2].y_velocity +
                                            shapeM3 * Bodies[bodyM].nodes[nodeM3].y_velocity ) ;
        vgapn = vgapx * xnorm + vgapy * ynorm ;
        vgapt = vgapx * xtan + vgapy * ytan ;
        contact_elements[icontact].gapn = gapn ;
        contact_elements[icontact].vgapn = vgapn ;
        contact_elements[icontact].gapt = gapt ;
        contact_elements[icontact].vgapt = vgapt ;
        contact_elements[icontact].xsi = xsi ;
        contact_elements[icontact].xnorm = xnorm ;
        contact_elements[icontact].ynorm = ynorm ;
        contact_elements[icontact].xtan = xtan ;
        contact_elements[icontact].ytan = ytan ;
        contact_elements[icontact].nodeM0 = nodeM0 ;
        contact_elements[icontact].shapeM0 = shapeM0 ;
        contact_elements[icontact].shiftM0 = shiftM0 ;
        contact_elements[icontact].nodeM1 = nodeM1 ;
        contact_elements[icontact].shapeM1 = shapeM1 ;
        contact_elements[icontact].shiftM1 = shiftM1 ;
        contact_elements[icontact].nodeM2 = nodeM2 ;
        contact_elements[icontact].shapeM2 = shapeM2 ;
        contact_elements[icontact].shiftM2 = shiftM2 ;
        contact_elements[icontact].nodeM3 = nodeM3 ;
        contact_elements[icontact].shapeM3 = shapeM3 ;
        contact_elements[icontact].shiftM3 = shiftM3 ;
    }
}



//********************************************//
//** UPDATE CONTACT FORCES *******************//
//********************************************//

void Body::Update_contact_forces( double Deltat, vector<Body>& Bodies, int Nb_contact_laws, vector<Contact_law>& Contact_laws, vector<vector<int>>& Contacts_Table, double xmin, double xmax )
{
    int bodyM, shiftM, borderS, border_nodeS, nodeS, nodeM0, nodeM1, nodeM2, nodeM3, Contact_law_index; // neig ;
    double shapeM0, shapeM1, shapeM2, shapeM3 ;
    double  gapn, vgapn, gapt, vgapt, xnorm, ynorm, xtan, ytan, effective_mass ;
    double Pn, Pt, Px, Py, Fsx(0), Fsy(0), length, initial_length, dx, dy, W ;
    string contact_law_type, length_evolution, material_nameM ;
    vector<double> parameters ;
    string material1, material2 ;
    int number_influencing_nodes ;// , flag_no_law ;
    vector<int> influencing_nodes ;
    vector<double> shape_functions ;
    double period = xmax - xmin ;
    for (int icontact(0) ; icontact<nb_contact_elements ; icontact++)
    {
        W = 0. ;
        borderS = contact_elements[icontact].borderS ;
        border_nodeS = contact_elements[icontact].border_nodeS ;
        length = borders[borderS].length[border_nodeS] ;
        initial_length = borders[borderS].initial_length[border_nodeS] ;
        gapn = contact_elements[icontact].gapn ;
        //if (gapn > length * 0.2 || gapn < -length * 0.05) continue ;
        vgapn = contact_elements[icontact].vgapn ;
        nodeS = contact_elements[icontact].nodeS ;
        bodyM = contact_elements[icontact].bodyM ;
        shiftM = contact_elements[icontact].shiftM ;
        gapt = contact_elements[icontact].gapt ;
        vgapt = contact_elements[icontact].vgapt ;
        xnorm = contact_elements[icontact].xnorm ;
        ynorm = contact_elements[icontact].ynorm ;
        xtan = contact_elements[icontact].xtan ;
        ytan = contact_elements[icontact].ytan ;
        effective_mass = contact_elements[icontact].effective_mass ;
        material_nameM = Bodies[bodyM].material_name ;
        /*
        flag_no_law = 1 ;
        for (int i(0) ; i<Nb_contact_laws ; i++)
        {
        	material1 = Contact_laws[i].material1 ;
        	material2 = Contact_laws[i].material2 ;
        	if (((material1==material_name) && (material2==material_nameM)) ||
        	    ((material2==material_name) && (material1==material_nameM)))
        	{
        		contact_law_type = Contact_laws[i].type ;
        		parameters = Contact_laws[i].parameters ;
        		length_evolution = Contact_laws[i].length_evolution ;
        		flag_no_law = 0 ;
        		break ;
        	}
        }
        if (flag_no_law ==1) continue ;
        */
        Contact_law_index = Contacts_Table[material_index][Bodies[bodyM].material_index] ;
        if (Contact_law_index > -1)
        {
            contact_law_type = Contact_laws[Contact_law_index].type ;
            parameters = Contact_laws[Contact_law_index].parameters ;
            length_evolution = Contact_laws[Contact_law_index].length_evolution ;
        }
        else
            continue ;
        if (contact_law_type == "Frictionless")
        {
            double kn = parameters[0] ;
            Apply_Frictionless( kn, gapn, gapt, Pn, Pt) ;
        }
        else if (contact_law_type == "Cohesion")
        {
            double kn = parameters[0] ;
            double kt = parameters[1] ;
            double gtang = parameters[2] ;
            double gtens = parameters[3] ;
            Apply_Cohesion( kn, kt, gtang*0.5, gtens*0.5, gapn, gapt, Pn, Pt) ;
        }
        else if (contact_law_type == "MohrCoulomb")
        {
            double kn = parameters[0] ;
            double kt = parameters[1] ;
            double fric = parameters[2] ;
            double coh = parameters[3] ;
            double tens = parameters[4] ;
            Apply_Mohr_Coulomb( kn, kt, fric, coh*0.5, tens*0.5, gapn, gapt, Pn, Pt) ;
        }
        else if (contact_law_type == "DampedMohrCoulomb")
        {
            double kn = parameters[0] ;
            double kt = parameters[1] ;
            double fric = parameters[2] ;
            double coh = parameters[3] ;
            double tens = parameters[4] ;
            double damp = parameters[5] ;
            Apply_Damped_Mohr_Coulomb(kn, kt, fric, coh*0.5, tens*0.5, damp, effective_mass, initial_length, gapn, vgapn, gapt, vgapt, Pn, Pt, W) ;
        }
        else if (contact_law_type == "TwoSlopesMohrCoulomb")
        {
            double kn = parameters[0] ;
            double kt = parameters[1] ;
            double fric = parameters[2] ;
            double coh = parameters[3] ;
            double tens = parameters[4] ;
            double rati = parameters[5] ;
            Apply_Two_Slopes_Mohr_Coulomb(kn, kt, fric, coh*0.5, tens*0.5, rati, gapn, vgapn, gapt, vgapt, Pn, Pt) ;
        }
        else if (contact_law_type == "BondedMohrCoulomb")
        {
            double kn = parameters[0] ;
            double kt = parameters[1] ;
            double kbond = parameters[2] ;
            double fricbond = parameters[3] ;
            double cohbond = parameters[4] ;
            double tensbond = parameters[5] ;
            double fricfree = parameters[6] ;
            double cohfree = parameters[7] ;
            double tensfree = parameters[8] ;
            double damp = parameters[9] ;
            int nb_internal = contact_elements[icontact].nb_internal ;
            double Damage = contact_elements[icontact].internal[0] ;
            Pn = contact_elements[icontact].internal[1] ;
            Pt = contact_elements[icontact].internal[2] ;
            if (nb_internal==0) // NB this would mean that the contact was just created and should be initialized to the state "broken bond"
            {
                Damage = 1. ;
                contact_elements[icontact].nb_internal = 3 ;
                contact_elements[icontact].internal = {1., 0., 0.} ;
            }
            Apply_Bonded_Mohr_Coulomb(kn, kt, kbond, fricbond, cohbond*0.5, tensbond*0.5, fricfree, cohfree*0.5, tensfree*0.5, damp, effective_mass, initial_length, gapn, vgapn, gapt, vgapt, Damage, Pn, Pt) ;
            contact_elements[icontact].internal[0] = Damage ;
            contact_elements[icontact].internal[1] = Pn ;
            contact_elements[icontact].internal[2] = Pt ;
        }
        else if (contact_law_type == "CZMlinear")
        {
            double kini = parameters[0] ;
            double Pnlim = parameters[1] ;
            double gapnlim = parameters[2] ;
            double Pnres = parameters[3] ;
            double NTratio = parameters[4] ;
            int nb_internal = contact_elements[icontact].nb_internal ;
            double Damage = contact_elements[icontact].internal[0] ;
            if (nb_internal==0) // NB: this would mean that the contact was just created and should be initialized to the state "broken bond"
            {
                Damage = 1. ;
                contact_elements[icontact].nb_internal = 1 ;
                contact_elements[icontact].internal = {1.} ;
            }
            Apply_CZM_linear(kini, Pnlim*0.5, gapnlim, Pnres*0.5, NTratio, gapn, gapt, Damage, Pn, Pt) ;
            contact_elements[icontact].internal[0] = Damage ;
        }
        else if (contact_law_type == "CZMfatigue")
        {
            double kini = parameters[0] ;
            double Pnlim = parameters[1] ;
            double gapnlim = parameters[2] ;
            double Pnres = parameters[3] ;
            double NTratio = parameters[4] ;
            double Pnfat = parameters[5] ;
            double Raten = parameters[6] ;
            double Ratet = parameters[7] ;
            int nb_internal = contact_elements[icontact].nb_internal ;
            double Damage = contact_elements[icontact].internal[0] ;
            Pn = contact_elements[icontact].internal[1] ;
            Pt = contact_elements[icontact].internal[2] ;
            if (nb_internal==0) // NB this would mean that the contact was just created and should be initialized to the state "broken bond"
            {
                Damage = 1. ;
                contact_elements[icontact].nb_internal = 3 ;
                contact_elements[icontact].internal = {1., 0., 0.} ;
            }
            //cout << "before " << gapn << ' ' << gapt << ' ' << Damage << ' ' << Pn << ' ' << Pt << endl ;
            Apply_CZM_fatigue(kini, Pnlim*0.5, gapnlim, Pnres*0.5, NTratio, Pnfat, Raten, Ratet, gapn, gapt, Damage, Pn, Pt) ;
            //cout << "after  " << gapn << ' ' << gapt << ' ' << Damage << ' ' << Pn << ' ' << Pt << endl ;
            contact_elements[icontact].internal[0] = Damage ;
            contact_elements[icontact].internal[1] = Pn ;
            contact_elements[icontact].internal[2] = Pt ;
        }
        // OTHER CONTACT LAWS ??

        nodeM0 = contact_elements[icontact].nodeM0 ;
        shapeM0 = contact_elements[icontact].shapeM0 ;
        nodeM1 = contact_elements[icontact].nodeM1 ;
        shapeM1 = contact_elements[icontact].shapeM1 ;
        nodeM2 = contact_elements[icontact].nodeM2 ;
        shapeM2 = contact_elements[icontact].shapeM2 ;
        nodeM3 = contact_elements[icontact].nodeM3 ;
        shapeM3 = contact_elements[icontact].shapeM3 ;
        Px = xnorm * Pn + xtan * Pt ;
        Py = ynorm * Pn + ytan * Pt ;
        if (length_evolution == "Fixed")
        {
            Fsx = Px * initial_length ;
            Fsy = Py * initial_length ;
            W *= initial_length * Deltat ;
        }
        else if (length_evolution == "Evolutive")
        {
            Fsx = Px * length ;
            Fsy = Py * length ;
            W *= length * Deltat ;
        }
        contact_elements[icontact].gapn = gapn ;
        contact_elements[icontact].gapt = gapt ;
        contact_elements[icontact].length = length ;
        contact_elements[icontact].fx = Fsx ;
        contact_elements[icontact].fy = Fsy ;

        contact_work += 0.5 * W ;
        Bodies[bodyM].contact_work += 0.5 * W ; // Possible writing conflict !


        //if (gapn>1.e-6 & abs(Pn)>0.) cout << index << ' ' << contact_elements[icontact].borderS << ' ' << contact_elements[icontact].border_nodeS << ' ' << contact_elements[icontact].bodyM << ' ' << gapn << ' ' << gapt << ' ' << Pn << ' ' << Pt << ' ' << Px << ' ' << Py << endl ;
        //borders[borderS].x_contact_pressure[border_nodeS] += Px ;
        //borders[borderS].y_contact_pressure[border_nodeS] += Py ;

        if (type=="deformable")
        {
            number_influencing_nodes = nodes[nodeS].number_influencing_nodes ;
            influencing_nodes = nodes[nodeS].influencing_nodes ;
            shape_functions = nodes[nodeS].shape_functions ;
            for (int j(0) ; j<number_influencing_nodes ; j++)
            {
                nodes[influencing_nodes[j]].x_contact_force += Fsx * shape_functions[j] ;
                nodes[influencing_nodes[j]].y_contact_force += Fsy * shape_functions[j] ;
            }
        }
        else if (type=="rigid")
        {
            dx = nodes[nodeS].x_current - x_current ;
            dy = nodes[nodeS].y_current - y_current ;
            //cout << endl ;
            //cout << "slave " << index << ' ' << dx << ' ' << dy << ' ' << sqrt(dx*dx+dy*dy) << endl ;
            x_contact_force += Fsx ;
            y_contact_force += Fsy ;
            r_contact_force += -Fsx * dy + Fsy * dx ;
        }

//        for (int j(0) ; j<Bodies[bodyM].nb_neighbours ; j++)
//        {
//            if (Bodies[bodyM].neighbours[j][0] == index)
//            {
//                neig = j ;
//                break ;
//            }
//        }

        if (Bodies[bodyM].type=="deformable")
        {
            number_influencing_nodes = Bodies[bodyM].nodes[nodeM0].number_influencing_nodes ;
            influencing_nodes = Bodies[bodyM].nodes[nodeM0].influencing_nodes ;
            shape_functions = Bodies[bodyM].nodes[nodeM0].shape_functions ;
            for (int j(0) ; j<number_influencing_nodes ; j++)
            {
                //Bodies[bodyM].nodes[influencing_nodes[j]].x_master_contact_forces[neig] -= Fsx * shapeM0 * shape_functions[j] ;
                //Bodies[bodyM].nodes[influencing_nodes[j]].y_master_contact_forces[neig] -= Fsy * shapeM0 * shape_functions[j] ;
                contact_forces_to_send.push_back({ -Fsx * shapeM0 * shape_functions[j], -Fsy * shapeM0 * shape_functions[j] }) ;
                contact_forces_to_send_to.push_back({ bodyM, influencing_nodes[j] }) ;
            }
            number_influencing_nodes = Bodies[bodyM].nodes[nodeM1].number_influencing_nodes ;
            influencing_nodes = Bodies[bodyM].nodes[nodeM1].influencing_nodes ;
            shape_functions = Bodies[bodyM].nodes[nodeM1].shape_functions ;
            for (int j(0) ; j<number_influencing_nodes ; j++)
            {
                //Bodies[bodyM].nodes[influencing_nodes[j]].x_master_contact_forces[neig] -= Fsx * shapeM1 * shape_functions[j] ;
                //Bodies[bodyM].nodes[influencing_nodes[j]].y_master_contact_forces[neig] -= Fsy * shapeM1 * shape_functions[j] ;
                contact_forces_to_send.push_back({ -Fsx * shapeM1 * shape_functions[j], -Fsy * shapeM1 * shape_functions[j] }) ;
                contact_forces_to_send_to.push_back({ bodyM, influencing_nodes[j] }) ;
            }
            number_influencing_nodes = Bodies[bodyM].nodes[nodeM2].number_influencing_nodes ;
            influencing_nodes = Bodies[bodyM].nodes[nodeM2].influencing_nodes ;
            shape_functions = Bodies[bodyM].nodes[nodeM2].shape_functions ;
            for (int j(0) ; j<number_influencing_nodes ; j++)
            {
                //Bodies[bodyM].nodes[influencing_nodes[j]].x_master_contact_forces[neig] -= Fsx * shapeM2 * shape_functions[j] ;
                //Bodies[bodyM].nodes[influencing_nodes[j]].y_master_contact_forces[neig] -= Fsy * shapeM2 * shape_functions[j] ;
                contact_forces_to_send.push_back({ -Fsx * shapeM2 * shape_functions[j], -Fsy * shapeM2 * shape_functions[j] }) ;
                contact_forces_to_send_to.push_back({ bodyM, influencing_nodes[j] }) ;
            }
            number_influencing_nodes = Bodies[bodyM].nodes[nodeM3].number_influencing_nodes ;
            influencing_nodes = Bodies[bodyM].nodes[nodeM3].influencing_nodes ;
            shape_functions = Bodies[bodyM].nodes[nodeM3].shape_functions ;
            for (int j(0) ; j<number_influencing_nodes ; j++)
            {
                //Bodies[bodyM].nodes[influencing_nodes[j]].x_master_contact_forces[neig] -= Fsx * shapeM3 * shape_functions[j] ;
                //Bodies[bodyM].nodes[influencing_nodes[j]].y_master_contact_forces[neig] -= Fsy * shapeM3 * shape_functions[j] ;
                contact_forces_to_send.push_back({ -Fsx * shapeM3 * shape_functions[j], -Fsy * shapeM3 * shape_functions[j] }) ;
                contact_forces_to_send_to.push_back({ bodyM, influencing_nodes[j] }) ;
            }
        }
        else if (Bodies[bodyM].type=="rigid")
        {
            dx = nodes[nodeS].x_current - Bodies[bodyM].x_current - shiftM * period ;
            dy = nodes[nodeS].y_current - Bodies[bodyM].y_current ;
            contact_forces_to_send.push_back({ -Fsx, -Fsy, Fsx * dy - Fsy * dx }) ;
            contact_forces_to_send_to.push_back({ bodyM, -1 }) ;
            //Bodies[bodyM].x_master_contact_forces[neig] -= Fsx ;
            //Bodies[bodyM].y_master_contact_forces[neig] -= Fsy ;
            //Bodies[bodyM].r_master_contact_forces[neig] -= -Fsx * dy + Fsy * dx ;
            //cout << "master " << bodyM << ' ' << dx << ' ' << dy << ' ' << sqrt(dx*dx+dy*dy) << endl ;
        }
    }
}


//********************************************//
//** SEND CONTACT FORCES *********************//
//********************************************//

void Body::Send_contact_forces(vector<Body>& Bodies)
{
    int b, n ;
    for (int j=0 ; j<(int)contact_forces_to_send_to.size() ; j++)
    {
        b = contact_forces_to_send_to[j][0] ;
        n = contact_forces_to_send_to[j][1] ;
        if (n == -1)
        {
            Bodies[b].x_contact_force += contact_forces_to_send[j][0] ;
            Bodies[b].y_contact_force += contact_forces_to_send[j][1] ;
            Bodies[b].r_contact_force += contact_forces_to_send[j][2] ;
        }
        else
        {
            Bodies[b].nodes[n].x_contact_force += contact_forces_to_send[j][0] ;
            Bodies[b].nodes[n].y_contact_force += contact_forces_to_send[j][1] ;
        }
    }
}



//********************************************//
//** SUM UP FORCES ***************************//
//********************************************//

void Body::Sum_up_forces()
{
    if (type=="deformable")
    {
        for (int i(0) ; i<nb_nodes ; i++)
        {
            nodes[i].Sum_up_forces() ;
        }
    }
    else if (type=="rigid")
    {
        //for (int i(0) ; i<nb_neighbours ; i++)
        //{
        //    x_contact_force += x_master_contact_forces[i] ;
        //    y_contact_force += y_master_contact_forces[i] ;
        //    r_contact_force += r_master_contact_forces[i] ;
        //    x_master_contact_forces[i] = 0. ;
        //    y_master_contact_forces[i] = 0. ;
        //    r_master_contact_forces[i] = 0. ;
        //}
        x_force = x_contact_force + x_body_force + x_dirichlet_force + x_neumann_force + x_damping_force ;
        y_force = y_contact_force + y_body_force + y_dirichlet_force + y_neumann_force + y_damping_force ;
        r_force = r_contact_force + r_body_force + r_dirichlet_force + r_neumann_force + r_damping_force ;
        unbalanced = sqrt(x_force * x_force + y_force * y_force) ;/*/ (sqrt(x_contact_force * x_contact_force + y_contact_force * y_contact_force)
                                                                  + sqrt(x_body_force * x_body_force + y_body_force * y_body_force)
                                                                  + sqrt(x_dirichlet_force * x_dirichlet_force + y_dirichlet_force * y_dirichlet_force)
                                                                  + sqrt(x_neumann_force * x_neumann_force + y_neumann_force * y_neumann_force)
                                                                  + sqrt(x_damping_force * x_damping_force + y_damping_force * y_damping_force)) ;*/
        for (int i(0) ; i<nb_nodes ; i++)
        {
            nodes[i].x_contact_force = x_contact_force ;
            nodes[i].y_contact_force = y_contact_force ;
            nodes[i].x_body_force = x_body_force ;
            nodes[i].y_body_force = y_body_force ;
            nodes[i].x_dirichlet_force = x_dirichlet_force ;
            nodes[i].y_dirichlet_force = y_dirichlet_force ;
            nodes[i].x_neumann_force = x_neumann_force ;
            nodes[i].y_neumann_force = y_neumann_force ;
            nodes[i].x_damping_force = x_damping_force ;
            nodes[i].y_damping_force = y_damping_force ;
            nodes[i].x_force = x_force ;
            nodes[i].y_force = y_force ;
        }
    }
}



//********************************************//
//** APPLY NEWTON ****************************//
//********************************************//

void Body::Apply_Newton()
{
    if (type=="deformable")
    {
        for (int i(0) ; i<nb_nodes ; i++)
        {
            nodes[i].Apply_Newton() ;
        }
    }
    else if (type=="rigid")
    {
        x_acceleration = x_force * inverse_mass / factor_mass_scaling ;
        y_acceleration = y_force * inverse_mass / factor_mass_scaling ;
        r_acceleration = r_force * inverse_inertia / factor_mass_scaling ;
    }
    for (int i(0) ; i<nb_nodes ; i++)
    {
        nodes[i].x_acceleration = x_acceleration ;
        nodes[i].y_acceleration = y_acceleration ;
    }
}



//********************************************//
//** APPLY EULER ****************************//
//********************************************//

void Body::Apply_Euler(double Deltat)
{
    if (type=="deformable")
    {
        double vx, vy ;
        for (int i(0) ; i<nb_nodes ; i++)
        {
            vx = nodes[i].x_velocity_parameter + 0.5 * Deltat * nodes[i].x_acceleration_parameter ;
            vy = nodes[i].y_velocity_parameter + 0.5 * Deltat * nodes[i].y_acceleration_parameter ;
            internal_work += Deltat * ( vx * nodes[i].x_internal_force + vy * nodes[i].y_internal_force );
            //contact_work += Deltat * ( vx * nodes[i].x_contact_force + vy * nodes[i].y_contact_force ) ;
            body_work += Deltat * ( vx * nodes[i].x_body_force + vy * nodes[i].y_body_force ) ;
            dirichlet_work += Deltat * ( vx * nodes[i].x_dirichlet_force + vy * nodes[i].y_dirichlet_force ) ;
            neumann_work += Deltat * ( vx * nodes[i].x_neumann_force + vy * nodes[i].y_neumann_force ) ;
            damping_work += Deltat * ( vx * nodes[i].x_damping_force + vy * nodes[i].y_damping_force ) ;
            alid_work += Deltat * ( vx * nodes[i].x_alid_force + vy * nodes[i].y_alid_force ) ;
            nodes[i].Apply_Euler(Deltat) ;
        }
        temperature = -(contact_work + damping_work) / (mass * heat_capacity) ;
    }
    else if (type=="rigid")
    {
        double vx, vy, vr ;
        vx = x_velocity + 0.5 * Deltat * x_acceleration ;
        vy = y_velocity + 0.5 * Deltat * y_acceleration ;
        vr = r_velocity + 0.5 * Deltat * r_acceleration ;
        //contact_work += Deltat * ( vx * x_contact_force + vy * y_contact_force + vr * r_contact_force ) ;
        body_work += Deltat * ( vx * x_body_force + vy * y_body_force + vr * r_body_force ) ;
        dirichlet_work += Deltat * ( vx * x_dirichlet_force + vy * y_dirichlet_force + vr * r_dirichlet_force ) ;
        neumann_work += Deltat * ( vx * x_neumann_force + vy * y_neumann_force + vr * r_neumann_force ) ;
        damping_work += Deltat * ( vx * x_damping_force + vy * y_damping_force + vr * r_damping_force ) ;
        temperature = -(contact_work + damping_work) / (mass * heat_capacity) ;
        if (drivendof[0] == 0.)
        {
            x_velocity += Deltat * x_acceleration ;
            x_displacement = x_displacement + Deltat * x_velocity ;
        }
        else if (drivendof[0] == 1.)
        {
            x_acceleration = drivendof[3] ;//x_acceleration = ( x_velocity - ( drivendof[1] - x_displacement ) / Deltat ) / Deltat ;
            x_velocity = drivendof[2] ;//x_velocity = ( drivendof[1] - x_displacement ) / Deltat ;
            x_displacement = drivendof[1] ;
            x_dirichlet_force = mass * x_acceleration - x_force ;
        }
        if (drivendof[4] == 0.)
        {
            y_velocity += Deltat * y_acceleration ;
            y_displacement = y_displacement + Deltat * y_velocity ;
        }
        else if (drivendof[4] == 1.)
        {
            y_acceleration = drivendof[7] ;//y_acceleration = ( y_velocity - ( drivendof[5] - y_displacement ) / Deltat ) / Deltat ;
            y_velocity = drivendof[6] ;//y_velocity = ( drivendof[5] - y_displacement ) / Deltat ;
            y_displacement = drivendof[5] ;
            y_dirichlet_force = mass * y_acceleration - y_force ;
        }
        if (drivendof[0] == 0. && drivendof[4] == 0. && periodicity != "Periodic")
        {
            r_velocity += Deltat * r_acceleration ;
            r_displacement = r_displacement + Deltat * r_velocity ;
        }
        else
        {
            r_velocity = 0. ;
            r_displacement = 0. ;
        }
    }
}



//********************************************//
//** APPLY EULER TEMPORARY *******************//
//********************************************//

void Body::Apply_Euler_temporary(double Deltat)
{
    if (type=="deformable")
    {
        for (int i(0) ; i<nb_nodes ; i++)
        {
            nodes[i].Apply_Euler_temporary(Deltat) ;
        }
    }
    else if (type=="rigid")
    {
        if (drivendof[0] == 0.)
        {
            x_velocity_temporary = x_velocity + Deltat * x_acceleration ;
            x_displacement_temporary = x_displacement + Deltat * x_velocity_temporary ;
        }
        else if (drivendof[0] == 1.)
        {
            x_velocity_temporary = drivendof[2] ;//x_velocity_temporary = ( drivendof[1] - x_displacement ) / Deltat ;
            x_displacement_temporary = drivendof[1] ;
        }
        if (drivendof[4] == 0.)
        {
            y_velocity_temporary = y_velocity + Deltat * y_acceleration ;
            y_displacement_temporary = y_displacement + Deltat * y_velocity_temporary ;
        }
        else if (drivendof[4] == 1.)
        {
            y_velocity_temporary = drivendof[6] ;//y_velocity_temporary = ( drivendof[5] - y_displacement ) / Deltat ;
            y_displacement_temporary = drivendof[5] ;
        }
        if (drivendof[0] == 0. && drivendof[4] == 0. && periodicity != "Periodic")
        {
            r_velocity_temporary = r_velocity + Deltat * r_acceleration ;
            r_displacement_temporary = r_displacement + Deltat * r_velocity_temporary ;
        }
        else
        {
            r_velocity_temporary = 0. ;
            r_displacement_temporary = 0. ;
        }
    }
}



//********************************************//
//** UPDATE KINEMATICS ****************************//
//********************************************//

void Body::Update_kinematics()
{
    double dx, dy, dr, vx, vy, ax, ay ;
    if (type=="deformable")
    {
        for (int i(0) ; i<nb_nodes ; i++)
        {
            dx = 0. ;
            dy = 0. ;
            vx = 0. ;
            vy = 0. ;
            ax = 0. ;
            ay = 0. ;
            for (int j(0) ; j<nodes[i].number_influencing_nodes ; j++)
            {
                dx += nodes[i].shape_functions[j] * nodes[nodes[i].influencing_nodes[j]].x_displacement_parameter ;
                dy += nodes[i].shape_functions[j] * nodes[nodes[i].influencing_nodes[j]].y_displacement_parameter ;
                vx += nodes[i].shape_functions[j] * nodes[nodes[i].influencing_nodes[j]].x_velocity_parameter ;
                vy += nodes[i].shape_functions[j] * nodes[nodes[i].influencing_nodes[j]].y_velocity_parameter ;
                ax += nodes[i].shape_functions[j] * nodes[nodes[i].influencing_nodes[j]].x_acceleration_parameter ;
                ay += nodes[i].shape_functions[j] * nodes[nodes[i].influencing_nodes[j]].y_acceleration_parameter ;
            }
            //internal_work += ( dx - nodes[i].x_displacement ) * nodes[i].x_internal_force + ( dy - nodes[i].y_displacement ) * nodes[i].y_internal_force ;
            //contact_work += ( dx - nodes[i].x_displacement ) * nodes[i].x_contact_force + ( dy - nodes[i].y_displacement ) * nodes[i].y_contact_force ;
            //body_work += ( dx - nodes[i].x_displacement ) * nodes[i].x_body_force + ( dy - nodes[i].y_displacement ) * nodes[i].y_body_force ;
            //dirichlet_work += ( dx - nodes[i].x_displacement ) * nodes[i].x_dirichlet_force + ( dy - nodes[i].y_displacement ) * nodes[i].y_dirichlet_force ;
            //neumann_work += ( dx - nodes[i].x_displacement ) * nodes[i].x_neumann_force + ( dy - nodes[i].y_displacement ) * nodes[i].y_neumann_force ;
            //damping_work += ( dx - nodes[i].x_displacement ) * nodes[i].x_damping_force + ( dy - nodes[i].y_displacement ) * nodes[i].y_damping_force ;
            //alid_work += ( dx - nodes[i].x_displacement ) * nodes[i].x_alid_force + ( dy - nodes[i].y_displacement ) * nodes[i].y_alid_force ;
            // ENERGIE SUR LES DISPLACEMENT PARAMETERS, NON ?
            nodes[i].x_displacement = dx ;
            nodes[i].y_displacement = dy ;
            nodes[i].x_velocity = vx ;
            nodes[i].y_velocity = vy ;
            nodes[i].x_acceleration = ax ;
            nodes[i].y_acceleration = ay ;
        }
    }
    else if (type=="rigid")
    {
        for (int i(0) ; i<nb_nodes ; i++)
        {
            dx = nodes[i].x_initial - x_initial ;
            dy = nodes[i].y_initial - y_initial ;
            dr = r_displacement - r_initial ;
            nodes[i].x_displacement = x_displacement + dx * cos( dr ) - dy * sin( dr ) - dx ;
            nodes[i].y_displacement = y_displacement + dx * sin( dr ) + dy * cos( dr ) - dy ;
            nodes[i].x_velocity = x_velocity ;
            nodes[i].y_velocity = y_velocity ;
        }
    }
}



//********************************************//
//** UPDATE KINEMATICS TEMPORARY *************//
//********************************************//

void Body::Update_kinematics_temporary()
{
    double dx, dy, dr; // vx, vy, ax, ay ;
    if (type=="deformable")
    {
        for (int i(0) ; i<nb_nodes ; i++)
        {
            dx = 0. ;
            dy = 0. ;
            for (int j(0) ; j<nodes[i].number_influencing_nodes ; j++)
            {
                dx += nodes[i].shape_functions[j] * nodes[nodes[i].influencing_nodes[j]].x_displacement_parameter_temporary ;
                dy += nodes[i].shape_functions[j] * nodes[nodes[i].influencing_nodes[j]].y_displacement_parameter_temporary ;
            }
            nodes[i].x_displacement_temporary = dx ;
            nodes[i].y_displacement_temporary = dy ;
        }
    }
    else if (type=="rigid")
    {
        for (int i(0) ; i<nb_nodes ; i++)
        {
            dx = nodes[i].x_initial - x_initial ;
            dy = nodes[i].y_initial - y_initial ;
            dr = r_displacement_temporary - r_initial ;
            nodes[i].x_displacement_temporary = x_displacement_temporary + dx * cos( dr ) - dy * sin( dr ) - dx ;
            nodes[i].y_displacement_temporary = y_displacement_temporary + dx * sin( dr ) + dy * cos( dr ) - dy ;
        }
    }
}



//********************************************//
//** UPDATE CURRENT POSITIONS ****************//
//********************************************//

void Body::Update_current_positions()
{
    if (type=="deformable")
    {
        for (int i(0) ; i<nb_nodes ; i++)
        {
            nodes[i].Update_current_positions() ;
        }
    }
    else if (type=="rigid")
    {
        double dx, dy, dr ;
        x_current = x_initial + x_displacement ;
        y_current = y_initial + y_displacement ;
        r_current = r_initial + r_displacement ;
        for (int i(0) ; i<nb_nodes ; i++)
        {
            dx = nodes[i].x_initial - x_initial ;
            dy = nodes[i].y_initial - y_initial ;
            dr = r_displacement - r_initial ;
            nodes[i].x_current = x_current + dx * cos( dr ) - dy * sin( dr ) ;
            nodes[i].y_current = y_current + dx * sin( dr ) + dy * cos( dr ) ;
            //nodes[i].x_displacement = nodes[i].x_current - nodes[i].x_initial ;
            //nodes[i].y_displacement = nodes[i].y_current - nodes[i].y_initial ;
            dx = nodes[i].x_current - x_current ;
            dy = nodes[i].y_current - y_current ;
            nodes[i].x_velocity = -dy * r_velocity + x_velocity ;
            nodes[i].y_velocity = dx * r_velocity + y_velocity ;
        }
    }
}



//********************************************//
//** COMPUTE NODAL STRESSES ******************//
//********************************************//

void Body::Compute_nodal_stresses( int Nb_materials, vector<Material>& Materials )
{
    if (type=="deformable")
    {
        string material_type ;
        vector<double> parameters ;
        double SigmaA, SigmaB ;
        /*
        for (int i(0) ; i<Nb_materials ; i++)
        {
            if (Materials[i].name == material_name)
            {
                material_type=Materials[i].type ;
                parameters=Materials[i].parameters ;
                break ;
            }
        }
        */
        material_type=Materials[material_index].type ;
        parameters=Materials[material_index].parameters ;
        int number_influencing_nodes ;
        vector<int> influencing_nodes ;
        vector<double> shape_xderiv, shape_yderiv ;
        double J, Sigmaxx(0), Sigmayy(0), Sigmaxy(0), Sigmazz(0), Mu, Lambda, Kappa, energy, Exx(0), Eyy(0), Exy(0), NormE(0) ;
        //vector<vector<double>> S={{0,0,0},{0,0,0},{0,0,0}};
        //vector<vector<double>> E={{0,0,0},{0,0,0},{0,0,0}};
        //vector<vector<double>> Eref={{0,0,0},{0,0,0},{0,0,0}};
        //vector<vector<double>> Fref={{1,0,0},{0,1,0},{0,0,1}};
        //vector<vector<double>> F ;
        double E11, E12, E22 ;
        double S11, S12, S22, S33 ;
        double F11, F12, F21, F22 ;
        for (int i(0) ; i<nb_nodes ; i++)
        {
            number_influencing_nodes = nodes[i].number_influencing_nodes ;
            influencing_nodes = nodes[i].influencing_nodes ;
            shape_xderiv = nodes[i].shape_xderiv ;
            shape_yderiv = nodes[i].shape_yderiv ;
            if (material_type == "ElasticLinear")
            {
                //E = Eref ;
                E11 = 0 ;
                E12 = 0 ;
                E22 = 0 ;
                for (int j(0) ; j<number_influencing_nodes ; j++)
                {
                    E11 += shape_xderiv[j] * nodes[influencing_nodes[j]].x_displacement_parameter ;
                    E12 += shape_xderiv[j] * nodes[influencing_nodes[j]].y_displacement_parameter + shape_yderiv[j] * nodes[influencing_nodes[j]].x_displacement_parameter ;
                    E22 += shape_yderiv[j] * nodes[influencing_nodes[j]].y_displacement_parameter ;
                }
                E12 = 0.5 * E12 ;
                //E[1][0] = E[0][1] ;
                Mu = parameters[3] ;
                Lambda = parameters[4] ;
                //Apply_Elastic_Linear( S , E , Mu , Lambda , J , energy ) ;
                Apply_Elastic_Linear( S11, S12, S22, S33, E11, E12, E22, Mu, Lambda, J, energy ) ;
                if (J<0.)
                {
                    S11 = 0. ;
                    S12 = 0. ;
                    S22 = 0. ;
                    S33 = 0. ;
                }
                Sigmaxx = S11 ;
                Sigmayy = S22 ;
                Sigmaxy = S12 ;
                Sigmazz = S33 ;
                Exx = E11 ;
                Eyy = E22 ;
                Exy = E12 ;
                NormE = pow ( Exx * Exx + Eyy * Eyy + 2. * Exy * Exy, 0.5 ) ;
            }
            else if (material_type == "NeoHookean")
            {
                //F = Fref ;
                F11 = 1. ;
                F12 = 0. ;
                F21 = 0. ;
                F22 = 1. ;
                for (int j(0) ; j<number_influencing_nodes ; j++)
                {
                    F11 += shape_xderiv[j] * nodes[influencing_nodes[j]].x_displacement_parameter ;
                    F12 += shape_yderiv[j] * nodes[influencing_nodes[j]].x_displacement_parameter ;
                    F21 += shape_xderiv[j] * nodes[influencing_nodes[j]].y_displacement_parameter ;
                    F22 += shape_yderiv[j] * nodes[influencing_nodes[j]].y_displacement_parameter ;
                }
                Mu = parameters[3] ;
                Kappa = parameters[4] ;
                Apply_NeoHookean( S11, S12, S22, S33, F11, F12, F21, F22, Mu, Kappa, J, energy ) ;
                if (J<0.)
                {
                    F11 = 0. ;
                    F12 = 0. ;
                    F21 = 0. ;
                    F22 = 0. ;
                    S11 = 0. ;
                    S12 = 0. ;
                    S22 = 0. ;
                    S33 = 0. ;
                }
                //Sigmaxx = F[0][0] * ( S[0][0] * F[0][0] + S[0][1] * F[0][1] ) + F[0][1] * ( S[0][1] * F[0][0] + S[1][1] * F[0][1] ) ;
                //Sigmayy = F[1][0] * ( S[0][0] * F[1][0] + S[0][1] * F[1][1] ) + F[1][1] * ( S[0][1] * F[1][0] + S[1][1] * F[1][1] ) ;
                //Sigmaxy = F[1][0] * ( S[0][0] * F[0][0] + S[0][1] * F[0][1] ) + F[1][1] * ( S[0][1] * F[0][0] + S[1][1] * F[0][1] ) ;
                //Sigmazz = S[2][2] ;
                //Exx = 0.5 * ( F[0][0] * F[0][0] + F[1][0] * F[1][0] - 1. ) ;
                //Eyy = 0.5 * ( F[0][1] * F[0][1] + F[1][1] * F[1][1] - 1. ) ;
                //Exy = 0.5 * ( F[0][1] * F[0][0] + F[1][1] * F[1][0] ) ;
                Sigmaxx = F11 * ( S11 * F11 + S12 * F12 ) + F12 * ( S12 * F11 + S22 * F12 ) ;
                Sigmayy = F21 * ( S11 * F21 + S12 * F22 ) + F22 * ( S12 * F21 + S22 * F22 ) ;
                Sigmaxy = F21 * ( S11 * F11 + S12 * F12 ) + F22 * ( S12 * F11 + S22 * F12 ) ;
                Sigmazz = S33 ;
                Exx = 0.5 * ( F11 * F11 + F21 * F21 - 1. ) ;
                Eyy = 0.5 * ( F12 * F12 + F22 * F22 - 1. ) ;
                Exy = 0.5 * ( F12 * F11 + F22 * F21 ) ;
                NormE = pow ( Exx * Exx + Eyy * Eyy + 2. * Exy * Exy, 0.5 ) ;
            }
            //
            // OTHER MATERIALS ??
            //
            SigmaA = 0.5 * ( Sigmaxx + Sigmayy + sqrt( ( Sigmaxx + Sigmayy ) * ( Sigmaxx + Sigmayy ) - 4. * ( Sigmaxx * Sigmayy - Sigmaxy * Sigmaxy ) ) ) ;
            SigmaB = 0.5 * ( Sigmaxx + Sigmayy - sqrt( ( Sigmaxx + Sigmayy ) * ( Sigmaxx + Sigmayy ) - 4. * ( Sigmaxx * Sigmayy - Sigmaxy * Sigmaxy ) ) ) ;
            if ( (SigmaA>SigmaB) && (SigmaA>Sigmazz) )
            {
                nodes[i].SigmaI = SigmaA ;
                if (SigmaB>Sigmazz)
                {
                    nodes[i].SigmaII = SigmaB ;
                    nodes[i].SigmaIII = Sigmazz ;
                }
                else
                {
                    nodes[i].SigmaII = Sigmazz ;
                    nodes[i].SigmaIII = SigmaB ;
                }
            }
            else if ( (SigmaB>SigmaA) && (SigmaB>Sigmazz) )
            {
                nodes[i].SigmaI = SigmaB ;
                if (SigmaA>Sigmazz)
                {
                    nodes[i].SigmaII = SigmaA ;
                    nodes[i].SigmaIII = Sigmazz ;
                }
                else
                {
                    nodes[i].SigmaII = Sigmazz ;
                    nodes[i].SigmaIII = SigmaA ;
                }
            }
            else
            {
                nodes[i].SigmaI = Sigmazz ;
                if (SigmaB>SigmaA)
                {
                    nodes[i].SigmaII = SigmaB ;
                    nodes[i].SigmaIII = SigmaA ;
                }
                else
                {
                    nodes[i].SigmaII = SigmaA ;
                    nodes[i].SigmaIII = SigmaB ;
                }
            }
            nodes[i].jacobian = J ;
            nodes[i].Sxx = S11 ;
            nodes[i].Syy = S22 ;
            nodes[i].Sxy = S12 ;
            nodes[i].Szz = S33 ;
            nodes[i].Sigmaxx = Sigmaxx ;
            nodes[i].Sigmayy = Sigmayy ;
            nodes[i].Sigmaxy = Sigmaxy ;
            nodes[i].Sigmazz = Sigmazz ;
            nodes[i].SigmaTresca = nodes[i].SigmaI - nodes[i].SigmaIII ;
            nodes[i].SigmaVM = sqrt( (nodes[i].Sigmaxx - nodes[i].Sigmayy) * (nodes[i].Sigmaxx - nodes[i].Sigmayy) +
                                     (nodes[i].Sigmaxx - nodes[i].Sigmazz) * (nodes[i].Sigmaxx - nodes[i].Sigmazz) +
                                     (nodes[i].Sigmazz - nodes[i].Sigmayy) * (nodes[i].Sigmazz - nodes[i].Sigmayy) +
                                     6. * nodes[i].Sigmaxy * nodes[i].Sigmaxy ) * 0.70710678118655 ;
            nodes[i].SigmaSph = 0.333333333 * ( Sigmaxx + Sigmayy + Sigmazz ) ;
            nodes[i].Exx = Exx ;
            nodes[i].Eyy = Eyy ;
            nodes[i].Exy = Exy ;
            nodes[i].NormE = NormE ;
        }
    }
    else if (type=="rigid")
    {
        for (int i(0) ; i<nb_nodes ; i++)
        {
            nodes[i].SigmaI = 0. ;
            nodes[i].SigmaII = 0. ;
            nodes[i].SigmaIII = 0. ;
            nodes[i].jacobian = 1. ;
            nodes[i].Sxx = 0. ;
            nodes[i].Syy = 0. ;
            nodes[i].Sxy = 0. ;
            nodes[i].Szz = 0. ;
            nodes[i].Sigmaxx = 0. ;
            nodes[i].Sigmayy = 0. ;
            nodes[i].Sigmaxy = 0. ;
            nodes[i].Sigmazz = 0. ;
            nodes[i].SigmaTresca = 0. ;
            nodes[i].SigmaVM = 0. ;
            nodes[i].SigmaSph = 0. ;
            nodes[i].Exx = 0. ;
            nodes[i].Eyy = 0. ;
            nodes[i].Exy = 0. ;
            nodes[i].NormE = 0. ;
        }
    }
}




//********************************************//
//** COMPUTE ERROR ***************************//
//********************************************//

void Body::Compute_error()
{
    total_error = 0. ;
    if (type=="deformable")
    {
        max_error = 0. ;
        node_for_max_error = 0 ;
        for (int i(0) ; i<nb_nodes ; i++)
            nodes[i].Compute_error(total_error, max_error, node_for_max_error, drivendof[0], drivendof[4], i) ;
        //if (std::isnan(total_error)) status = "inactive" ;
    }
    else if (type=="rigid")
    {
        if (drivendof[0] != 1.)
            total_error += ( x_displacement - x_displacement_temporary ) * ( x_displacement - x_displacement_temporary ) ;
        if (drivendof[4] != 1.)
            total_error += ( y_displacement - y_displacement_temporary ) * ( y_displacement - y_displacement_temporary ) ;
        total_error = 0.1 * sqrt( total_error ) / nodal_distance ;
        if (total_error==0.)
            total_error = 1.e-16 ;
        max_error = total_error ;
        node_for_max_error = -1 ;
        for (int i(0) ; i<nb_nodes ; i++)
            nodes[i].error_norm = total_error ;
    }
}

//********************************************//
//** COMPUTE MASS SCALING ********************//
//********************************************//

void Body::Compute_mass_scaling(double Target_error, double Inv_Target_error, double Control_parameter_mass_scaling, double Max_mass_scaling,
                                double Error_factor_mass_scaling, double Accepted_ratio, double Decrease_factor_mass_scaling)
{
    mass_mass_scaling = 0. ;

    if (type == "deformable")
    {
        max_factor_mass_scaling = 0. ;
        for (int i(0) ; i<nb_nodes ; i++)
            nodes[i].Compute_mass_scaling(Target_error, Inv_Target_error, Control_parameter_mass_scaling, Max_mass_scaling, Error_factor_mass_scaling,
                                          Accepted_ratio, Decrease_factor_mass_scaling, max_factor_mass_scaling, index, i, mass_mass_scaling) ;

    }
    else if (type == "rigid")
    {

        //ADDITIVE
        if (total_error >= Target_error * Error_factor_mass_scaling)
        {
            delta_factor_mass_scaling = pow((total_error - Target_error * Error_factor_mass_scaling), Control_parameter_mass_scaling) ;
            factor_mass_scaling += delta_factor_mass_scaling ;
            if (factor_mass_scaling >= Max_mass_scaling)
                factor_mass_scaling = Max_mass_scaling ;
        }
        else
        {
            if(factor_mass_scaling != 1)
                factor_mass_scaling *= Decrease_factor_mass_scaling ;
            if(factor_mass_scaling < 1)
                factor_mass_scaling = 1;
        }

        mass_mass_scaling = mass * factor_mass_scaling ;
    }

}


//********************************************//
//** STORE ***********************************//
//********************************************//

void Body::Store()
{
    stored_nodes = nodes ;
    stored_borders = borders ;
    stored_contact_elements = contact_elements ;
    stored_internal_work = internal_work ;
    stored_contact_work = contact_work ;
    stored_body_work = body_work ;
    stored_dirichlet_work = dirichlet_work ;
    stored_neumann_work = neumann_work ;
    stored_damping_work = damping_work ;
    stored_alid_work = alid_work ;
    stored_delta_factor_mass_scaling = delta_factor_mass_scaling ;
    if (type=="rigid")
    {
        stored_x_current = x_current ;
        stored_y_current = y_current ;
        stored_r_current = r_current ;
        stored_x_displacement = x_displacement ;
        stored_y_displacement = y_displacement ;
        stored_r_displacement = r_displacement ;
        stored_x_velocity = x_velocity ;
        stored_y_velocity = y_velocity ;
        stored_r_velocity = r_velocity ;
    }
}



//********************************************//
//** RESTORE *********************************//
//********************************************//

void Body::Restore()
{
    nodes = stored_nodes ;
    borders = stored_borders ;
    contact_elements = stored_contact_elements ;
    internal_work = stored_internal_work ;
    contact_work = stored_contact_work ;
    body_work = stored_body_work ;
    dirichlet_work = stored_dirichlet_work ;
    neumann_work = stored_neumann_work ;
    damping_work = stored_damping_work ;
    alid_work = stored_alid_work ;
    delta_factor_mass_scaling = stored_delta_factor_mass_scaling ;
    if (type=="rigid")
    {
        x_current = stored_x_current ;
        y_current = stored_y_current ;
        r_current = stored_r_current ;
        x_displacement = stored_x_displacement ;
        y_displacement = stored_y_displacement ;
        r_displacement = stored_r_displacement ;
        x_velocity = stored_x_velocity ;
        y_velocity = stored_y_velocity ;
        r_velocity = stored_r_velocity ;
    }
}



//********************************************//
//** UPDATE DAMAGE ***************************//
//********************************************//

void Body::Update_damage()
{
    double total_nodes = 0. ;
    double d ;
    double damage_previous = damage ;
    damage = 0. ;
    for (int i(0) ; i<nb_borders ; i++ )
    {
        total_nodes += borders[i].number_border_nodes ;
        for (int j(0) ; j<borders[i].number_border_nodes ; j++)
        {
            d = 1. ;
            for (int k(0) ; k<nb_contact_elements ; k++)
            {
                if (contact_elements[k].borderS == i && contact_elements[k].border_nodeS == j )
                {
                    if (contact_elements[k].nb_internal > 0)
                    {
                        if (contact_elements[k].internal[0] < d )
                            d = contact_elements[k].internal[0] ;
                    }
                }
            }
            damage += d ;
            //cout << damage << endl ;
        }
    }
    damage = damage / total_nodes ;
    //cout << damage << endl ;
    if (damage < damage_previous)
        damage = damage_previous ;
    //cout << damage << endl ;
    //cout << ' ' << endl ;
}



//********************************************//
//** UPDATE MATERIAL *************************//
//********************************************//

void Body::Update_material( int Nb_materials, vector<Material>& Materials, vector<int> flags )
{
    //cout << "Updating Material Properties of Body " << index << endl ;
    double Rho(0), Alpha(0), Beta(0), Mu(0), Lambda(0), Kappa, E(0), Nu(0) ;
    string material_type ;
    vector<double> parameters ;
    int dof0, dof1 ;
    /*
    for (int i(0) ; i<Nb_materials ; i++)
    {
        if (Materials[i].name == material_name)
        {
            material_type=Materials[i].type ;
            parameters=Materials[i].parameters ;
            break ;
        }
    }
    */
    material_type=Materials[material_index].type ;
    parameters=Materials[material_index].parameters ;
    if (material_type == "ElasticLinear")
    {
        Rho = parameters[0] ;
        Alpha = parameters[1] ;
        Beta = parameters[2] ;
        Mu = parameters[3] ;
        Lambda = parameters[4] ;
        E = Mu * (3. * Lambda + 2. * Mu) / (Lambda + Mu) ;
        Nu = Mu / (2. * Lambda + 2. * Mu) ;
    }
    else if (material_type == "NeoHookean")
    {
        Rho = parameters[0] ;
        Alpha = parameters[1] ;
        Beta = parameters[2] ;
        Mu = parameters[3] ;
        Kappa = parameters[4] ;
        E = 9. * Kappa * Mu / (3. * Kappa + Mu) ;
        Nu = (3. * Kappa - 2. * Mu) / (6. * Kappa + 2. * Mu) ;
    }
    //
    // OTHER MATERIALS ??
    //
    if (type=="deformable")
    {
        if (flags[4] == 1)
        {
            double dm ;
            mass = 0. ;
            for (int i(0) ; i<nb_nodes ; i++)
            {
                nodes[i].x_mass = 0. ;
                nodes[i].y_mass = 0. ;
            }
            for (int i(0) ; i<nb_gauss ; i++)
            {
                for (int j(0) ; j<gpoints[i].number_influencing_nodes ; j++)
                {
                    for (int k(0) ; k<gpoints[i].number_influencing_nodes ; k++)
                    {
                        dm = gpoints[i].weight * gpoints[i].jacobian * gpoints[i].shape_functions[j] * gpoints[i].shape_functions[k] ;
                        nodes[gpoints[i].influencing_nodes[j]].x_mass += dm ;
                        nodes[gpoints[i].influencing_nodes[j]].y_mass += dm ;
                    }
                }
            }
            for (int i(0) ; i<nb_nodes ; i++)
            {
                nodes[i].x_mass *= Rho ;
                nodes[i].y_mass *= Rho ;
                nodes[i].x_inverse_mass = 1. / nodes[i].x_mass ;
                nodes[i].y_inverse_mass = 1. / nodes[i].y_mass ;
                mass += sqrt(nodes[i].x_mass * nodes[i].x_mass + nodes[i].y_mass * nodes[i].y_mass);
                inverse_mass = 1/mass;
            }
        }
        if (flags[5] == 1)
        {
            vector<vector<double>> Elasticity_Matrix = {{1.-Nu,Nu,0},{Nu,1.-Nu,0},{0,0,1.-2.*Nu}};
            for (int i(0) ; i<3 ; i++)
                for (int j(0) ; j<3 ; j++)
                    Elasticity_Matrix[i][j] *= E / ((1.+Nu) * (1.-2.*Nu)) ;
            for (int i(0) ; i<nb_matrix ; i++)
                stiffness[i] = 0. ;
            for (int g(0) ; g<nb_gauss ; g++)
            {
                vector<double> Empty(2*gpoints[g].number_influencing_nodes) ;
                vector<vector<double>> Local_Stiffness_Matrix(2*gpoints[g].number_influencing_nodes) ;
                for (int j(0) ; j<2*gpoints[g].number_influencing_nodes ; j++)
                    Local_Stiffness_Matrix[j] = Empty ;
                for (int ii(0) ; ii<gpoints[g].number_influencing_nodes ; ii++)
                {
                    int i = 2 * ii ;
                    for (int jj(0) ; jj<gpoints[g].number_influencing_nodes ; jj++)
                    {
                        int j = 2 * jj ;
                        Local_Stiffness_Matrix[i][j] += gpoints[g].shape_xderiv[ii]*Elasticity_Matrix[0][0]*gpoints[g].shape_xderiv[jj];
                        Local_Stiffness_Matrix[i][j] += gpoints[g].shape_xderiv[ii]*Elasticity_Matrix[0][2]*gpoints[g].shape_yderiv[jj];
                        Local_Stiffness_Matrix[i][j] += gpoints[g].shape_yderiv[ii]*Elasticity_Matrix[2][0]*gpoints[g].shape_xderiv[jj];
                        Local_Stiffness_Matrix[i][j] += gpoints[g].shape_yderiv[ii]*Elasticity_Matrix[2][2]*gpoints[g].shape_yderiv[jj];
                    }
                    for (int jj(0) ; jj<gpoints[g].number_influencing_nodes ; jj++)
                    {
                        int j = 2 * jj + 1 ;
                        Local_Stiffness_Matrix[i][j] += gpoints[g].shape_xderiv[ii]*Elasticity_Matrix[0][1]*gpoints[g].shape_yderiv[jj];
                        Local_Stiffness_Matrix[i][j] += gpoints[g].shape_xderiv[ii]*Elasticity_Matrix[0][2]*gpoints[g].shape_xderiv[jj];
                        Local_Stiffness_Matrix[i][j] += gpoints[g].shape_yderiv[ii]*Elasticity_Matrix[2][1]*gpoints[g].shape_yderiv[jj];
                        Local_Stiffness_Matrix[i][j] += gpoints[g].shape_yderiv[ii]*Elasticity_Matrix[2][2]*gpoints[g].shape_xderiv[jj];
                    }
                }
                for (int ii(0) ; ii<gpoints[g].number_influencing_nodes ; ii++)
                {
                    int i = 2 * ii + 1 ;
                    for (int jj(0) ; jj<gpoints[g].number_influencing_nodes ; jj++)
                    {
                        int j = 2 * jj ;
                        Local_Stiffness_Matrix[i][j] += gpoints[g].shape_yderiv[ii]*Elasticity_Matrix[1][0]*gpoints[g].shape_xderiv[jj];
                        Local_Stiffness_Matrix[i][j] += gpoints[g].shape_yderiv[ii]*Elasticity_Matrix[1][2]*gpoints[g].shape_yderiv[jj];
                        Local_Stiffness_Matrix[i][j] += gpoints[g].shape_xderiv[ii]*Elasticity_Matrix[2][0]*gpoints[g].shape_xderiv[jj];
                        Local_Stiffness_Matrix[i][j] += gpoints[g].shape_xderiv[ii]*Elasticity_Matrix[2][2]*gpoints[g].shape_yderiv[jj];
                    }
                    for (int jj(0) ; jj<gpoints[g].number_influencing_nodes ; jj++)
                    {
                        int j = 2 * jj + 1 ;
                        Local_Stiffness_Matrix[i][j] += gpoints[g].shape_yderiv[ii]*Elasticity_Matrix[1][1]*gpoints[g].shape_yderiv[jj];
                        Local_Stiffness_Matrix[i][j] += gpoints[g].shape_yderiv[ii]*Elasticity_Matrix[1][2]*gpoints[g].shape_xderiv[jj];
                        Local_Stiffness_Matrix[i][j] += gpoints[g].shape_xderiv[ii]*Elasticity_Matrix[2][1]*gpoints[g].shape_yderiv[jj];
                        Local_Stiffness_Matrix[i][j] += gpoints[g].shape_xderiv[ii]*Elasticity_Matrix[2][2]*gpoints[g].shape_xderiv[jj];
                    }
                }

                for (int i(0) ; i < 2 * gpoints[g].number_influencing_nodes ; i++)
                {
                    dof0 = 2 * gpoints[g].influencing_nodes[i/2] + (i - 2 * (i/2)) ;
                    for (int j(0) ; j < 2 * gpoints[g].number_influencing_nodes ; j++)
                    {
                        dof1 = 2 * gpoints[g].influencing_nodes[j/2] + (j - 2 * (j/2)) ;
                        if (dof1<dof0)
                            continue ;
                        for (int n(0) ; n<nb_matrix ; n++)
                        {
                            if (matrix_coordinates[n][0] == dof0)
                            {
                                if (matrix_coordinates[n][1] == dof1)
                                {
                                    stiffness[n] += Local_Stiffness_Matrix[i][j] * gpoints[g].weight * gpoints[g].jacobian ;
                                    break ;
                                }
                            }
                        }
                    }
                }
            }
        }
        if (flags[6] == 1)
        {
            int i0 ;
            for (int i(0) ; i<nb_matrix ; i++)
            {
                if (matrix_coordinates[i][0] == matrix_coordinates[i][1])
                {
                    i0 = matrix_coordinates[i][0] / 2 ;
                    damping[i] = Alpha * stiffness[i] + Beta * nodes[i0].x_mass ;
                }
                else
                    damping[i] = Alpha * stiffness[i] ;
            }
        }
    }
    else if (type=="rigid")
    {
        mass = 0. ;
        inertia = 0. ;
        int n0, n1, n2 ;
        double x0, y0, x1, y1, x2, y2, xc, yc, surf ;
        for (int i(0) ; i<nb_cells ; i++)
        {
            n0 = triangulation[i][0] ;
            n1 = triangulation[i][1] ;
            n2 = triangulation[i][2] ;
            x0 = nodes[n0].x_initial ;
            y0 = nodes[n0].y_initial ;
            x1 = nodes[n1].x_initial ;
            y1 = nodes[n1].y_initial ;
            x2 = nodes[n2].x_initial ;
            y2 = nodes[n2].y_initial ;
            xc = ( x0 + x1 + x2 ) / 3. ;
            yc = ( y0 + y1 + y2 ) / 3. ;
            surf = Triangle_surface( x0, y0, x1, y1, x2, y2 ) ;
            mass += surf ;
            inertia += surf * ( (xc-x_initial)*(xc-x_initial) + (yc-y_initial)*(yc-y_initial) ) ;
            inertia += Triangle_inertia( x0, y0, x1, y1, x2, y2, xc, yc ) ;
        }
        mass *= Rho ;
        inertia *= Rho ;
        inverse_mass = 1. / mass ;
        inverse_inertia = 1. / inertia ;
    }
}

#endif
