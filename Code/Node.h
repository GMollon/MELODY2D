#ifndef DEF_NODE
#define DEF_NODE

// Header of the class

class Node
{
public :

    // Static attributes
    double x_initial ;
    double y_initial ;
    double domain_size ;
    double x_mass ;
    double y_mass ;
    double x_inverse_mass ;
    double y_inverse_mass ;
    int    number_influencing_nodes ;
    vector<int> influencing_nodes ;
    vector<double> shape_functions ;
    vector<double> shape_xderiv ;
    vector<double> shape_yderiv ;

    // Dynamic attributes
    double x_current ;
    double y_current ;
    double x_displacement ;
    double y_displacement ;
    double x_velocity ;
    double y_velocity ;
    double x_acceleration ;
    double y_acceleration ;
    double x_displacement_parameter ;
    double y_displacement_parameter ;
    double x_velocity_parameter ;
    double y_velocity_parameter ;
    double x_acceleration_parameter ;
    double y_acceleration_parameter ;
    double x_internal_force ;
    vector<double> x_regions_internal_forces ;
    double y_internal_force ;
    vector<double> y_regions_internal_forces ;
    double x_contact_force ;
    double y_contact_force ;
    vector<double> x_master_contact_forces ;
    vector<double> y_master_contact_forces ;
    double x_self_contact_force ;
    double y_self_contact_force ;
    double x_body_force ;
    double y_body_force ;
    double x_dirichlet_force ;
    double y_dirichlet_force ;
    double x_neumann_force ;
    double y_neumann_force ;
    double x_damping_force ;
    double y_damping_force ;
    double x_alid_force ;
    double y_alid_force ;
    double x_force ;
    double y_force ;
    double jacobian ;
    double Sxx ;
    double Syy ;
    double Sxy ;
    double Szz ;
    double Sigmaxx ;
    double Sigmayy ;
    double Sigmaxy ;
    double Sigmazz ;
    double SigmaI ;
    double SigmaII ;
    double SigmaIII ;
    double SigmaTresca ;
    double SigmaVM ;
    double SigmaSph ;
    double Exx ;
    double Eyy ;
    double Exy ;
    double NormE ;
    double x_displacement_parameter_temporary ;
    double y_displacement_parameter_temporary ;
    double x_displacement_temporary ;
    double y_displacement_temporary ;
    double inverse_distance_estimator ;
    double x_error ;
    double y_error ;
    double error_norm ;
    double alid_alpha ;
    double alid_beta ;

    // Mass scaling attributes
    double delta_factor_mass_scaling ;
    double factor_mass_scaling ;
    double x_factor_mass_scaling ;
    double y_factor_mass_scaling ;
    double delta_x_factor_mass_scaling ;
    double delta_y_factor_mass_scaling ;
    double max_factor_mass_scaling;
    double mass_mass_scaling;

    // Constructor and Destructor
    Node(double x, double y, double d, double mx, double my, double imx, double imy) ;
    ~Node() ;

    // Accessors

    // Modifiers

    // Methods
    void Sum_up_forces() ;
    void Apply_Newton() ;
    void Update_current_positions() ;
    void Apply_Euler(double dt) ;
    void Apply_Euler_temporary(double dt) ;
    void Compute_error(double& total_error, double& max_error, int& node_for_max_error, double& flagdrivenx, double& flagdriveny, int index) ;
    void Compute_mass_scaling(double Target_error, double Inv_Target_error, double Control_parameter_mass_scaling, double Max_mass_scaling,
                              double Error_factor_mass_scaling, double Accepted_ratio, double Decrease_factor_mass_scaling, double& max_factor_mass_scaling, int index, int num_node, double& total_mass_mass_scaling) ;
} ;



// Content of the class

// Constructor and Destructor
Node::Node (double x, double y, double d, double mx, double my, double imx, double imy)
{
    x_initial = x ;
    y_initial = y ;
    domain_size = d ;
    x_mass = mx ;
    y_mass = my ;
    x_inverse_mass = imx ;
    y_inverse_mass = imy ;

    number_influencing_nodes = 0 ;
    x_current = 0. ;
    y_current = 0. ;
    x_displacement = 0. ;
    y_displacement = 0. ;
    x_velocity = 0. ;
    y_velocity = 0. ;
    x_acceleration = 0. ;
    y_acceleration = 0. ;
    x_displacement_parameter = 0. ;
    y_displacement_parameter = 0. ;
    x_velocity_parameter = 0. ;
    y_velocity_parameter = 0. ;
    x_acceleration_parameter = 0. ;
    y_acceleration_parameter = 0. ;
    x_internal_force = 0. ;
    x_regions_internal_forces = {0.} ;
    y_internal_force = 0. ;
    y_regions_internal_forces = {0.} ;
    x_contact_force = 0. ;
    y_contact_force = 0. ;
    x_master_contact_forces = {0.} ;
    y_master_contact_forces = {0.} ;
    x_self_contact_force = 0. ;
    y_self_contact_force = 0. ;
    x_body_force = 0. ;
    y_body_force = 0. ;
    x_dirichlet_force = 0. ;
    y_dirichlet_force = 0. ;
    x_neumann_force = 0. ;
    y_neumann_force = 0. ;
    x_damping_force = 0. ;
    y_damping_force = 0. ;
    x_alid_force = 0. ;
    y_alid_force = 0. ;
    x_force = 0. ;
    y_force = 0. ;
    jacobian = 1. ;
    Sxx = 0. ;
    Syy = 0. ;
    Sxy = 0. ;
    Szz = 0. ;
    Sigmaxx = 0. ;
    Sigmayy = 0. ;
    Sigmaxy = 0. ;
    Sigmazz = 0. ;
    SigmaI = 0. ;
    SigmaII = 0. ;
    SigmaIII = 0. ;
    SigmaTresca = 0. ;
    SigmaVM = 0. ;
    x_displacement_parameter_temporary = 0. ;
    y_displacement_parameter_temporary = 0. ;
    x_displacement_temporary = 0. ;
    y_displacement_temporary = 0. ;
    inverse_distance_estimator = 0. ;
    x_error = 0. ;
    y_error = 0. ;
    error_norm = 0. ;

    factor_mass_scaling = 1. ;
    x_factor_mass_scaling = 1. ;
    y_factor_mass_scaling = 1. ;
    delta_x_factor_mass_scaling = 1. ;
    delta_y_factor_mass_scaling = 1. ;
    mass_mass_scaling = 0. ;

}

Node::~Node()
{
}

// Accessors

// Modifiers

// Methods
void Node::Sum_up_forces()
{
    x_internal_force = 0. ;
    y_internal_force = 0. ;
    for (int i(0) ; i<(int)x_regions_internal_forces.size() ; i++)
    {
        x_internal_force += x_regions_internal_forces[i] ;
        y_internal_force += y_regions_internal_forces[i] ;
        x_regions_internal_forces[i] = 0. ;
        y_regions_internal_forces[i] = 0. ;
    }
    for (int i(0) ; i<(int)x_master_contact_forces.size() ; i++)
    {
        x_contact_force += x_master_contact_forces[i] ;
        y_contact_force += y_master_contact_forces[i] ;
        x_master_contact_forces[i] = 0. ;
        y_master_contact_forces[i] = 0. ;
    }
    x_force = x_internal_force + x_contact_force + x_self_contact_force + x_body_force + x_dirichlet_force + x_neumann_force + x_damping_force + x_alid_force ;
    y_force = y_internal_force + y_contact_force + y_self_contact_force + y_body_force + y_dirichlet_force + y_neumann_force + y_damping_force + y_alid_force ;
}

void Node::Apply_Newton()
{
    x_acceleration_parameter = x_force * x_inverse_mass / x_factor_mass_scaling ;
    y_acceleration_parameter = y_force * y_inverse_mass / y_factor_mass_scaling ;
}

void Node::Update_current_positions()
{
    x_current = x_initial + x_displacement ;
    y_current = y_initial + y_displacement ;
}

void Node::Apply_Euler(double Deltat)
{
    x_velocity_parameter = x_velocity_parameter + Deltat * x_acceleration_parameter ;
    y_velocity_parameter = y_velocity_parameter + Deltat * y_acceleration_parameter ;
    x_displacement_parameter = x_displacement_parameter + Deltat * x_velocity_parameter ;
    y_displacement_parameter = y_displacement_parameter + Deltat * y_velocity_parameter ;
}

void Node::Apply_Euler_temporary(double Deltat)
{
    x_displacement_parameter_temporary = x_displacement_parameter + Deltat * (x_velocity_parameter + Deltat * x_acceleration_parameter) ;
    y_displacement_parameter_temporary = y_displacement_parameter + Deltat * (y_velocity_parameter + Deltat * y_acceleration_parameter) ;
}

void Node::Compute_error(double& total_error, double& max_error, int& node_for_max_error, double& flagdrivenx, double& flagdriveny, int index)
{
    if (flagdrivenx == 1.)
        x_error = 0. ;
    else
        x_error = x_displacement - x_displacement_temporary ;
    if (flagdriveny == 1.)
        y_error = 0. ;
    else
        y_error = y_displacement - y_displacement_temporary ;
    error_norm = sqrt( x_error * x_error + y_error * y_error ) * inverse_distance_estimator ;
    total_error += error_norm ;
    if (error_norm==0.)
        error_norm = 1.e-16 ;
    if (error_norm > max_error)
    {
        max_error = error_norm ;
        node_for_max_error = index ;
    }
}

void Node::Compute_mass_scaling(double Target_error, double Inv_Target_error, double Control_parameter_mass_scaling, double Max_mass_scaling,
                                double Error_factor_mass_scaling, double Accepted_ratio, double Decrease_factor_mass_scaling, double& max_factor_mass_scaling, int index, int num_node, double& total_mass_mass_scaling)
{
    if (abs(x_error * inverse_distance_estimator) >= Target_error * Error_factor_mass_scaling)
    {
        delta_x_factor_mass_scaling = pow((abs(x_error * inverse_distance_estimator) * Inv_Target_error / Error_factor_mass_scaling), Control_parameter_mass_scaling) ;
        x_factor_mass_scaling *= delta_x_factor_mass_scaling ;
        if (x_factor_mass_scaling > Max_mass_scaling)
            x_factor_mass_scaling = Max_mass_scaling ;
    }
    else
    {
        if (x_factor_mass_scaling != 1) x_factor_mass_scaling *= Decrease_factor_mass_scaling ;
        if (x_factor_mass_scaling < 1) x_factor_mass_scaling = 1;
    }

    if (abs(y_error * inverse_distance_estimator) >= Target_error * Error_factor_mass_scaling)
    {
        delta_y_factor_mass_scaling = pow((abs(y_error * inverse_distance_estimator) * Inv_Target_error / Error_factor_mass_scaling), Control_parameter_mass_scaling) ;
        y_factor_mass_scaling *= delta_y_factor_mass_scaling ;
        if (y_factor_mass_scaling > Max_mass_scaling)
            y_factor_mass_scaling = Max_mass_scaling ;
    }
    else
    {
        if (y_factor_mass_scaling != 1) y_factor_mass_scaling *= Decrease_factor_mass_scaling ;
        if (y_factor_mass_scaling < 1) y_factor_mass_scaling = 1;
    }

    factor_mass_scaling = (x_factor_mass_scaling + y_factor_mass_scaling) * 0.5 ;

    if(factor_mass_scaling > max_factor_mass_scaling)
    {
        max_factor_mass_scaling = factor_mass_scaling ;
    }

    mass_mass_scaling = sqrt(x_mass * x_mass * x_factor_mass_scaling * x_factor_mass_scaling + y_mass * y_mass * y_factor_mass_scaling * y_factor_mass_scaling) ;
    total_mass_mass_scaling += mass_mass_scaling ;
}

#endif
