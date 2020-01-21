#ifndef DEF_PROXIMITY
#define DEF_PROXIMITY

//********************************************//
//** DETECT SEGMENT PROXIMITY ****************//
//********************************************//

int Detect_Segment_Proximity(double Xs, double Ys, double X1, double Y1, double X2, double Y2, double detection_distance)
{
    // OPTIMIZED
    detection_distance = detection_distance * detection_distance ;
    if ((Xs-X1) * (Xs-X1) + (Ys-Y1) * (Ys-Y1) < detection_distance)
        return 1 ;
    if ((Xs-X2) * (Xs-X2) + (Ys-Y2) * (Ys-Y2) < detection_distance)
        return 1 ;
    double dx = X2-X1 ;
    double dy = Y2-Y1 ;
    double a = dx * dx + dy * dy ;
    double Xsi = -1. + 2. * ( dx * (Xs-X1) + dy * (Ys-Y1) ) / a ;
    if ( (Xsi<1.) && (Xsi>-1.) )
    {
        double b = dy * ( Xs - 0.5 * (X1+X2) - Xsi * 0.5 * dx ) - dx * ( Ys - 0.5 * (Y1+Y2) - Xsi * 0.5 * dy ) ;
        if (b * b < a * detection_distance)
            return 1 ;
    }
    return 0 ;
}
/*
{
    if ((Xs-X1) * (Xs-X1) + (Ys-Y1) * (Ys-Y1) < detection_distance * detection_distance) return 1 ;
    if ((Xs-X2) * (Xs-X2) + (Ys-Y2) * (Ys-Y2) < detection_distance * detection_distance) return 1 ;
    double Xsi = -1. + 2. * ( (X2-X1) * (Xs-X1) + (Y2-Y1) * (Ys-Y1) ) / ( (X2-X1) * (X2-X1) + (Y2-Y1) * (Y2-Y1) ) ;
    if ( (Xsi<1.) && (Xsi>-1.) )
    {
        if (abs((Y2-Y1)*(Xs-0.5*(X1+X2)-Xsi*0.5*(X2-X1))-(X2-X1)*(Ys-0.5*(Y1+Y2)-Xsi*0.5*(Y2-Y1)))<sqrt( (X2-X1) * (X2-X1) + (Y2-Y1) * (Y2-Y1) )*detection_distance)
            return 1 ;
    }
    return 0 ;
}
*/



//********************************************//
//** UPDATE PROXIMITY ************************//
//********************************************//

void Update_proximity(
    int Nb_bodies,
    vector<Body>& Bodies,
    double Xmin_period,
    double Xmax_period,
    vector<int>& flags )
{
    //cout << endl ;
    if( flags[7] == 0 )
        cout << "Updating proximities" << endl ;
    double period = Xmax_period - Xmin_period ;
    double invperiod = 1. / period ;
    int flagok ;

    // Determining the spatial extension of each border of each body
    #pragma omp parallel
    {
        #pragma omp for schedule(dynamic)
        for (int i=0 ; i<Nb_bodies ; i++)
        {
            double xmin, xmax, ymin, ymax, x, y, d ;
            int shiftmin, shiftmax, nb_border_nodes ;
            vector<int> border_nodes ;
            Bodies[i].neighbours.clear() ;
            Bodies[i].nb_neighbours = 0 ;
            if (Bodies[i].status == "inactive")
                continue ;
            for (int j(0) ; j<Bodies[i].nb_borders ; j++)
            {
                xmin = 1.e100 ;
                xmax = -1.e100 ;
                ymin = 1.e100 ;
                ymax = -1.e100 ;
                nb_border_nodes=Bodies[i].borders[j].number_border_nodes ;
                border_nodes=Bodies[i].borders[j].border_nodes ;
                for (int k(0) ; k<nb_border_nodes ; k++)
                {
                    x = Bodies[i].nodes[border_nodes[k]].x_current ;
                    y = Bodies[i].nodes[border_nodes[k]].y_current ;
                    //if (Bodies[i].type=="deformable") d = Bodies[i].nodes[border_nodes[k]].domain_size * 0.25 ;
                    //else if (Bodies[i].type=="rigid") d = Bodies[i].nodal_distance * 0.5 ;
                    d = Bodies[i].detection_distance ;
                    //
                    if (x - d < xmin)
                        xmin = x - d ;
                    if (x + d > xmax)
                        xmax = x + d ;
                    if (y - d < ymin)
                        ymin = y - d ;
                    if (y + d > ymax)
                        ymax = y + d ;
                }
                shiftmin = floor ( ( xmin - Xmin_period ) * invperiod ) ;
                shiftmax = floor ( ( xmax - Xmin_period ) * invperiod ) ;
                Bodies[i].borders[j].xmin_prox = xmin ;
                Bodies[i].borders[j].xmax_prox = xmax ;
                Bodies[i].borders[j].ymin_prox = ymin ;
                Bodies[i].borders[j].ymax_prox = ymax ;
                Bodies[i].borders[j].shiftmin = shiftmin ;
                Bodies[i].borders[j].shiftmax = shiftmax ;
            }
        }

        // Computing global proximities
        #pragma omp for schedule(dynamic)
        for (int body1=0 ; body1<Nb_bodies ; body1++)
        {
            if (Bodies[body1].status == "inactive")
                continue ;
            double xmin1, xmax1, ymin1, ymax1, xmin2, xmax2, ymin2, ymax2 ;
            int shiftmin1, shiftmax1, shiftmin2, shiftmax2 ;
            for (int border1(0) ; border1<Bodies[body1].nb_borders ; border1++)
            {
                xmin1 = Bodies[body1].borders[border1].xmin_prox ;
                xmax1 = Bodies[body1].borders[border1].xmax_prox ;
                ymin1 = Bodies[body1].borders[border1].ymin_prox ;
                ymax1 = Bodies[body1].borders[border1].ymax_prox ;
                shiftmin1 = Bodies[body1].borders[border1].shiftmin ;
                shiftmax1 = Bodies[body1].borders[border1].shiftmax ;
                for (int body2(0) ; body2<Nb_bodies ; body2++)
                {
                    if (Bodies[body2].status == "inactive")
                        continue ;
                    if (body1==body2)
                        continue ;
                    for (int border2(0) ; border2<Bodies[body2].nb_borders ; border2++)
                    {
                        xmin2 = Bodies[body2].borders[border2].xmin_prox ;
                        xmax2 = Bodies[body2].borders[border2].xmax_prox ;
                        ymin2 = Bodies[body2].borders[border2].ymin_prox ;
                        ymax2 = Bodies[body2].borders[border2].ymax_prox ;
                        shiftmin2 = Bodies[body2].borders[border2].shiftmin ;
                        shiftmax2 = Bodies[body2].borders[border2].shiftmax ;
                        if ( (ymax2<ymin1) || (ymax1<ymin2) )
                        {
                            continue ;
                        }
                        for (int shift(shiftmin1 - shiftmax2) ; shift <= shiftmax1 - shiftmin2 ; shift++)
                        {
                            if ( (xmax2+shift*period>xmin1) && (xmin2+shift*period<xmax1) )
                            {
                                flagok = 1 ;
                                for (int i(0) ; i<Bodies[body1].nb_neighbours ; i++)
                                {
                                    if (Bodies[body1].neighbours[i][0]==body2 && Bodies[body1].neighbours[i][1]==border2 && Bodies[body1].neighbours[i][2]==shift)
                                    {
                                        flagok = 0 ;
                                        break ;
                                    }
                                }
                                if (flagok ==1)
                                {
                                    Bodies[body1].neighbours.push_back({ body2, border2, shift }) ;
                                    Bodies[body1].nb_neighbours++ ;
                                }

                            }
                        }
                    }
                }
            }
        }

        // Initilizing master contact forces
        /*
        #pragma omp for schedule(dynamic)
        for (int body=0 ; body<Nb_bodies ; body++)
        {
            vector<double> v(Bodies[body].nb_neighbours) ;
            for (int i(0) ; i<Bodies[body].nb_neighbours ; i++)
                v[i] = 0. ;
            if (Bodies[body].type=="deformable")
            {
                for (int i(0) ; i<Bodies[body].nb_nodes ; i++)
                {
                    Bodies[body].nodes[i].x_master_contact_forces = v ;
                    Bodies[body].nodes[i].y_master_contact_forces = v ;
                }
            }
            else if (Bodies[body].type=="rigid")
            {
                Bodies[body].x_master_contact_forces = v ;
                Bodies[body].y_master_contact_forces = v ;
                Bodies[body].r_master_contact_forces = v ;
            }
        }
        */

        // Computing local proximities
        #pragma omp for schedule(dynamic)
        for (int bodyS=0 ; bodyS<Nb_bodies ; bodyS++)
        {
            int flag_detect ;
            int nb_border_nodesS, nb_border_nodesM, bodyM, borderM, shiftM ;
            int flag_exists ;
            int NodeS, NodeM0, ShiftM0, NodeM1, ShiftM1, NodeM2, ShiftM2, NodeM3, ShiftM3 ;
            double Xs, Ys, X0, Y0, X1, Y1, X2, Y2, X3, Y3, N0, N1, N2, N3 ;
            double Gapn, Xsi, Xnorm, Ynorm, Xtan, Ytan, Xclosest, Yclosest, length, effective_mass ;
            vector<int> border_nodesS, NodesM0, ShiftsM0, NodesM1, ShiftsM1, NodesM2, ShiftsM2, NodesM3, ShiftsM3, Proximity_previous ;
            vector<double> Xsis ;
            vector<double> internal ;
            string interpolantM ;
            int nb_contact_elements = 0 ;
            Contact_element new_contact ;
            vector<Contact_element> contact_elements ;
            for (int borderS(0) ; borderS<Bodies[bodyS].nb_borders ; borderS++)
            {
                nb_border_nodesS=Bodies[bodyS].borders[borderS].number_border_nodes ;
                for (int neighbour(0) ; neighbour<Bodies[bodyS].nb_neighbours ; neighbour++)
                {
                    bodyM = Bodies[bodyS].neighbours[neighbour][0] ;
                    borderM = Bodies[bodyS].neighbours[neighbour][1] ;
                    interpolantM = Bodies[bodyM].borders[borderM].interpolant ;
                    shiftM = Bodies[bodyS].neighbours[neighbour][2] ;
                    NodesM0 = Bodies[bodyM].borders[borderM].node0 ;
                    ShiftsM0 = Bodies[bodyM].borders[borderM].shift0 ;
                    NodesM1 = Bodies[bodyM].borders[borderM].node1 ;
                    ShiftsM1 = Bodies[bodyM].borders[borderM].shift1 ;
                    NodesM2 = Bodies[bodyM].borders[borderM].node2 ;
                    ShiftsM2 = Bodies[bodyM].borders[borderM].shift2 ;
                    NodesM3 = Bodies[bodyM].borders[borderM].node3 ;
                    ShiftsM3 = Bodies[bodyM].borders[borderM].shift3 ;
                    nb_border_nodesM = Bodies[bodyM].borders[borderM].number_border_nodes ;
                    border_nodesS = Bodies[bodyS].borders[borderS].border_nodes ;
                    if      ( Bodies[bodyM].mass == 0. )
                        effective_mass = 0. ;
                    else if ( Bodies[bodyS].mass == 0. )
                        effective_mass = 0. ;
                    else
                        effective_mass = 1. / (Bodies[bodyM].inverse_mass + Bodies[bodyS].inverse_mass) ;
                    for (int border_nodeS(0) ; border_nodeS<nb_border_nodesS ; border_nodeS++)
                    {
                        NodeS = border_nodesS[border_nodeS] ;
                        length = Bodies[bodyS].borders[borderS].length[border_nodeS] ;
                        Xs = Bodies[bodyS].nodes[NodeS].x_current ;
                        Ys = Bodies[bodyS].nodes[NodeS].y_current ;
                        for (int segmentM(0) ; segmentM<nb_border_nodesM-1 ; segmentM++)
                        {
                            NodeM1 = NodesM1[segmentM] ;
                            ShiftM1 = ShiftsM1[segmentM] ;
                            X1 = Bodies[bodyM].nodes[NodeM1].x_current + period * ( shiftM + ShiftM1 ) ;
                            Y1 = Bodies[bodyM].nodes[NodeM1].y_current ;
                            NodeM2 = NodesM2[segmentM] ;
                            ShiftM2 = ShiftsM2[segmentM] ;
                            X2 = Bodies[bodyM].nodes[NodeM2].x_current + period * ( shiftM + ShiftM2 ) ;
                            Y2 = Bodies[bodyM].nodes[NodeM2].y_current ;
                            if (Detect_Segment_Proximity(Xs,Ys,X1,Y1,X2,Y2,Bodies[bodyS].detection_distance)==0)
                                continue ;
                            NodeM0 = NodesM0[segmentM] ;
                            ShiftM0 = ShiftsM0[segmentM] ;
                            X0 = Bodies[bodyM].nodes[NodeM0].x_current + period * ( shiftM + ShiftM0 ) ;
                            Y0 = Bodies[bodyM].nodes[NodeM0].y_current ;
                            NodeM3 = NodesM3[segmentM] ;
                            ShiftM3 = ShiftsM3[segmentM] ;
                            X3 = Bodies[bodyM].nodes[NodeM3].x_current + period * ( shiftM + ShiftM3 ) ;
                            Y3 = Bodies[bodyM].nodes[NodeM3].y_current ;
                            Xsi = 0. ;
                            Closest_2_segment(interpolantM, Xs, Ys, X0, Y0, X1, Y1, X2, Y2, X3, Y3,
                                              1.e-16, Gapn, Xsi,
                                              Xnorm, Ynorm, Xtan, Ytan, Xclosest, Yclosest,
                                              N0, N1, N2, N3, flag_detect) ;
                            if ( (flag_detect==-2 || flag_detect==-1 || flag_detect==0) && (Gapn>-Bodies[bodyS].contact_distance) )
                            {
                                flag_exists = 0 ;
                                for (int k(0) ; k<Bodies[bodyS].nb_contact_elements ; k++)
                                {
                                    if ((Bodies[bodyS].contact_elements[k].borderS == borderS) &&
                                            (Bodies[bodyS].contact_elements[k].border_nodeS == border_nodeS) &&
                                            (Bodies[bodyS].contact_elements[k].bodyM == bodyM) &&
                                            (Bodies[bodyS].contact_elements[k].borderM == borderM) &&
                                            (abs(Bodies[bodyS].contact_elements[k].segmentM - segmentM) <= 1.))
                                    {
                                        new_contact = Bodies[bodyS].contact_elements[k] ;
                                        flag_exists = 1 ;
                                        break ;
                                    }
                                }
                                if (flag_exists == 0)
                                {
                                    new_contact.borderS = borderS ;
                                    new_contact.border_nodeS = border_nodeS ;
                                    new_contact.nodeS = Bodies[bodyS].borders[borderS].border_nodes[border_nodeS] ;
                                    new_contact.bodyM = bodyM ;
                                    new_contact.borderM = borderM ;
                                    new_contact.shiftM = shiftM ;
                                    new_contact.segmentM = segmentM ;
                                    new_contact.nodeM0 = Bodies[bodyM].borders[borderM].node0[segmentM] ;
                                    new_contact.shiftM0  = Bodies[bodyM].borders[borderM].shift0[segmentM] ;
                                    new_contact.shapeM0 = N0 ;
                                    new_contact.nodeM1 = Bodies[bodyM].borders[borderM].node1[segmentM] ; ;
                                    new_contact.shiftM1  = Bodies[bodyM].borders[borderM].shift1[segmentM] ;
                                    new_contact.shapeM1 = N1 ;
                                    new_contact.nodeM2 = Bodies[bodyM].borders[borderM].node2[segmentM] ; ;
                                    new_contact.shiftM2  = Bodies[bodyM].borders[borderM].shift2[segmentM] ;
                                    new_contact.shapeM2 = N2 ;
                                    new_contact.nodeM3 = Bodies[bodyM].borders[borderM].node3[segmentM] ; ;
                                    new_contact.shiftM3  = Bodies[bodyM].borders[borderM].shift3[segmentM] ;
                                    new_contact.shapeM3 = N3 ;
                                    new_contact.gapn = Gapn ;
                                    new_contact.gapt = 0. ;
                                    new_contact.xsi = Xsi ;
                                    new_contact.xnorm = Xnorm ;
                                    new_contact.ynorm = Ynorm ;
                                    new_contact.xtan = Xtan ;
                                    new_contact.ytan = Ytan ;
                                    new_contact.effective_mass = effective_mass ;
                                    new_contact.length = length ;
                                    new_contact.fx = 0. ;
                                    new_contact.fy = 0. ;
                                    new_contact.nb_internal = 0 ;
                                    new_contact.internal = {0.} ;
                                }
                                contact_elements.push_back(new_contact) ;
                                nb_contact_elements++ ;
                            }
                        }
                    }
                }
            }
            Bodies[bodyS].nb_contact_elements = nb_contact_elements ;
            Bodies[bodyS].contact_elements = contact_elements ;
        }

        // Computing self proximities
        if ( flags[8] == 0 )
        {
            #pragma omp for schedule(dynamic)
            for (int bodyS=0 ; bodyS<Nb_bodies ; bodyS++)
            {
                if (Bodies[bodyS].type=="rigid" || Bodies[bodyS].periodicity=="Periodic")
                    continue ;
                int flag_detect ;
                int nb_border_nodesS, nb_border_nodesM; //nb_internal ; //borderM,
                int flag_exists ;
                int NodeS, NodeM0, NodeM1, NodeM2, NodeM3 ;
                double Xs, Ys, X0, Y0, X1, Y1, X2, Y2, X3, Y3, N0, N1, N2, N3 ;
                double Gapn, Xsi, Xnorm, Ynorm, Xtan, Ytan, Xclosest, Yclosest, length, effective_mass ; //Gapt,
                vector<int> border_nodesS, NodesM0, NodesM1, NodesM2, NodesM3, Proximity_previous ;
                vector<double> Xsis ;
                vector<double> internal ;
                string interpolantM ;
                int nb_contact_elements = Bodies[bodyS].nb_contact_elements ;
                Contact_element new_contact ;
                vector<Contact_element> contact_elements = Bodies[bodyS].contact_elements ;
                effective_mass = Bodies[bodyS].mass ;
                for (int borderS(0) ; borderS<Bodies[bodyS].nb_borders ; borderS++)
                {
                    border_nodesS = Bodies[bodyS].borders[borderS].border_nodes ;
                    nb_border_nodesS=Bodies[bodyS].borders[borderS].number_border_nodes ;
                    for (int border_nodeS(0) ; border_nodeS<nb_border_nodesS ; border_nodeS++)
                    {
                        NodeS = border_nodesS[border_nodeS] ;
                        length = Bodies[bodyS].borders[borderS].length[border_nodeS] ;
                        Xs = Bodies[bodyS].nodes[NodeS].x_current ;
                        Ys = Bodies[bodyS].nodes[NodeS].y_current ;
                        for (int borderM(0) ; borderM<Bodies[bodyS].nb_borders ; borderM++)
                        {
                            interpolantM = Bodies[bodyS].borders[borderM].interpolant ;
                            NodesM0 = Bodies[bodyS].borders[borderM].node0 ;
                            NodesM1 = Bodies[bodyS].borders[borderM].node1 ;
                            NodesM2 = Bodies[bodyS].borders[borderM].node2 ;
                            NodesM3 = Bodies[bodyS].borders[borderM].node3 ;
                            nb_border_nodesM = Bodies[bodyS].borders[borderM].number_border_nodes ;
                            for (int segmentM(0) ; segmentM<nb_border_nodesM-1 ; segmentM++)
                            {
                                NodeM1 = NodesM1[segmentM] ;
                                if (NodeM1==NodeS)
                                    continue ;
                                X1 = Bodies[bodyS].nodes[NodeM1].x_current ;
                                Y1 = Bodies[bodyS].nodes[NodeM1].y_current ;
                                NodeM2 = NodesM2[segmentM] ;
                                if (NodeM2==NodeS)
                                    continue ;
                                X2 = Bodies[bodyS].nodes[NodeM2].x_current ;
                                Y2 = Bodies[bodyS].nodes[NodeM2].y_current ;
                                if (Detect_Segment_Proximity(Xs,Ys,X1,Y1,X2,Y2,Bodies[bodyS].detection_distance)==0)
                                    continue ;
                                NodeM0 = NodesM0[segmentM] ;
                                X0 = Bodies[bodyS].nodes[NodeM0].x_current ;
                                Y0 = Bodies[bodyS].nodes[NodeM0].y_current ;
                                NodeM3 = NodesM3[segmentM] ;
                                X3 = Bodies[bodyS].nodes[NodeM3].x_current ;
                                Y3 = Bodies[bodyS].nodes[NodeM3].y_current ;
                                Xsi = 0. ;
                                if (NodeM0==NodeS || NodeM3==NodeS)
                                {
                                    Closest_2_segment_self(Xs, Ys, X1, Y1, X2, Y2,
                                                           1.e-16, Gapn, Xsi,
                                                           Xnorm, Ynorm, Xtan, Ytan, Xclosest, Yclosest,
                                                           N0, N1, N2, N3, flag_detect) ;
                                }
                                else
                                {
                                    Closest_2_segment(interpolantM, Xs, Ys, X0, Y0, X1, Y1, X2, Y2, X3, Y3,
                                                      1.e-16, Gapn, Xsi,
                                                      Xnorm, Ynorm, Xtan, Ytan, Xclosest, Yclosest,
                                                      N0, N1, N2, N3, flag_detect) ;
                                }
                                //if ( (-1.<=Xsi) && (Xsi<1.) && (Gapn>-Bodies[bodyS].contact_distance) )
                                if ( flag_detect==0 && (Gapn>-Bodies[bodyS].contact_distance) )
                                {
                                    flag_exists = 0 ;
                                    for (int k(0) ; k<Bodies[bodyS].nb_contact_elements ; k++)
                                    {
                                        //if ((Bodies[bodyS].contact_elements[k].borderS == borderS) &&
                                        //    (Bodies[bodyS].contact_elements[k].border_nodeS == border_nodeS) &&
                                        //    (Bodies[bodyS].contact_elements[k].bodyM == bodyS))
                                        if ((Bodies[bodyS].contact_elements[k].borderS == borderS) &&
                                                (Bodies[bodyS].contact_elements[k].border_nodeS == border_nodeS) &&
                                                (Bodies[bodyS].contact_elements[k].bodyM == bodyS) &&
                                                (Bodies[bodyS].contact_elements[k].borderM == borderM) &&
                                                (abs(Bodies[bodyS].contact_elements[k].segmentM - segmentM) <= 1.))
                                        {
                                            new_contact = Bodies[bodyS].contact_elements[k] ;
                                            flag_exists = 1 ;
                                            break ;
                                        }
                                    }
                                    if (flag_exists == 0)
                                    {
                                        new_contact.borderS = borderS ;
                                        new_contact.border_nodeS = border_nodeS ;
                                        new_contact.nodeS = Bodies[bodyS].borders[borderS].border_nodes[border_nodeS] ;
                                        new_contact.bodyM = bodyS ;
                                        new_contact.borderM = borderM ;
                                        new_contact.shiftM = 0 ;
                                        new_contact.segmentM = segmentM ;
                                        new_contact.nodeM0 = Bodies[bodyS].borders[borderM].node0[segmentM] ;
                                        new_contact.shiftM0  = 0 ;
                                        new_contact.shapeM0 = N0 ;
                                        new_contact.nodeM1 = Bodies[bodyS].borders[borderM].node1[segmentM] ; ;
                                        new_contact.shiftM1  = 0 ;
                                        new_contact.shapeM1 = N1 ;
                                        new_contact.nodeM2 = Bodies[bodyS].borders[borderM].node2[segmentM] ; ;
                                        new_contact.shiftM2  = 0 ;
                                        new_contact.shapeM2 = N2 ;
                                        new_contact.nodeM3 = Bodies[bodyS].borders[borderM].node3[segmentM] ; ;
                                        new_contact.shiftM3  = 0 ;
                                        new_contact.shapeM3 = N3 ;
                                        new_contact.gapn = Gapn ;
                                        new_contact.gapt = 0. ;
                                        new_contact.xsi = Xsi ;
                                        new_contact.xnorm = Xnorm ;
                                        new_contact.ynorm = Ynorm ;
                                        new_contact.xtan = Xtan ;
                                        new_contact.ytan = Ytan ;
                                        new_contact.effective_mass = effective_mass ;
                                        new_contact.length = length ;
                                        new_contact.fx = 0. ;
                                        new_contact.fy = 0. ;
                                        new_contact.nb_internal = 0 ;
                                        new_contact.internal = {0.} ;
                                    }
                                    contact_elements.push_back(new_contact) ;
                                    nb_contact_elements++ ;
                                }
                            }
                        }
                    }
                }
                Bodies[bodyS].nb_contact_elements = nb_contact_elements ;
                Bodies[bodyS].contact_elements = contact_elements ;
            }
        }
    }
}



//********************************************//
//** INITIALIZE CZM **************************//
//********************************************//

void Initialize_CZM(int Nb_bodies, vector<Body>& Bodies, int Nb_contact_laws, vector<Contact_law>& Contact_laws, vector<int>& flags, double Xmin_period, double Xmax_period)
{
    cout << "Initializing cohesive bonds" << endl ;
    for (int i=0 ; i<Nb_bodies ; i++)
        Bodies[i].Update_contacts(Bodies, Xmin_period, Xmax_period) ;
    flags[3] = 0 ;
    int bodyM ;
    string contact_law_type, material_nameM, material_nameS, material1, material2 ;
    for (int bodyS(0) ; bodyS<Nb_bodies ; bodyS++)
    {
        for (int icontact(0) ; icontact<Bodies[bodyS].nb_contact_elements ; icontact++)
        {
            bodyM = Bodies[bodyS].contact_elements[icontact].bodyM ;
            material_nameM = Bodies[bodyM].material_name ;
            material_nameS = Bodies[bodyS].material_name ;
            vector<double> parameters ;
            for (int i(0) ; i<Nb_contact_laws ; i++)
            {
                material1 = Contact_laws[i].material1 ;
                material2 = Contact_laws[i].material2 ;
                if (((material1==material_nameS) && (material2==material_nameM)) ||
                        ((material2==material_nameS) && (material1==material_nameM)))
                {
                    contact_law_type = Contact_laws[i].type ;
                    parameters = Contact_laws[i].parameters ;
                    break ;
                }
            }
            if (contact_law_type=="CZMlinear")
            {
                double gapinit = parameters[5] ;
                Bodies[bodyS].contact_elements[icontact].nb_internal = 1 ;
                if ( Bodies[bodyS].contact_elements[icontact].gapn < gapinit )
                    Bodies[bodyS].contact_elements[icontact].internal = {0.} ;
                else
                    Bodies[bodyS].contact_elements[icontact].internal = {1.} ;
            }
            else if (contact_law_type=="CZMfatigue")
            {
                double gapinit = parameters[8] ;
                Bodies[bodyS].contact_elements[icontact].nb_internal = 3 ;
                if ( Bodies[bodyS].contact_elements[icontact].gapn < gapinit )
                    Bodies[bodyS].contact_elements[icontact].internal = {0., 0., 0.} ;
                else
                    Bodies[bodyS].contact_elements[icontact].internal = {1., 0., 0.} ;
            }
            else if (contact_law_type=="BondedMohrCoulomb")
            {
                double gapinit = parameters[10] ;
                Bodies[bodyS].contact_elements[icontact].nb_internal = 3 ;
                if ( Bodies[bodyS].contact_elements[icontact].gapn < gapinit )
                    Bodies[bodyS].contact_elements[icontact].internal = {0., 0., 0.} ;
                else
                    Bodies[bodyS].contact_elements[icontact].internal = {1., 0., 0.} ;
            }

            // NB : other contact laws to initialize ?
        }
    }

}




//********************************************//
//** UPDATE STATUS ***************************//
//********************************************//

void Update_status(
    int Nb_bodies,
    vector<Body>& Bodies,
    int Nb_deactivated,
    vector<vector<double>>& Deactivated)
{
    double x ;
    double y ;
    double r ;
    double m ;
    for (int i(0) ; i < Nb_deactivated ; i++)
    {
        if (Deactivated[i][0] == 0)
        {
            //cout << Desactivated[i][0] << Desactivated[i][1] << Desactivated[i][2] << Desactivated[i][3] << endl ;
            if (Deactivated[i][3] >=0 )
            {
                x = Bodies[Deactivated[i][2]].nodes[Deactivated[i][3]].x_current ;
                y = Bodies[Deactivated[i][2]].nodes[Deactivated[i][3]].y_current ;
            }
            else
            {
                if (Bodies[Deactivated[i][2]].type == "deformable")
                {
                    x = 0. ;
                    y = 0. ;
                    m = 0. ;
                    for (int j(0) ; j<Bodies[Deactivated[i][2]].nb_nodes ; j++)
                    {
                        x += Bodies[Deactivated[i][2]].nodes[j].x_current * Bodies[Deactivated[i][2]].nodes[j].x_mass ;
                        y += Bodies[Deactivated[i][2]].nodes[j].y_current * Bodies[Deactivated[i][2]].nodes[j].y_mass ;
                        m += Bodies[Deactivated[i][2]].nodes[j].x_mass ;
                    }
                    x /= m ;
                    y /= m ;
                }
                else if (Bodies[Deactivated[i][2]].type == "rigid")
                {
                    x = Bodies[Deactivated[i][2]].x_current ;
                    y = Bodies[Deactivated[i][2]].y_current ;
                    r = Bodies[Deactivated[i][2]].r_current ;
                }
            }
            if (Deactivated[i][1] == 0)
            {
                if ((Deactivated[i][4] == 0) & (x <= Deactivated[i][5]))
                    Bodies[Deactivated[i][2]].status = "inactive" ;
                else if ((Deactivated[i][4] == 1) & (x == Deactivated[i][5]))
                    Bodies[Deactivated[i][2]].status = "inactive" ;
                else if ((Deactivated[i][4] == 2) & (x >= Deactivated[i][5]))
                    Bodies[Deactivated[i][2]].status = "inactive" ;
            }
            else if (Deactivated[i][1] == 1)
            {
                if ((Deactivated[i][4] == 0) & (y <= Deactivated[i][5]))
                    Bodies[Deactivated[i][2]].status = "inactive" ;
                else if ((Deactivated[i][4] == 1) & (y == Deactivated[i][5]))
                    Bodies[Deactivated[i][2]].status = "inactive" ;
                else if ((Deactivated[i][4] == 2) & (y >= Deactivated[i][5]))
                    Bodies[Deactivated[i][2]].status = "inactive" ;
            }
            else if (Deactivated[i][1] == 2)
            {
                if ((Deactivated[i][4] == 0) & (r <= Deactivated[i][5]))
                    Bodies[Deactivated[i][2]].status = "inactive" ;
                else if ((Deactivated[i][4] == 1) & (r == Deactivated[i][5]))
                    Bodies[Deactivated[i][2]].status = "inactive" ;
                else if ((Deactivated[i][4] == 2) & (r >= Deactivated[i][5]))
                    Bodies[Deactivated[i][2]].status = "inactive" ;
            }
        }
        else if (Deactivated[i][0] == 1)
        {
            //cout << Desactivated[i][0] << Desactivated[i][1] << Desactivated[i][2] << Desactivated[i][3] << endl ;
            if (Deactivated[i][3] >=0 )
            {
                x = Bodies[Deactivated[i][2]].nodes[Deactivated[i][3]].x_displacement ;
                y = Bodies[Deactivated[i][2]].nodes[Deactivated[i][3]].y_displacement ;
            }
            else
            {
                if (Bodies[Deactivated[i][2]].type == "deformable")
                {
                    x = 0. ;
                    y = 0. ;
                    m = 0. ;
                    for (int j(0) ; j<Bodies[Deactivated[i][2]].nb_nodes ; j++)
                    {
                        x += Bodies[Deactivated[i][2]].nodes[j].x_displacement * Bodies[Deactivated[i][2]].nodes[j].x_mass ;
                        y += Bodies[Deactivated[i][2]].nodes[j].y_displacement * Bodies[Deactivated[i][2]].nodes[j].y_mass ;
                        m += Bodies[Deactivated[i][2]].nodes[j].x_mass ;
                    }
                    x /= m ;
                    y /= m ;
                }
                else if (Bodies[Deactivated[i][2]].type == "rigid")
                {
                    x = Bodies[Deactivated[i][2]].x_displacement ;
                    y = Bodies[Deactivated[i][2]].y_displacement ;
                    r = Bodies[Deactivated[i][2]].r_displacement ;
                }
            }
            if (Deactivated[i][1] == 0)
            {
                if ((Deactivated[i][4] == 0) & (x <= Deactivated[i][5]))
                    Bodies[Deactivated[i][2]].status = "inactive" ;
                else if ((Deactivated[i][4] == 1) & (x == Deactivated[i][5]))
                    Bodies[Deactivated[i][2]].status = "inactive" ;
                else if ((Deactivated[i][4] == 2) & (x >= Deactivated[i][5]))
                    Bodies[Deactivated[i][2]].status = "inactive" ;
            }
            else if (Deactivated[i][1] == 1)
            {
                if ((Deactivated[i][4] == 0) & (y <= Deactivated[i][5]))
                    Bodies[Deactivated[i][2]].status = "inactive" ;
                else if ((Deactivated[i][4] == 1) & (y == Deactivated[i][5]))
                    Bodies[Deactivated[i][2]].status = "inactive" ;
                else if ((Deactivated[i][4] == 2) & (y >= Deactivated[i][5]))
                    Bodies[Deactivated[i][2]].status = "inactive" ;
            }
            else if (Deactivated[i][1] == 2)
            {
                if ((Deactivated[i][4] == 0) & (pow(x*x+y+y,0.5) <= Deactivated[i][5]))
                    Bodies[Deactivated[i][2]].status = "inactive" ;
                else if ((Deactivated[i][4] == 1) & (pow(x*x+y+y,0.5) == Deactivated[i][5]))
                    Bodies[Deactivated[i][2]].status = "inactive" ;
                else if ((Deactivated[i][4] == 2) & (pow(x*x+y+y,0.5) >= Deactivated[i][5]))
                    Bodies[Deactivated[i][2]].status = "inactive" ;
            }
            else if (Deactivated[i][1] == 3)
            {
                if ((Deactivated[i][4] == 0) & (r <= Deactivated[i][5]))
                    Bodies[Deactivated[i][2]].status = "inactive" ;
                else if ((Deactivated[i][4] == 1) & (r == Deactivated[i][5]))
                    Bodies[Deactivated[i][2]].status = "inactive" ;
                else if ((Deactivated[i][4] == 2) & (r >= Deactivated[i][5]))
                    Bodies[Deactivated[i][2]].status = "inactive" ;
            }
        }
        else if (Deactivated[i][0] == 2)
        {
            //cout << Desactivated[i][0] << Desactivated[i][1] << Desactivated[i][2] << Desactivated[i][3] << endl ;
            if (Deactivated[i][3] >=0 )
            {
                x = Bodies[Deactivated[i][2]].nodes[Deactivated[i][3]].x_velocity ;
                y = Bodies[Deactivated[i][2]].nodes[Deactivated[i][3]].y_velocity ;
            }
            else
            {
                if (Bodies[Deactivated[i][2]].type == "deformable")
                {
                    x = 0. ;
                    y = 0. ;
                    m = 0. ;
                    for (int j(0) ; j<Bodies[Deactivated[i][2]].nb_nodes ; j++)
                    {
                        x += Bodies[Deactivated[i][2]].nodes[j].x_velocity * Bodies[Deactivated[i][2]].nodes[j].x_mass ;
                        y += Bodies[Deactivated[i][2]].nodes[j].y_velocity * Bodies[Deactivated[i][2]].nodes[j].y_mass ;
                        m += Bodies[Deactivated[i][2]].nodes[j].x_mass ;
                    }
                    x /= m ;
                    y /= m ;
                }
                else if (Bodies[Deactivated[i][2]].type == "rigid")
                {
                    x = Bodies[Deactivated[i][2]].x_velocity ;
                    y = Bodies[Deactivated[i][2]].y_velocity ;
                    r = Bodies[Deactivated[i][2]].r_velocity ;
                }
            }
            if (Deactivated[i][1] == 0)
            {
                if ((Deactivated[i][4] == 0) & (x <= Deactivated[i][5]))
                    Bodies[Deactivated[i][2]].status = "inactive" ;
                else if ((Deactivated[i][4] == 1) & (x == Deactivated[i][5]))
                    Bodies[Deactivated[i][2]].status = "inactive" ;
                else if ((Deactivated[i][4] == 2) & (x >= Deactivated[i][5]))
                    Bodies[Deactivated[i][2]].status = "inactive" ;
            }
            else if (Deactivated[i][1] == 1)
            {
                if ((Deactivated[i][4] == 0) & (y <= Deactivated[i][5]))
                    Bodies[Deactivated[i][2]].status = "inactive" ;
                else if ((Deactivated[i][4] == 1) & (y == Deactivated[i][5]))
                    Bodies[Deactivated[i][2]].status = "inactive" ;
                else if ((Deactivated[i][4] == 2) & (y >= Deactivated[i][5]))
                    Bodies[Deactivated[i][2]].status = "inactive" ;
            }
            else if (Deactivated[i][1] == 2)
            {
                if ((Deactivated[i][4] == 0) & (pow(x*x+y+y,0.5) <= Deactivated[i][5]))
                    Bodies[Deactivated[i][2]].status = "inactive" ;
                else if ((Deactivated[i][4] == 1) & (pow(x*x+y+y,0.5) == Deactivated[i][5]))
                    Bodies[Deactivated[i][2]].status = "inactive" ;
                else if ((Deactivated[i][4] == 2) & (pow(x*x+y+y,0.5) >= Deactivated[i][5]))
                    Bodies[Deactivated[i][2]].status = "inactive" ;
            }
            else if (Deactivated[i][1] == 3)
            {
                if ((Deactivated[i][4] == 0) & (r <= Deactivated[i][5]))
                    Bodies[Deactivated[i][2]].status = "inactive" ;
                else if ((Deactivated[i][4] == 1) & (r == Deactivated[i][5]))
                    Bodies[Deactivated[i][2]].status = "inactive" ;
                else if ((Deactivated[i][4] == 2) & (r >= Deactivated[i][5]))
                    Bodies[Deactivated[i][2]].status = "inactive" ;
            }
        }
        else if (Deactivated[i][0] == 3)
        {
            //cout << Desactivated[i][0] << Desactivated[i][1] << Desactivated[i][2] << Desactivated[i][3] << endl ;
            if (Deactivated[i][3] >=0 )
            {
                x = Bodies[Deactivated[i][2]].nodes[Deactivated[i][3]].x_acceleration ;
                y = Bodies[Deactivated[i][2]].nodes[Deactivated[i][3]].y_acceleration ;
            }
            else
            {
                if (Bodies[Deactivated[i][2]].type == "deformable")
                {
                    x = 0. ;
                    y = 0. ;
                    m = 0. ;
                    for (int j(0) ; j<Bodies[Deactivated[i][2]].nb_nodes ; j++)
                    {
                        x += Bodies[Deactivated[i][2]].nodes[j].x_acceleration * Bodies[Deactivated[i][2]].nodes[j].x_mass ;
                        y += Bodies[Deactivated[i][2]].nodes[j].y_acceleration * Bodies[Deactivated[i][2]].nodes[j].y_mass ;
                        m += Bodies[Deactivated[i][2]].nodes[j].x_mass ;
                    }
                    x /= m ;
                    y /= m ;
                }
                else if (Bodies[Deactivated[i][2]].type == "rigid")
                {
                    x = Bodies[Deactivated[i][2]].x_acceleration ;
                    y = Bodies[Deactivated[i][2]].y_acceleration ;
                    r = Bodies[Deactivated[i][2]].r_acceleration ;
                }
            }
            if (Deactivated[i][1] == 0)
            {
                if ((Deactivated[i][4] == 0) & (x <= Deactivated[i][5]))
                    Bodies[Deactivated[i][2]].status = "inactive" ;
                else if ((Deactivated[i][4] == 1) & (x == Deactivated[i][5]))
                    Bodies[Deactivated[i][2]].status = "inactive" ;
                else if ((Deactivated[i][4] == 2) & (x >= Deactivated[i][5]))
                    Bodies[Deactivated[i][2]].status = "inactive" ;
            }
            else if (Deactivated[i][1] == 1)
            {
                if ((Deactivated[i][4] == 0) & (y <= Deactivated[i][5]))
                    Bodies[Deactivated[i][2]].status = "inactive" ;
                else if ((Deactivated[i][4] == 1) & (y == Deactivated[i][5]))
                    Bodies[Deactivated[i][2]].status = "inactive" ;
                else if ((Deactivated[i][4] == 2) & (y >= Deactivated[i][5]))
                    Bodies[Deactivated[i][2]].status = "inactive" ;
            }
            else if (Deactivated[i][1] == 2)
            {
                if ((Deactivated[i][4] == 0) & (pow(x*x+y+y,0.5) <= Deactivated[i][5]))
                    Bodies[Deactivated[i][2]].status = "inactive" ;
                else if ((Deactivated[i][4] == 1) & (pow(x*x+y+y,0.5) == Deactivated[i][5]))
                    Bodies[Deactivated[i][2]].status = "inactive" ;
                else if ((Deactivated[i][4] == 2) & (pow(x*x+y+y,0.5) >= Deactivated[i][5]))
                    Bodies[Deactivated[i][2]].status = "inactive" ;
            }
            else if (Deactivated[i][1] == 3)
            {
                if ((Deactivated[i][4] == 0) & (r <= Deactivated[i][5]))
                    Bodies[Deactivated[i][2]].status = "inactive" ;
                else if ((Deactivated[i][4] == 1) & (r == Deactivated[i][5]))
                    Bodies[Deactivated[i][2]].status = "inactive" ;
                else if ((Deactivated[i][4] == 2) & (r >= Deactivated[i][5]))
                    Bodies[Deactivated[i][2]].status = "inactive" ;
            }
        }
        else if (Deactivated[i][0] == 4)
        {
            //cout << Desactivated[i][0] << Desactivated[i][1] << Desactivated[i][2] << Desactivated[i][3] << Desactivated[i][4] << endl ;
            if (Deactivated[i][4] >=0 )
            {
                if (Deactivated[i][2] == 0)
                {
                    x = Bodies[Deactivated[i][3]].nodes[Deactivated[i][4]].x_force ;
                    y = Bodies[Deactivated[i][3]].nodes[Deactivated[i][4]].y_force ;
                }
                else if (Deactivated[i][2] == 1)
                {
                    x = Bodies[Deactivated[i][3]].nodes[Deactivated[i][4]].x_internal_force ;
                    y = Bodies[Deactivated[i][3]].nodes[Deactivated[i][4]].y_internal_force ;
                }
                else if (Deactivated[i][2] == 2)
                {
                    x = Bodies[Deactivated[i][3]].nodes[Deactivated[i][4]].x_contact_force ;
                    y = Bodies[Deactivated[i][3]].nodes[Deactivated[i][4]].y_contact_force ;
                }
                else if (Deactivated[i][2] == 3)
                {
                    x = Bodies[Deactivated[i][3]].nodes[Deactivated[i][4]].x_body_force ;
                    y = Bodies[Deactivated[i][3]].nodes[Deactivated[i][4]].y_body_force ;
                }
                else if (Deactivated[i][2] == 4)
                {
                    x = Bodies[Deactivated[i][3]].nodes[Deactivated[i][4]].x_dirichlet_force ;
                    y = Bodies[Deactivated[i][3]].nodes[Deactivated[i][4]].y_dirichlet_force ;
                }
                else if (Deactivated[i][2] == 5)
                {
                    x = Bodies[Deactivated[i][3]].nodes[Deactivated[i][4]].x_neumann_force ;
                    y = Bodies[Deactivated[i][3]].nodes[Deactivated[i][4]].y_neumann_force ;
                }
                else if (Deactivated[i][2] == 6)
                {
                    x = Bodies[Deactivated[i][3]].nodes[Deactivated[i][4]].x_damping_force ;
                    y = Bodies[Deactivated[i][3]].nodes[Deactivated[i][4]].y_damping_force ;
                }
            }
            else
            {
                if (Bodies[Deactivated[i][2]].type == "deformable")
                {
                    x = 0. ;
                    y = 0. ;
                    for (int j(0) ; j<Bodies[Deactivated[i][2]].nb_nodes ; j++)
                    {
                        if (Deactivated[i][2] == 0)
                        {
                            x += Bodies[Deactivated[i][3]].nodes[j].x_force ;
                            y += Bodies[Deactivated[i][3]].nodes[j].y_force ;
                        }
                        else if (Deactivated[i][2] == 1)
                        {
                            x += Bodies[Deactivated[i][3]].nodes[j].x_internal_force ;
                            y += Bodies[Deactivated[i][3]].nodes[j].y_internal_force ;
                        }
                        else if (Deactivated[i][2] == 2)
                        {
                            x += Bodies[Deactivated[i][3]].nodes[j].x_contact_force ;
                            y += Bodies[Deactivated[i][3]].nodes[j].y_contact_force ;
                        }
                        else if (Deactivated[i][2] == 3)
                        {
                            x += Bodies[Deactivated[i][3]].nodes[j].x_body_force ;
                            y += Bodies[Deactivated[i][3]].nodes[j].y_body_force ;
                        }
                        else if (Deactivated[i][2] == 4)
                        {
                            x += Bodies[Deactivated[i][3]].nodes[j].x_dirichlet_force ;
                            y += Bodies[Deactivated[i][3]].nodes[j].y_dirichlet_force ;
                        }
                        else if (Deactivated[i][2] == 5)
                        {
                            x += Bodies[Deactivated[i][3]].nodes[j].x_neumann_force ;
                            y += Bodies[Deactivated[i][3]].nodes[j].y_neumann_force ;
                        }
                        else if (Deactivated[i][2] == 6)
                        {
                            x += Bodies[Deactivated[i][3]].nodes[j].x_damping_force ;
                            y += Bodies[Deactivated[i][3]].nodes[j].y_damping_force ;
                        }
                    }
                }
                else if (Bodies[Deactivated[i][2]].type == "rigid")
                {
                    if (Deactivated[i][2] == 0)
                    {
                        x = Bodies[Deactivated[i][3]].x_force ;
                        y = Bodies[Deactivated[i][3]].y_force ;
                        r = Bodies[Deactivated[i][3]].r_force ;
                    }
                    else if (Deactivated[i][2] == 2)
                    {
                        x = Bodies[Deactivated[i][3]].x_contact_force ;
                        y = Bodies[Deactivated[i][3]].y_contact_force ;
                        r = Bodies[Deactivated[i][3]].r_contact_force ;
                    }
                    else if (Deactivated[i][2] == 3)
                    {
                        x = Bodies[Deactivated[i][3]].x_body_force ;
                        y = Bodies[Deactivated[i][3]].y_body_force ;
                        r = Bodies[Deactivated[i][3]].r_body_force ;
                    }
                    else if (Deactivated[i][2] == 4)
                    {
                        x = Bodies[Deactivated[i][3]].x_dirichlet_force ;
                        y = Bodies[Deactivated[i][3]].y_dirichlet_force ;
                        r = Bodies[Deactivated[i][3]].r_dirichlet_force ;
                    }
                    else if (Deactivated[i][2] == 5)
                    {
                        x = Bodies[Deactivated[i][3]].x_neumann_force ;
                        y = Bodies[Deactivated[i][3]].y_neumann_force ;
                        r = Bodies[Deactivated[i][3]].r_neumann_force ;
                    }
                    else if (Deactivated[i][2] == 6)
                    {
                        x = Bodies[Deactivated[i][3]].x_damping_force ;
                        y = Bodies[Deactivated[i][3]].y_damping_force ;
                        r = Bodies[Deactivated[i][3]].r_damping_force ;
                    }
                }
            }
            if (Deactivated[i][1] == 0)
            {
                if ((Deactivated[i][5] == 0) & (x <= Deactivated[i][6]))
                    Bodies[Deactivated[i][3]].status = "inactive" ;
                else if ((Deactivated[i][5] == 1) & (x == Deactivated[i][6]))
                    Bodies[Deactivated[i][3]].status = "inactive" ;
                else if ((Deactivated[i][5] == 2) & (x >= Deactivated[i][6]))
                    Bodies[Deactivated[i][3]].status = "inactive" ;
            }
            else if (Deactivated[i][1] == 1)
            {
                if ((Deactivated[i][5] == 0) & (y <= Deactivated[i][6]))
                    Bodies[Deactivated[i][3]].status = "inactive" ;
                else if ((Deactivated[i][5] == 1) & (y == Deactivated[i][6]))
                    Bodies[Deactivated[i][3]].status = "inactive" ;
                else if ((Deactivated[i][5] == 2) & (y >= Deactivated[i][6]))
                    Bodies[Deactivated[i][3]].status = "inactive" ;
            }
            else if (Deactivated[i][1] == 2)
            {
                if ((Deactivated[i][5] == 0) & (pow(x*x+y+y,0.5) <= Deactivated[i][6]))
                    Bodies[Deactivated[i][3]].status = "inactive" ;
                else if ((Deactivated[i][5] == 1) & (pow(x*x+y+y,0.5) == Deactivated[i][6]))
                    Bodies[Deactivated[i][3]].status = "inactive" ;
                else if ((Deactivated[i][5] == 2) & (pow(x*x+y+y,0.5) >= Deactivated[i][6]))
                    Bodies[Deactivated[i][3]].status = "inactive" ;
            }
            else if (Deactivated[i][1] == 3)
            {
                if ((Deactivated[i][5] == 0) & (r <= Deactivated[i][6]))
                    Bodies[Deactivated[i][3]].status = "inactive" ;
                else if ((Deactivated[i][5] == 1) & (r == Deactivated[i][6]))
                    Bodies[Deactivated[i][3]].status = "inactive" ;
                else if ((Deactivated[i][5] == 2) & (r >= Deactivated[i][6]))
                    Bodies[Deactivated[i][3]].status = "inactive" ;
            }
        }
        else if (Deactivated[i][0] == 8)
        {
            //cout << Desactivated[i][0] << Desactivated[i][1] << endl ;
            Bodies[Deactivated[i][1]].Update_damage() ;
            //current_monitoring.push_back(Bodies[Desactivated[i][1]].damage) ;
            if ((Deactivated[i][2] == 0) & (Bodies[Deactivated[i][1]].damage <= Deactivated[i][3]))
                Bodies[Deactivated[i][1]].status = "inactive" ;
            else if ((Deactivated[i][2] == 1) & (Bodies[Deactivated[i][1]].damage == Deactivated[i][3]))
                Bodies[Deactivated[i][1]].status = "inactive" ;
            else if ((Deactivated[i][2] == 2) & (Bodies[Deactivated[i][1]].damage >= Deactivated[i][3]))
                Bodies[Deactivated[i][1]].status = "inactive" ;
        }
    }






}


//********************************************//
//** UPDATE CONTACT PRESSURES ****************//
//********************************************//

void Update_contact_pressures(
    int Nb_bodies,
    vector<Body>& Bodies )
{
    // Initializing contact pressures
    #pragma omp parallel
    {
        #pragma omp for schedule(dynamic)
        for (int i=0 ; i<Nb_bodies ; i++)
        {
            for (int j(0) ; j<Bodies[i].nb_borders ; j++)
            {
                for (int k(0) ; k<Bodies[i].borders[j].number_border_nodes ; k++)
                {
                    Bodies[i].borders[j].x_contact_pressure[k] = 0. ;
                    Bodies[i].borders[j].y_contact_pressure[k] = 0. ;
                }
            }
        }
    }

    // Computing contact pressures
    int bodyM, borderM, nodeM1, nodeM2, borderS, border_nodeS ;
    double fx, fy, invlength, shapeM1, shapeM2 ;
    // NB : possible conflict in memory writing, could only be parallelized with two-steps writing
    //#pragma omp parallel
    {
        //#pragma omp for schedule(dynamic)
        for (int i=0 ; i<Nb_bodies ; i++)
        {
            for (int j(0) ; j<Bodies[i].nb_contact_elements ; j++)
            {
                bodyM = Bodies[i].contact_elements[j].bodyM ;
                borderM = Bodies[i].contact_elements[j].borderM ;
                //nodeM1 = Bodies[i].contact_elements[j].nodeM1 ; // non, il faudrait border_node
                //nodeM2 = Bodies[i].contact_elements[j].nodeM2 ; // non, il faudrait border_node
                nodeM1 = Bodies[i].contact_elements[j].segmentM ;
                nodeM2 = Bodies[i].contact_elements[j].segmentM + 1 ;
                //
                borderS = Bodies[i].contact_elements[j].borderS ;
                border_nodeS = Bodies[i].contact_elements[j].border_nodeS ;
                fx = Bodies[i].contact_elements[j].fx ;
                fy = Bodies[i].contact_elements[j].fy ;
                invlength = 1. / Bodies[i].contact_elements[j].length ;
                shapeM1 = Bodies[i].contact_elements[j].shapeM1 ;
                shapeM2 = Bodies[i].contact_elements[j].shapeM2 ;
                Bodies[i].borders[borderS].x_contact_pressure[border_nodeS] += fx * invlength ;
                Bodies[i].borders[borderS].y_contact_pressure[border_nodeS] += fy * invlength ;
                Bodies[bodyM].borders[borderM].x_contact_pressure[nodeM1] -= fx * invlength * shapeM1 ;
                Bodies[bodyM].borders[borderM].y_contact_pressure[nodeM1] -= fy * invlength * shapeM1 ;
                Bodies[bodyM].borders[borderM].x_contact_pressure[nodeM2] -= fx * invlength * shapeM2 ;
                Bodies[bodyM].borders[borderM].y_contact_pressure[nodeM2] -= fy * invlength * shapeM2 ;
            }
        }
    }
}









#endif
