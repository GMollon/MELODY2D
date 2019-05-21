#ifndef DEF_TOOLS
#define DEF_TOOLS

//********************************************//
//** C1 ELEMENT ******************************//
//********************************************//

void C1_element(double Xsi ,
				double X0 ,
				double Y0 ,
				double X1 ,
				double Y1 ,
				double X2 ,
				double Y2 ,
				double X3 ,
				double Y3 ,
				double& X ,
				double& Y ,
				double& Xp ,
				double& Yp ,
				double& Xpp ,
				double& Ypp ,
				double& N0 ,
				double& N1 ,
				double& N2 ,
				double& N3)
{
	N0 = - ( Xsi * Xsi - 1. ) * ( Xsi - 1. ) * 0.0625 ;
	N1 = ( 1. - Xsi ) * 0.5 + ( Xsi * Xsi - 1. ) * ( Xsi - 1. ) * 0.125 + ( Xsi * Xsi - 1. ) * ( Xsi + 1. ) * 0.0625 ;
	N2 = ( 1. + Xsi ) * 0.5 - ( Xsi * Xsi - 1. ) * ( Xsi - 1. ) * 0.0625 - ( Xsi * Xsi - 1. ) * ( Xsi + 1. ) * 0.125 ;
	N3 = ( Xsi * Xsi - 1. ) * ( Xsi + 1. ) * 0.0625 ;

	double N0p = - ( 3. * Xsi * Xsi - 2. * Xsi - 1. ) * 0.0625 ;
	double N1p = - 0.5 + ( 3. * Xsi * Xsi - 2. * Xsi - 1. ) * 0.125 + ( 3. * Xsi * Xsi + 2. * Xsi - 1. ) * 0.0625 ;
	double N2p = 0.5 - ( 3. * Xsi * Xsi - 2. * Xsi - 1. ) * 0.0625 - ( 3. * Xsi * Xsi + 2. * Xsi - 1. ) * 0.125 ;
	double N3p = ( 3. * Xsi * Xsi + 2. * Xsi - 1. ) * 0.0625 ;

	double N0pp = - ( 3. * Xsi - 1. ) * 0.125 ;
	double N1pp = ( 3. * Xsi - 1. ) * 0.25 + ( 3. * Xsi + 1. ) * 0.125 ;
	double N2pp = - ( 3. * Xsi - 1. ) * 0.125 - ( 3. * Xsi + 1. ) * 0.25 ;
	double N3pp = ( 3. * Xsi + 1. ) * 0.125 ;

	X = X0 * N0 + X1 * N1 + X2 * N2 + X3 * N3 ;
	Y = Y0 * N0 + Y1 * N1 + Y2 * N2 + Y3 * N3 ;
	Xp = X0 * N0p + X1 * N1p + X2 * N2p + X3 * N3p ;
	Yp = Y0 * N0p + Y1 * N1p + Y2 * N2p + Y3 * N3p ;
	Xpp = X0 * N0pp + X1 * N1pp + X2 * N2pp + X3 * N3pp ;
	Ypp = Y0 * N0pp + Y1 * N1pp + Y2 * N2pp + Y3 * N3pp ;
}


//********************************************//
//** CLOSEST 2 SEGMENT ***********************//
//********************************************//

void Closest_2_segment(string interpolant ,
                       double Xs ,
					   double Ys ,
					   double X0 ,
					   double Y0 ,
					   double X1 ,
					   double Y1 ,
					   double X2 ,
					   double Y2 ,
					   double X3 ,
					   double Y3 ,
					   double Tolerance ,
					   double& Gap ,
					   double& Xsi ,
					   double& Xnorm ,
					   double& Ynorm ,
					   double& Xtan ,
					   double& Ytan ,
					   double& Xclosest ,
					   double& Yclosest ,
					   double& N0 ,
					   double& N1 ,
					   double& N2 ,
					   double& N3 ,
					   int& flag_detect)
{
    if (interpolant == "Spline")
    {
        int max_step = 1000 ;
        double R ;
        int Step = 0 ;
        double Xp , Yp , Xpp , Ypp , Rp , invl ;
        do
        {
            Step++ ;
            C1_element(Xsi , X0 , Y0 , X1 , Y1 , X2 , Y2 , X3 , Y3 ,
                       Xclosest , Yclosest , Xp , Yp , Xpp , Ypp ,
                       N0 , N1 , N2 , N3 ) ;
            R = ( Xclosest - Xs ) * Xp + ( Yclosest - Ys ) * Yp ;
            Rp = Xp * Xp + ( Xclosest - Xs ) * Xpp + Yp * Yp + ( Yclosest - Ys ) * Ypp ;
            Xsi = Xsi - R / Rp ;
        } while ( (abs(R) > Tolerance) & (Step < max_step) ) ;
        invl = pow ( Xp * Xp + Yp * Yp , -0.5 ) ;
        Xtan = Xp * invl ;
        Ytan = Yp * invl ;
        Xnorm = Ytan ;
        Ynorm = -Xtan ;
        Gap = Xnorm * ( Xs - Xclosest) + Ynorm * ( Ys - Yclosest) ;
    }
    else if (interpolant == "Linear")
    {
        double ilength = pow( (X2-X1) * (X2-X1) + (Y2-Y1) * (Y2-Y1) , -0.5 ) ;
        double ilengthm = pow( (X1-X0) * (X1-X0) + (Y1-Y0) * (Y1-Y0) , -0.5 ) ;
        double ilengthp = pow( (X3-X2) * (X3-X2) + (Y3-Y2) * (Y3-Y2) , -0.5 ) ;
        Xtan = ( X2 - X1 ) * ilength ;
        Ytan = ( Y2 - Y1 ) * ilength ;
        Xnorm = Ytan ;
        Ynorm = -Xtan ;
        double Xtanm = ( X1 - X0 ) * ilengthm ;
        double Ytanm = ( Y1 - Y0 ) * ilengthm ;
        double Xtanp = ( X3 - X2 ) * ilengthp ;
        double Ytanp = ( Y3 - Y2 ) * ilengthp ;
        if (((Xs - X1) * (Xtanm + Xtan) + (Ys - Y1) * (Ytanm + Ytan)) * ( Xtan * (Xtanm + Xtan) + Ytan * (Ytanm + Ytan)) >= 0.)
        {
            if (((Xs - X2) * (Xtan + Xtanp) + (Ys - Y2) * (Ytan + Ytanp)) * ( Xtan * (Xtan + Xtanp) + Ytan * (Ytan + Ytanp)) <= 0.)
            {
                Xsi = -1. + 2. * ( (X2-X1) * (Xs-X1) + (Y2-Y1) * (Ys-Y1) ) * ilength * ilength ;
                if ( (Xsi<1.) && (Xsi>-1.) )
                {
                    Xnorm = Ytan ;
                    Ynorm = -Xtan ;
                    Xclosest = 0.5 * (X1+X2) + Xsi * 0.5 * (X2-X1) ;
                    Yclosest = 0.5 * (Y1+Y2) + Xsi * 0.5 * (Y2-Y1) ;
                    Gap = Xnorm * ( Xs - Xclosest) + Ynorm * ( Ys - Yclosest) ;
                    N0 = 0. ;
                    N1 = 0.5 * ( 1. - Xsi ) ;
                    N2 = 0.5 * ( 1. + Xsi ) ;
                    N3 = 0. ;
                    flag_detect = 0 ;
                }
                else if (Xsi==-1.)
                {
                    Xnorm = Ytan ;
                    Ynorm = -Xtan ;
                    Xclosest = X1 ;
                    Yclosest = Y1 ;
                    Gap = Xnorm * ( Xs - Xclosest) + Ynorm * ( Ys - Yclosest) ;
                    N0 = 0. ;
                    N1 = 1. ;
                    N2 = 0. ;
                    N3 = 0. ;
                    flag_detect = -1 ;
                }
                else if (Xsi<-1.)
                {
                    double invdistm = pow( (Xs-X1) * (Xs-X1) + (Ys-Y1) * (Ys-Y1) , -0.5 ) ;
                    Xnorm = ( Xs - X1 ) * invdistm ;
                    Ynorm = ( Ys - Y1 ) * invdistm ;
                    if ( -Ynorm * (Xtan + Xtanm) + Xnorm * (Ytan + Ytanm) < 0. )
                    {
                        Xnorm = -Xnorm ;
                        Ynorm = -Ynorm ;
                    }
                    Xtan = -Ynorm ;
                    Ytan = Xnorm ;
                    Xclosest = X1 ;
                    Yclosest = Y1 ;
                    Gap = Xnorm * ( Xs - Xclosest) + Ynorm * ( Ys - Yclosest) ;
                    N0 = 0. ;
                    N1 = 1. ;
                    N2 = 0. ;
                    N3 = 0. ;
                    Xsi = -1. ;
                    flag_detect = -2 ;
                }
                else if (Xsi==1.)
                {
                    Xnorm = Ytan ;
                    Ynorm = -Xtan ;
                    Xclosest = X2 ;
                    Yclosest = Y2 ;
                    Gap = Xnorm * ( Xs - Xclosest) + Ynorm * ( Ys - Yclosest) ;
                    N0 = 0. ;
                    N1 = 0. ;
                    N2 = 1. ;
                    N3 = 0. ;
                    flag_detect = 1 ;
                }
                else if (Xsi>1.)
                {
                    double invdistp = pow( (Xs-X2) * (Xs-X2) + (Ys-Y2) * (Ys-Y2) , -0.5 ) ;
                    Xnorm = ( Xs - X2 ) * invdistp ;
                    Ynorm = ( Ys - Y2 ) * invdistp ;
                    if ( -Ynorm * (Xtan + Xtanp) + Xnorm * (Ytan + Ytanp) < 0. )
                    {
                        Xnorm = -Xnorm ;
                        Ynorm = -Ynorm ;
                    }
                    Xtan = -Ynorm ;
                    Ytan = Xnorm ;
                    Xclosest = X2 ;
                    Yclosest = Y2 ;
                    Gap = Xnorm * ( Xs - Xclosest) + Ynorm * ( Ys - Yclosest) ;
                    N0 = 0. ;
                    N1 = 0. ;
                    N2 = 1. ;
                    N3 = 0. ;
                    Xsi = 1. ;
                    flag_detect = 2 ;
                }
            }
            else
            {
                flag_detect = 3 ;
                Xsi = 1. ;
            }
        }
        else
        {
            flag_detect = -3 ;
            Xsi = -1. ;
        }


        /*
        Xsi = -1. + 2. * ( (X2-X1) * (Xs-X1) + (Y2-Y1) * (Ys-Y1) ) * ilength * ilength ;
        if ( (Xsi<=1.) && (Xsi>=-1.) )
        {
            Xnorm = Ytan ;
            Ynorm = -Xtan ;
            Xclosest = 0.5 * (X1+X2) + Xsi * 0.5 * (X2-X1) ;
            Yclosest = 0.5 * (Y1+Y2) + Xsi * 0.5 * (Y2-Y1) ;
            Gap = Xnorm * ( Xs - Xclosest) + Ynorm * ( Ys - Yclosest) ;
            N0 = 0. ;
            N1 = 0.5 * ( 1. - Xsi ) ;
            N2 = 0.5 * ( 1. + Xsi ) ;
            N3 = 0. ;
            //
            if ( Gap < 0. )
            {
                double ilengthm = pow( (X1-X0) * (X1-X0) + (Y1-Y0) * (Y1-Y0) , -0.5 ) ;
                if ( Gap < -( Xsi + 1. ) * 0.5  * ( (X1-X0)*ilengthm * (X2-X1)*ilength + (Y1-Y0)*ilengthm * (Y2-Y1)*ilength + 1. ) / ( ilength * ( (X1-X0)*ilengthm * (Y2-Y1)*ilength + (Y1-Y0)*ilengthm * (X1-X2)*ilength ) ) )
                {
                    Xsi = -2. ;
                    return ;
                }
                double ilengthp = pow( (X3-X2) * (X3-X2) + (Y3-Y2) * (Y3-Y2) , -0.5 ) ;
                if ( Gap < -( 1. - Xsi ) * 0.5  * ( (X3-X2)*ilengthp * (X2-X1)*ilength + (Y3-Y2)*ilengthp * (Y2-Y1)*ilength + 1. ) / ( ilength * ( (X3-X2)*ilengthp * (Y2-Y1)*ilength + (Y3-Y2)*ilengthp * (X1-X2)*ilength ) ) )
                {
                    Xsi = 2. ;
                    return ;
                }
            }
            //
        }
        else if (Xsi>1.)
        {
            if ( -1. + 2. * ( (X3-X2) * (Xs-X2) + (Y3-Y2) * (Ys-Y2) ) / ( (X3-X2) * (X3-X2) + (Y3-Y2) * (Y3-Y2) ) < -1. )
            {
                double ilengthp = pow( (X3-X2) * (X3-X2) + (Y3-Y2) * (Y3-Y2) , -0.5 ) ;
                double invdistp = pow( (Xs-X2) * (Xs-X2) + (Ys-Y2) * (Ys-Y2) , -0.5 ) ;
                double Xtanp = ( X3 - X2 ) * ilengthp ;
                double Ytanp = ( Y3 - Y2 ) * ilengthp ;
                if ( (Xs-X2)*Xtanp+(Ys-Y2)*Ytanp < -(Xs-X2)*Xtan-(Ys-Y2)*Ytan ) Xsi = 0.999999999999 ;
                else Xsi = 1.000000000001 ;
                Xnorm = ( Xs - X2 ) * invdistp ;
                Ynorm = ( Ys - Y2 ) * invdistp ;
                if ( -Ynorm * (Xtan + Xtanp) + Xnorm * (Ytan + Ytanp) < 0. )
                {
                    Xnorm = -Xnorm ;
                    Ynorm = -Ynorm ;
                }
                Xtan = -Ynorm ;
                Ytan = Xnorm ;
                Xclosest = X2 ;
                Yclosest = Y2 ;
                Gap = Xnorm * ( Xs - Xclosest) + Ynorm * ( Ys - Yclosest) ;
                N0 = 0. ;
                N1 = 0. ;
                N2 = 1. ;
                N3 = 0. ;
            }
        }
        else if (Xsi<-1.)
        {
            if ( -1. + 2. * ( (X1-X0) * (Xs-X0) + (Y1-Y0) * (Ys-Y0) ) / ( (X1-X0) * (X1-X0) + (Y1-Y0) * (Y1-Y0) ) > 1. )
            {
                double ilengthm = pow( (X1-X0) * (X1-X0) + (Y1-Y0) * (Y1-Y0) , -0.5 ) ;
                double invdistm = pow( (Xs-X1) * (Xs-X1) + (Ys-Y1) * (Ys-Y1) , -0.5 ) ;
                double Xtanm = ( X1 - X0 ) * ilengthm ;
                double Ytanm = ( Y1 - Y0 ) * ilengthm ;
                if ( -(Xs-X1)*Xtanm-(Ys-Y1)*Ytanm < (Xs-X1)*Xtan+(Ys-Y1)*Ytan ) Xsi = -0.999999999999 ;
                else Xsi = -1.000000000001 ;
                Xnorm = ( Xs - X1 ) * invdistm ;
                Ynorm = ( Ys - Y1 ) * invdistm ;
                if ( -Ynorm * (Xtan + Xtanm) + Xnorm * (Ytan + Ytanm) < 0. )
                {
                    Xnorm = -Xnorm ;
                    Ynorm = -Ynorm ;
                }
                Xtan = -Ynorm ;
                Ytan = Xnorm ;
                Xclosest = X1 ;
                Yclosest = Y1 ;
                Gap = Xnorm * ( Xs - Xclosest) + Ynorm * ( Ys - Yclosest) ;
                N0 = 0. ;
                N1 = 1. ;
                N2 = 0. ;
                N3 = 0. ;
            }
        }
        */
    }
}




//********************************************//
//** CLOSEST 2 SEGMENT SELF ******************//
//********************************************//

void Closest_2_segment_self(double Xs ,
					       double Ys ,
					       double X1 ,
					       double Y1 ,
					       double X2 ,
					       double Y2 ,
					       double Tolerance ,
					       double& Gap ,
					       double& Xsi ,
					       double& Xnorm ,
					       double& Ynorm ,
					       double& Xtan ,
					       double& Ytan ,
					       double& Xclosest ,
					       double& Yclosest ,
					       double& N0 ,
					       double& N1 ,
					       double& N2 ,
					       double& N3 ,
                           int& flag_detect)
{
    double ilength = pow( (X2-X1) * (X2-X1) + (Y2-Y1) * (Y2-Y1) , -0.5 ) ;
    Xtan = ( X2 - X1 ) * ilength ;
    Ytan = ( Y2 - Y1 ) * ilength ;
    if (((Xs - X1) * Xtan + (Ys - Y1) * Ytan) * (Xtan * Xtan + Ytan * Ytan) >= 0.)
    {
        if (((Xs - X2) * Xtan + (Ys - Y2) * Ytan) * (Xtan * Xtan + Ytan * Ytan) <= 0.)
        {
            Xsi = -1. + 2. * ( (X2-X1) * (Xs-X1) + (Y2-Y1) * (Ys-Y1) ) * ilength * ilength ;
            if ( (Xsi<1.) && (Xsi>-1.) )
            {
                flag_detect = 0 ;
                Xnorm = Ytan ;
                Ynorm = -Xtan ;
                Xclosest = 0.5 * (X1+X2) + Xsi * 0.5 * (X2-X1) ;
                Yclosest = 0.5 * (Y1+Y2) + Xsi * 0.5 * (Y2-Y1) ;
                Gap = Xnorm * ( Xs - Xclosest) + Ynorm * ( Ys - Yclosest) ;
                N0 = 0. ;
                N1 = 0.5 * ( 1. - Xsi ) ;
                N2 = 0.5 * ( 1. + Xsi ) ;
                N3 = 0. ;
            }
        }
        else Xsi = 1. ;
        flag_detect = 3 ;
    }
    else Xsi = -1. ;
    flag_detect = -3 ;
}





//********************************************//
//** TRIANGLE SURFACE ************************//
//********************************************//

double Triangle_surface(double x0 , double y0 , double x1 , double y1 , double x2 , double y2)
{
	return 0.5 * abs( (x1-x0) * (y2-y0) - (x2-x0) * (y1-y0) ) ;
}



//********************************************//
//** TRIANGLE INERTIA ************************//
//********************************************//

double Triangle_inertia(double x0 , double y0 , double x1 , double y1 , double x2 , double y2 , double xc , double yc)
{
    double area , PP , PQ , QQ , inertia(0) ;

    area = Triangle_surface( xc , yc , x0 , y0 , x1 , y1 ) ;
    PP = (x0-xc) * (x0-xc) + (y0-yc) * (y0-yc) ;
    PQ = (x0-xc) * (x1-xc) + (y0-yc) * (y1-yc) ;
    QQ = (x1-xc) * (x1-xc) + (y1-yc) * (y1-yc) ;
    inertia += area * ( PP + PQ + QQ ) ;

    area = Triangle_surface( xc , yc , x1 , y1 , x2 , y2 ) ;
    PP = (x1-xc) * (x1-xc) + (y1-yc) * (y1-yc) ;
    PQ = (x1-xc) * (x2-xc) + (y1-yc) * (y2-yc) ;
    QQ = (x2-xc) * (x2-xc) + (y2-yc) * (y2-yc) ;
    inertia += area * ( PP + PQ + QQ ) ;

    area = Triangle_surface( xc , yc , x2 , y2 , x0 , y0 ) ;
    PP = (x2-xc) * (x2-xc) + (y2-yc) * (y2-yc) ;
    PQ = (x2-xc) * (x0-xc) + (y2-yc) * (y0-yc) ;
    QQ = (x0-xc) * (x0-xc) + (y0-yc) * (y0-yc) ;
    inertia += area * ( PP + PQ + QQ ) ;

	return abs( inertia / 6. ) ;
}






#endif
