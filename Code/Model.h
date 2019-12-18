#ifndef DEF_MODELS
#define DEF_MODELS


//********************************************//
//** APPLY ELASTIC LINEAR ********************//
//********************************************//

//void Apply_Elastic_Linear(vector<vector<double>>& Sig , vector<vector<double>>& Eps , double Mu , double Lambda , double& J , double& energy)
//{
//	J = 1. + Eps[0][0] + Eps[1][1] ;
//	Sig[0][0] = Lambda * (J - 1.) + 2. * Mu * Eps[0][0] ;
//	Sig[1][1] = Lambda * (J - 1.) + 2. * Mu * Eps[1][1] ;
//  Sig[0][1] = 2. * Mu * Eps[0][1] ;
//	Sig[1][0] = Sig[0][1] ;
//	Sig[2][2] = Lambda * (J - 1.) ;
//	energy = 0. ;
//}
void Apply_Elastic_Linear(double& S11 , double& S12 , double& S22 , double& S33 , double& E11 , double& E12 , double& E22 , double Mu , double Lambda , double& J , double& energy)
{
    // OPTIMIZED
	J = 1. + E11 + E22 ;
	S33 = Lambda * (J - 1.) ;
	Mu = 2 * Mu ;
	S11 = S33 + Mu * E11 ;
	S22 = S33 + Mu * E22 ;
    S12 = Mu * E12;
	energy = 0. ;
}



//********************************************//
//** APPLY NEOHOOKEAN ************************//
//********************************************//

//void Apply_NeoHookean(vector<vector<double>>& S , vector<vector<double>>& F , vector<vector<double>>& E , double Mu , double Kappa , double& J , double& energy)
//{
//	J = F[0][0] * F[1][1] - F[0][1] * F[1][0] ;
//	double MuJ53 = Mu * pow( J , -1.666666666667 ) ;
//	double trBover3 = 0.333333333333 * ( F[0][0] * F[0][0] + F[0][1] * F[0][1] + F[1][0] * F[1][0] + F[1][1] * F[1][1] + 1. ) ;
//	double sig11 = Kappa * (J - 1.) + MuJ53 * ( F[0][0] * F[0][0] + F[0][1] * F[0][1] - trBover3 ) ;
//	double sig22 = Kappa * (J - 1.) + MuJ53 * ( F[1][0] * F[1][0] + F[1][1] * F[1][1] - trBover3 ) ;
//	double sig12 = MuJ53 * ( F[0][0] * F[1][0] + F[0][1] * F[1][1] ) ;
//	S[0][0] = ( F[0][1] * ( F[0][1] * sig22 - F[1][1] * sig12 ) - F[1][1] * ( F[0][1] * sig12 - F[1][1] * sig11 ) ) / J ;
//	S[1][1] = ( F[0][0] * ( F[0][0] * sig22 - F[1][0] * sig12 ) - F[1][0] * ( F[0][0] * sig12 - F[1][0] * sig11 ) ) / J ;
//  S[0][1] = ( F[1][0] * ( F[0][1] * sig12 - F[1][1] * sig11 ) - F[0][0] * ( F[0][1] * sig22 - F[1][1] * sig12 ) ) / J ;
//	S[1][0] = S[0][1] ;
//	S[2][2] = ( Kappa * (J - 1.) + MuJ53 * ( 1. - trBover3 ) ) * J ;
//	energy = 0. ;
//}
void Apply_NeoHookean(double& S11 , double& S12 , double& S22 , double& S33 , double& F11 , double& F12 , double& F21 , double& F22 , double Mu , double Kappa , double& J , double& energy)
{
    // OPTIMIZED
    double F11F12 = F11 * F11 + F12 * F12 ;
    double F21F22 = F21 * F21 + F22 * F22 ;

	J = F11 * F22 - F12 * F21 ;
	double invJ = 1. / J ;
	double K = Kappa * (J - 1.) ;

	double MuJ53 = Mu * pow( J , -1.666666666667 ) ;
	double trBover3 = 0.333333333333 * ( F11F12 + F21F22 + 1. ) ;
	double sig11 = K + MuJ53 * ( F11F12 - trBover3 ) ;
	double sig22 = K + MuJ53 * ( F21F22 - trBover3 ) ;
	double sig12 = MuJ53 * ( F11 * F21 + F12 * F22 ) ;

	S11 = ( F12 * ( F12 * sig22 - F22 * sig12 ) - F22 * ( F12 * sig12 - F22 * sig11 ) ) * invJ ;
	S22 = ( F11 * ( F11 * sig22 - F21 * sig12 ) - F21 * ( F11 * sig12 - F21 * sig11 ) ) * invJ ;
    S12 = ( F21 * ( F12 * sig12 - F22 * sig11 ) - F11 * ( F12 * sig22 - F22 * sig12 ) ) * invJ ;
	S33 = ( K + MuJ53 * ( 1. - trBover3 ) ) * J ;
	energy = 0. ;
}



//********************************************//
//** APPLY FRICTIONLESS **********************//
//********************************************//

void Apply_Frictionless(double kn , double gapn , double& gapt , double& Pn , double& Pt)
{
	if (gapn>0.) Pn = 0. ;
	else Pn = -kn * gapn ;
	Pt = 0. ;
	gapt = 0. ;
}



//********************************************//
//** APPLY COHESION **************************//
//********************************************//

void Apply_Cohesion(double kn , double kt , double gtang , double gtens , double gapn , double& gapt , double& Pn , double& Pt)
{
	Pn = -kn * gapn ;
	Pt = -kt * gapt ;
	if (Pn<-gtens)
	{
		Pn = 0. ;
		Pt = 0. ;
		gapt = 0. ;
	}
	else if (abs(Pt)>gtang)
	{
		if (gapt>0.) Pt = -gtang ;
		else Pt = gtang ;
		gapt = -Pt / kt ;
	}
}



//********************************************//
//** APPLY MOHR COULOMB **********************//
//********************************************//

void Apply_Mohr_Coulomb(double kn , double kt , double fric , double coh , double tens , double gapn , double& gapt , double& Pn , double& Pt)
{
	Pn = -kn * gapn ;
	Pt = -kt * gapt ;
	if (Pn<-tens)
	{
		Pn = 0. ;
		Pt = 0. ;
		gapt = 0. ;
	}
	else if (abs(Pt)>coh+fric*Pn)
	{
		if (gapt>0.) Pt = -coh - fric * Pn ;
		else Pt = coh + fric * Pn ;
		gapt = -Pt / kt ;
	}
}



//********************************************//
//** APPLY DAMPED MOHR COULOMB ***************//
//********************************************//

void Apply_Damped_Mohr_Coulomb(double kn , double kt , double fric , double coh , double tens , double damp , double mass , double length , double gapn , double vgapn , double& gapt , double vgapt , double& Pn , double& Pt , double& W)
{
	Pn = 0. ;
	if (gapn<0.)
    {
        Pn = -damp * vgapn * 2. * sqrt(mass * kn / length) ;
        W += Pn * vgapn ;   // This is per unit time and per unit length : must be multiplied by length and deltat after returned
    }
    Pn -= kn * gapn ;
	if (Pn<-tens)
	{
		Pn = 0. ;
		Pt = 0. ;
		gapt = 0. ;
	}
	else
	{
	    double Pt_damp = - damp * vgapt * 2. * sqrt(mass * kt / length) ;
	    Pt = -kt * gapt + Pt_damp ;
	    if (abs(Pt)>coh+fric*Pn)
        {
            if (gapt>0.)    Pt = -coh - fric * Pn ;
            else            Pt = coh + fric * Pn ;
            gapt = -Pt / kt ;
            W += Pt * vgapt ;            // This is per unit time and per unit length : must be multiplied by length and deltat after returned
        }
        else    W += Pt_damp * vgapt ;   // This is per unit time and per unit length : must be multiplied by length and deltat after returned
	}
}

/*
void Apply_Damped_Mohr_Coulomb(double kn , double kt , double fric , double coh , double tens , double damp , double mass , double length , double gapn , double vgapn , double& gapt , double vgapt , double& Pn , double& Pt)
{
	Pn = -kn * gapn ;
	if (gapn<0.) Pn -= damp * vgapn * 2. * sqrt(mass * kn / length) ;
	if (Pn<-tens)
	{
		Pn = 0. ;
		Pt = 0. ;
		gapt = 0. ;
	}
	else
	{
	    Pt = -kt * gapt - damp * vgapt * 2. * sqrt(mass * kt / length) ;
	    if (abs(Pt)>coh+fric*Pn)
        {
            if (gapt>0.)    Pt = -coh - fric * Pn ;
            else            Pt = coh + fric * Pn ;
            gapt = -Pt / kt ;
        }
	}
}
*/



//************************************************//
//** APPLY TWO SLOPES MOHR COULOMB ***************//
//************************************************//

void Apply_Two_Slopes_Mohr_Coulomb(double kn , double kt , double fric , double coh , double tens , double rati , double gapn , double vgapn , double& gapt , double vgapt , double& Pn , double& Pt)
{
    if (vgapn <= 0.) Pn = -kn * gapn ;
    else             Pn = -kn * gapn * rati ;
	if (Pn < -tens)
	{
		Pn = 0. ;
		Pt = 0. ;
		gapt = 0. ;
	}
	else
	{
	    //if (vgapt * gapt > 0.) Pt = -kt * gapt ;
	    //else                   Pt = -kt * gapt * rati ;
	    Pt = -kt * gapt ;
	    //
        if (abs(Pt) > coh + fric * Pn)
        {
            if (gapt > 0.) Pt = -coh - fric * Pn ;
            else Pt = coh + fric * Pn ;
            //if (vgapt * gapt > 0.) gapt = -Pt / kt ;
            //else                   gapt = -Pt / (kt * rati) ;
            gapt = -Pt / kt ;
            //
        }
	}
}



//********************************************//
//** APPLY CZM LINEAR ************************//
//********************************************//

void Apply_CZM_linear(double kini , double Pnlim , double gapnlim , double Pnres , double NTratio , double gapn , double& gapt , double& Damage , double& Pn , double& Pt)
{
    if (Damage >= 1.)
    {
        // Bond already broken
        Pn = -kini * gapn ;
        Pt = -kini * gapt ;
        if (Pn<-Pnres)
        {
            Pn = 0. ;
            Pt = 0. ;
            gapt = 0. ;
        }
        else if (abs(Pt) > Pnres * NTratio)
        {
            if (gapt>0.) Pt = -Pnres * NTratio ;
            else Pt = Pnres * NTratio ;
            gapt = -Pt / kini ;
        }
    }
    else
    {
        double gapnmax = Pnlim / kini * ( 1. - Damage ) + gapnlim * Damage ;
        double gaptmax = gapnmax * NTratio ;
        double kmax ;
        if (gapn < gapnmax && abs(gapt) < gaptmax)
        {
            // Elastic part of the bond
            kmax = (Pnlim * (1. - Damage) + Pnres * Damage) / gapnmax ;
            if (gapn<0.)    Pn = -kini * gapn ;
            else            Pn = -kmax * gapn ;
            Pt = -kmax * gapt ;
        }
        else if (gapn > gapnmax && abs(gapt) < gaptmax)
        {
            // Damage in the normal direction only
            gapnmax = gapn ;
            gaptmax = gapnmax * NTratio ;
            Damage = (gapnmax - Pnlim / kini) / (gapnlim - Pnlim / kini) ;
            kmax = (Pnlim * (1. - Damage) + Pnres * Damage) / gapnmax ;
            if (gapn<0.)    Pn = -kini * gapn ;
            else            Pn = -kmax * gapn ;
            Pt = -kmax * gapt ;
        }
        else if (gapn < gapnmax && abs(gapt) > gaptmax)
        {
            // Damage in the tangential direction only
            gaptmax = abs(gapt) ;
            gapnmax = gaptmax / NTratio ;
            Damage = (gapnmax - Pnlim / kini) / (gapnlim - Pnlim / kini) ;
            kmax = (Pnlim * (1. - Damage) + Pnres * Damage) / gapnmax ;
            if (gapn<0.)    Pn = -kini * gapn ;
            else            Pn = -kmax * gapn ;
            Pt = -kmax * gapt ;
        }
        else if (gapn > gapnmax && abs(gapt) > gaptmax)
        {
            // Damage in both directions
            if (gapn/gapnmax > abs(gapt)/gaptmax)
            {
                gapnmax = gapn ;
                gaptmax = gapnmax * NTratio ;
            }
            else
            {
                gaptmax = abs(gapt) ;
                gapnmax = gaptmax / NTratio ;
            }
            Damage = (gapnmax - Pnlim / kini) / (gapnlim - Pnlim / kini) ;
            kmax = (Pnlim * (1. - Damage) + Pnres * Damage) / gapnmax ;
            if (gapn<0.)    Pn = -kini * gapn ;
            else            Pn = -kmax * gapn ;
            Pt = -kmax * gapt ;
        }
    }
}



//********************************************//
//** APPLY CZM FATIGUE ***********************//
//********************************************//

void Apply_CZM_fatigue(double kini , double Pnlim , double gapnres , double Pnres , double NTratio , double Pnfat, double Raten, double Ratet , double gapn , double& gapt , double& Damage , double& Pn , double& Pt)
{
    if (Damage < 1.)
    {
        double gapnmax = Pnlim / kini * ( 1. - Damage ) + gapnres * Damage ;
        double gaptmax = gapnmax * NTratio ;
        double kmax ;
        if (gapn < gapnmax && abs(gapt) < gaptmax)
        {
            // Elastic part of the bond
            double Pnm = Pn ;
            double Ptm = Pt ;
            kmax = (Pnlim * (1. - Damage) + Pnres * Damage) / gapnmax ;
            if (gapn<0.)    Pn = -kini * gapn ;
            else            Pn = -kmax * gapn ;
            Pt = -kmax * gapt ;
            double Pnmax = Pnlim * ( 1. - Damage ) + Pnres * Damage ;
            //double Pnfat = Degratio * Pnlim ;
            if ( ( abs(Pn) > Pnfat ) && ( abs(Pn) > abs(Pnm) ) )
            {
                double Ratiopn = ( abs(Pn) - Pnfat ) / ( Pnmax - Pnfat ) ;
                double Deltapn = abs(Pn) - abs(Pnm) ;
                Damage = Damage + Ratiopn * Raten * Deltapn ;
            }
            if ( ( abs(Pt) > NTratio * Pnfat ) && ( abs(Pt) > abs(Ptm) ) )
            {
                double Ratiopt = ( abs(Pt) - Pnfat * NTratio ) / ( NTratio * ( Pnmax - Pnfat ) ) ;
                double Deltapt = abs(Pt) - abs(Ptm) ;
                Damage = Damage + Ratiopt * Ratet * Deltapt ;
            }
        }
        else if (gapn > gapnmax && abs(gapt) < gaptmax)
        {
            // Damage in the normal direction only
            gapnmax = gapn ;
            gaptmax = gapnmax * NTratio ;
            Damage = (gapnmax - Pnlim / kini) / (gapnres - Pnlim / kini) ;
            kmax = (Pnlim * (1. - Damage) + Pnres * Damage) / gapnmax ;
            if (gapn<0.)    Pn = -kini * gapn ;
            else            Pn = -kmax * gapn ;
            Pt = -kmax * gapt ;
        }
        else if (gapn < gapnmax && abs(gapt) > gaptmax)
        {
            // Damage in the tangential direction only
            gaptmax = abs(gapt) ;
            gapnmax = gaptmax / NTratio ;
            Damage = (gapnmax - Pnlim / kini) / (gapnres - Pnlim / kini) ;
            kmax = (Pnlim * (1. - Damage) + Pnres * Damage) / gapnmax ;
            if (gapn<0.)    Pn = -kini * gapn ;
            else            Pn = -kmax * gapn ;
            Pt = -kmax * gapt ;
        }
        else if (gapn > gapnmax && abs(gapt) > gaptmax)
        {
            // Damage in both directions
            if (gapn/gapnmax > abs(gapt)/gaptmax)
            {
                gapnmax = gapn ;
                gaptmax = gapnmax * NTratio ;
            }
            else
            {
                gaptmax = abs(gapt) ;
                gapnmax = gaptmax / NTratio ;
            }
            Damage = (gapnmax - Pnlim / kini) / (gapnres - Pnlim / kini) ;
            kmax = (Pnlim * (1. - Damage) + Pnres * Damage) / gapnmax ;
            if (gapn<0.)    Pn = -kini * gapn ;
            else            Pn = -kmax * gapn ;
            Pt = -kmax * gapt ;
        }
    }
    if (Damage >= 1.)
    {
        // Bond broken
        Damage = 1. ;
        Pn = -kini * gapn ;
        Pt = -kini * gapt ;
        if (Pn<-Pnres)
        {
            Pn = 0. ;
            Pt = 0. ;
            gapt = 0. ;
        }
        else if (abs(Pt) > Pnres * NTratio)
        {
            if (gapt>0.) Pt = -Pnres * NTratio ;
            else Pt = Pnres * NTratio ;
            gapt = -Pt / kini ;
        }
    }
}



//********************************************//
//** APPLY BONDED MOHR COULOMB ***************//
//********************************************//

void Apply_Bonded_Mohr_Coulomb(double kn , double kt , double kbond , double fricbond , double cohbond , double tensbond , double fricfree , double cohfree , double tensfree , double damp , double mass , double length , double gapn , double vgapn , double& gapt , double vgapt , double& Damage , double& Pn , double& Pt)
{
    if (Damage >= 1.)
    {
        // Bond already broken
        Pn = -kn * gapn - damp * vgapn * 2. * sqrt(mass * kn / length) ;
        if (Pn<-tensfree)
        {
            Pn = 0. ;
            Pt = 0. ;
            gapt = 0. ;
        }
        else
        {
            Pt = -kt * gapt - damp * vgapt * 2. * sqrt(mass * kt / length) ;
        }
        if (abs(Pt)>cohfree+fricfree*Pn)
        {
            if (gapt>0.) Pt = -cohfree - fricfree * Pn ;
            else Pt = cohfree + fricfree * Pn ;
            gapt = -Pt / kt ;
        }
    }
    else
    {
        // Bond still active
        if (gapn<0.)    Pn = -kn * gapn - damp * vgapn * 2. * sqrt(mass * kn / length) ;
        else            Pn = -kbond * gapn - damp * vgapn * 2. * sqrt(mass * kn / length) ;
        Pt = -kbond * gapt - damp * vgapt * 2. * sqrt(mass * kt / length) ;
        if (Pn<-tensbond)
        {
            Damage = 1. ;
            if (Pn<-tensfree)
            {
                Pn = 0. ;
                Pt = 0. ;
                gapt = 0. ;
            }
        }
        if (abs(Pt)>cohbond+fricbond*Pn)
        {
            Damage = 1. ;
            if (gapt>0.) Pt = -cohfree - fricfree * Pn ;
            else Pt = cohfree + fricfree * Pn ;
            gapt = -Pt / kt ;
        }
    }
}






#endif
