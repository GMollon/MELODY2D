SigmaC=100000;     % Contrainte maximale caract�ristique du syst�me
dc=0.00001;        % Distance nodale minimale caract�ristique du syst�me
E0=1e6;         % Module d'Young maximal caract�ristique du syst�me
Rho=1000;      % Masse volumique minimale caract�ristique du syst�me

eps=SigmaC/E0;
Tv=dc*sqrt(Rho/(2*E0));
Alpha=16*Rho*dc^2/(SigmaC*Tv^2);
kn=Alpha*SigmaC/dc;

disp(' ')
disp(['Niveau de d�formation caract�ristique : ',num2str(eps)])
disp(['P�riode caract�ristique d un �l�ment de volume : ',num2str(Tv)])
disp(['Raideur kn pour avoir la m�me p�riode aux contacts : ',num2str(kn)])
disp(['Ratio d interp�n�tration correspondant : ',num2str(Alpha)])
disp(' ')