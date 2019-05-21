SigmaC=100000;     % Contrainte maximale caractéristique du système
dc=0.00001;        % Distance nodale minimale caractéristique du système
E0=1e6;         % Module d'Young maximal caractéristique du système
Rho=1000;      % Masse volumique minimale caractéristique du système

eps=SigmaC/E0;
Tv=dc*sqrt(Rho/(2*E0));
Alpha=16*Rho*dc^2/(SigmaC*Tv^2);
kn=Alpha*SigmaC/dc;

disp(' ')
disp(['Niveau de déformation caractéristique : ',num2str(eps)])
disp(['Période caractéristique d un élément de volume : ',num2str(Tv)])
disp(['Raideur kn pour avoir la même période aux contacts : ',num2str(kn)])
disp(['Ratio d interpénétration correspondant : ',num2str(Alpha)])
disp(' ')