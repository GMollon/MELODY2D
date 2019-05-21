function Shape_Functions=T3FEM(Poi_Position,Fnodes_Positions,List_Domain_Fnodes)
X1=Fnodes_Positions(List_Domain_Fnodes(1),1);
Y1=Fnodes_Positions(List_Domain_Fnodes(1),2);
X2=Fnodes_Positions(List_Domain_Fnodes(2),1);
Y2=Fnodes_Positions(List_Domain_Fnodes(2),2);
X3=Fnodes_Positions(List_Domain_Fnodes(3),1);
Y3=Fnodes_Positions(List_Domain_Fnodes(3),2);
X=Poi_Position(1);
Y=Poi_Position(2);

N1=cross([X2-X1,Y2-Y1,0-1],[X3-X1,Y3-Y1,0-1]);
N1=N1/norm(N1);
N2=cross([X2-X1,Y2-Y1,1-0],[X3-X1,Y3-Y1,0-0]);
N2=N2/norm(N2);
N3=cross([X2-X1,Y2-Y1,0-0],[X3-X1,Y3-Y1,1-0]);
N3=N3/norm(N3);

Phi1=(X1*N1(1)+Y1*N1(2)+1*N1(3)-X*N1(1)-Y*N1(2))/N1(3);
Phi2=(X2*N2(1)+Y2*N2(2)+1*N2(3)-X*N2(1)-Y*N2(2))/N2(3);
Phi3=(X3*N3(1)+Y3*N3(2)+1*N3(3)-X*N3(1)-Y*N3(2))/N3(3);

Shape_Functions=[Phi1,Phi2,Phi3;-N1(1)/N1(3),-N2(1)/N2(3),-N3(1)/N3(3);-N1(2)/N1(3),-N2(2)/N2(3),-N3(2)/N3(3)]';