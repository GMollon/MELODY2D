function Inertia=TriangleInertia(X1,Y1,X2,Y2,X3,Y3,Xc,Yc,Rho)
%
% Function returns the moment of inertia of the tiangle 123
% of density rho with respect to its centre of mass c.

Area=polyarea([Xc,X1,X2,Xc],[Yc,Y1,Y2,Yc]);
PP=(X1-Xc)*(X1-Xc)+(Y1-Yc)*(Y1-Yc);
PQ=(X1-Xc)*(X2-Xc)+(Y1-Yc)*(Y2-Yc);
QQ=(X2-Xc)*(X2-Xc)+(Y2-Yc)*(Y2-Yc);
Inertia=Rho*Area/6*(PP+PQ+QQ);

Area=polyarea([Xc,X2,X3,Xc],[Yc,Y2,Y3,Yc]);
PP=(X2-Xc)*(X2-Xc)+(Y2-Yc)*(Y2-Yc);
PQ=(X2-Xc)*(X3-Xc)+(Y2-Yc)*(Y3-Yc);
QQ=(X3-Xc)*(X3-Xc)+(Y3-Yc)*(Y3-Yc);
Inertia=Inertia+Rho*Area/6*(PP+PQ+QQ);

Area=polyarea([Xc,X3,X1,Xc],[Yc,Y3,Y1,Yc]);
PP=(X3-Xc)*(X3-Xc)+(Y3-Yc)*(Y3-Yc);
PQ=(X3-Xc)*(X1-Xc)+(Y3-Yc)*(Y1-Yc);
QQ=(X1-Xc)*(X1-Xc)+(Y1-Yc)*(Y1-Yc);
Inertia=Inertia+Rho*Area/6*(PP+PQ+QQ);

Inertia=abs(Inertia);

