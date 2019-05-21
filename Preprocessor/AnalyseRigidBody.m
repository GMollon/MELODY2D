function [MASSES,INVERSE_MASSES,CENTRE_OF_MASS]=AnalyseRigidBody(NUMBER_NODES,INITIAL_POSITIONS,NUMBER_CELLS,TRIANGULATION,MATERIAL,Materials)
CENTRE_OF_MASS=[0,0];
MASSES=[0,0];

for i=1:size(Materials,1)
    if strcmp(Materials{i,1},MATERIAL)
        Material_Number=i;
        break
    end
end
Material_Type=Materials{Material_Number,2};
Material_Parameters=Materials{Material_Number,3};
Rho=Material_Parameters(1,1);

for i=1:NUMBER_CELLS
    x1=INITIAL_POSITIONS(TRIANGULATION(i,1),1);
    y1=INITIAL_POSITIONS(TRIANGULATION(i,1),2);
    x2=INITIAL_POSITIONS(TRIANGULATION(i,2),1);
    y2=INITIAL_POSITIONS(TRIANGULATION(i,2),2);
    x3=INITIAL_POSITIONS(TRIANGULATION(i,3),1);
    y3=INITIAL_POSITIONS(TRIANGULATION(i,3),2);
    Area=polyarea([x1,x2,x3,x1],[y1,y2,y3,y1]);
    xc=(x1+x2+x3)/3;
    yc=(y1+y2+y3)/3;
    MASSES(1)=MASSES(1)+Area*Rho;
    CENTRE_OF_MASS(1)=CENTRE_OF_MASS(1)+Area*Rho*xc;
    CENTRE_OF_MASS(2)=CENTRE_OF_MASS(2)+Area*Rho*yc;
end
CENTRE_OF_MASS=CENTRE_OF_MASS/MASSES(1);

for i=1:NUMBER_CELLS
    x1=INITIAL_POSITIONS(TRIANGULATION(i,1),1);
    y1=INITIAL_POSITIONS(TRIANGULATION(i,1),2);
    x2=INITIAL_POSITIONS(TRIANGULATION(i,2),1);
    y2=INITIAL_POSITIONS(TRIANGULATION(i,2),2);
    x3=INITIAL_POSITIONS(TRIANGULATION(i,3),1);
    y3=INITIAL_POSITIONS(TRIANGULATION(i,3),2);
    Area=polyarea([x1,x2,x3,x1],[y1,y2,y3,y1]);
    xc=(x1+x2+x3)/3;
    yc=(y1+y2+y3)/3;
    MASSES(2)=MASSES(2)+Area*Rho*((xc-CENTRE_OF_MASS(1))^2+(yc-CENTRE_OF_MASS(2))^2);
    MASSES(2)=MASSES(2)+TriangleInertia(x1,y1,x2,y2,x3,y3,xc,yc,Rho);
    
end

INVERSE_MASSES=1./MASSES;