function RBF=Gauss_RBF(Fnode_Position,Poi_Position,Domain_Size,param)
x=Poi_Position(1)-Fnode_Position(1);
y=Poi_Position(2)-Fnode_Position(2);
d=norm([x,y]);
D=Domain_Size;
if param==1
    if d<D
        RBF=exp(-(d*3/D)^2);
    else
        RBF=0;
    end
elseif param==3
    RBF=[0,0,0];
    if d<D
        RBF(1)=exp(-(d*3/D)^2);
        RBF(2)=-2*(3/D)^2*exp(-(d*3/D)^2)*x;
        RBF(3)=-2*(3/D)^2*exp(-(d*3/D)^2)*y;
    end
end