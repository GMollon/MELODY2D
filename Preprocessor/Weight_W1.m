function W=Weight_W1(Poi_Relative_Positions,List_Domain_Fnodes,Domains_Dimensions,Domain_Type)
%
%Cubic spline weight function
%
%Input variables :
%    Poi_Relative_Position : Poi relative positions with respect to each field node
%    List_Domain_Fnodes (Nv) : Field nodes used in the support domain
%    Domains_Dimensions (Ds) : x and y dimensions of the support domain of each field node
%
%Output variables :
%    W, with the following columns :
%        w
%        dwdx
%        dwdy
%
Number_Fnodes_Domain=size(List_Domain_Fnodes,1);
W=zeros(Number_Fnodes_Domain,3);
if Domain_Type==1
    for i=1:Number_Fnodes_Domain
        n=List_Domain_Fnodes(i);
        x=Poi_Relative_Positions(i,1)/Domains_Dimensions(n,1);
        y=Poi_Relative_Positions(i,2)/Domains_Dimensions(n,1);
        r=norm(Poi_Relative_Positions(i,:))/Domains_Dimensions(n,1);
        dr=1/Domains_Dimensions(n,1);
        if r>0.5 & r<=1
            w=(4/3)-4*r+4*r^2-4/3*r^3;
            dwdx=((-4*x^3*r+8*x^3-4*x*y^2*r+8*x*y^2-4*x*r)/(r^2))*dr;
            dwdy=((-4*y^3*r+8*y^3-4*y*x^2*r+8*y*x^2-4*y*r)/(r^2))*dr;
        elseif r<=0.5
            w=(2/3)-4*r^2+4*r^3;
            dwdx=(12*x*r-8*x)*dr;
            dwdy=(12*y*r-8*y)*dr;
        else
            w=0;
            dwdx=0;
            dwdy=0;
        end
        W(i,1)=w;
        W(i,2)=dwdx;
        W(i,3)=dwdy;
    end
elseif Domain_Type==2
    for i=1:Number_Fnodes_Domain
        n=List_Domain_Fnodes(i);
        if abs(Poi_Relative_Positions(i,1))==0
            drdx=1/Domains_Dimensions(n,1);
            rx=0;
        elseif Poi_Relative_Positions(i,1)<0
            drdx=sign(Poi_Relative_Positions(i,1))/Domains_Dimensions(n,1);
            rx=abs(Poi_Relative_Positions(i,1))/Domains_Dimensions(n,1);
        elseif Poi_Relative_Positions(i,1)>0
            drdx=sign(Poi_Relative_Positions(i,1))/Domains_Dimensions(n,2);
            rx=abs(Poi_Relative_Positions(i,1))/Domains_Dimensions(n,2);
        end
        if rx>0.5 & rx<=1
            wx=(4/3)-4*rx+4*rx^2-4/3*rx^3;
            dwxdx=(-4+8*rx-4*rx^2)*drdx;
        elseif rx<=0.5
            wx=(2/3)-4*rx^2+4*rx^3;
            dwxdx=(-8*rx+12*rx^2)*drdx;
        else
            wx=0;
            dwxdx=0;
        end

        if abs(Poi_Relative_Positions(i,2))==0
            drdy=1/Domains_Dimensions(n,2);
            ry=0;
        elseif Poi_Relative_Positions(i,2)<0
            drdy=sign(Poi_Relative_Positions(i,2))/Domains_Dimensions(n,3);
            ry=abs(Poi_Relative_Positions(i,2))/Domains_Dimensions(n,3);
        elseif Poi_Relative_Positions(i,2)>0
            drdy=sign(Poi_Relative_Positions(i,2))/Domains_Dimensions(n,4);
            ry=abs(Poi_Relative_Positions(i,2))/Domains_Dimensions(n,4);
        end
        if ry>0.5 & ry<=1
            wy=(4/3)-4*ry+4*ry^2-4/3*ry^3;
            dwydy=(-4+8*ry-4*ry^2)*drdy;
        elseif ry<=0.5
            wy=(2/3)-4*ry^2+4*ry^3;
            dwydy=(-8*ry+12*ry^2)*drdy;
        else
            wy=0;
            dwydy=0;
        end
        W(i,1)=wx*wy;
        W(i,2)=wy*dwxdx;
        W(i,3)=wx*dwydy;
    end
end