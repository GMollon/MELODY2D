function W=WeightW1(Poi_Relative_Positions,List_Domain_Fnodes,Domains_Dimensions)
Number_Fnodes_Domain=size(List_Domain_Fnodes,1);
W=zeros(Number_Fnodes_Domain,1);
for i=1:Number_Fnodes_Domain
    n=List_Domain_Fnodes(i);
    d=Domains_Dimensions(n,1);
    x=Poi_Relative_Positions(i,1);
    y=Poi_Relative_Positions(i,2);
    r=sqrt(x^2+y^2)/d;
    if r>0.5 & r<=1
        w=(4/3)-4*r+4*r^2-4/3*r^3;
    elseif r<=0.5
        w=(2/3)-4*r^2+4*r^3;
    else
        w=0;
    end
    W(i,1)=w;
end