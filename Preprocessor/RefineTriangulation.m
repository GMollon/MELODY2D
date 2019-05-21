function [Tri,Pos]=RefineTriangulation(Triangulation,Positions)
Tri=zeros(4*size(Triangulation,1),3);
Pos=zeros(6*size(Triangulation,1),2);
for Facet=1:size(Triangulation,1)
    P1=Positions(Triangulation(Facet,1),:);
    P2=Positions(Triangulation(Facet,2),:);
    P3=Positions(Triangulation(Facet,3),:);
    Tri(4*(Facet-1)+1,:)=[6*(Facet-1)+1,6*(Facet-1)+4,6*(Facet-1)+6];
    Tri(4*(Facet-1)+2,:)=[6*(Facet-1)+4,6*(Facet-1)+2,6*(Facet-1)+5];
    Tri(4*(Facet-1)+3,:)=[6*(Facet-1)+6,6*(Facet-1)+5,6*(Facet-1)+3];
    Tri(4*(Facet-1)+4,:)=[6*(Facet-1)+4,6*(Facet-1)+5,6*(Facet-1)+6];
    Pos(6*(Facet-1)+1,:)=P1;
    Pos(6*(Facet-1)+2,:)=P2;
    Pos(6*(Facet-1)+3,:)=P3;
    Pos(6*(Facet-1)+4,:)=(P1+P2)/2;
    Pos(6*(Facet-1)+5,:)=(P2+P3)/2;
    Pos(6*(Facet-1)+6,:)=(P3+P1)/2;
end