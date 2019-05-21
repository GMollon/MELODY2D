function [INITIAL_POSITIONS,BORDERS,TRIANGULATION]=SetUnstructuredNodesForHoledContour(Contour,ContourIn,Interpolant,Average_Nodal_Distance,Mesh_Ratio,Activate_Plot)

% Close the contours if needed
if norm(Contour(1,:)-Contour(size(Contour,1),:))>1e-3*Average_Nodal_Distance
    Contour=[Contour;Contour(1,:)];
else
    Contour(size(Contour,1),:)=Contour(1,:);
end
if norm(ContourIn(1,:)-ContourIn(size(ContourIn,1),:))>1e-3*Average_Nodal_Distance
    ContourIn=[ContourIn;ContourIn(1,:)];
else
    ContourIn(size(ContourIn,1),:)=ContourIn(1,:);
end

% Compute regular boundary nodes
Ltot=0;
Lengths=zeros(size(Contour,1)+1,1);
for j=1:size(Contour,1)-1
    Ltot=Ltot+sqrt((Contour(j+1,1)-Contour(j,1))^2+(Contour(j+1,2)-Contour(j,2))^2);
    Lengths(j+1,1)=Ltot;
end
Ltot=Ltot+sqrt((Contour(j+1,1)-Contour(1,1))^2+(Contour(j+1,2)-Contour(1,2))^2);
Lengths(j+2,1)=Ltot;
Ninterp=ceil(Ltot/Average_Nodal_Distance);
DistInterp=Ltot/Ninterp;
InterpContour=zeros(Ninterp,2);
InterpContour(1,:)=Contour(1,1:2);
for i=2:Ninterp
    D=(i-1)*DistInterp;
    for j=1:size(Lengths,1)
        if D>=Lengths(j,1) & D<=Lengths(j+1,1)
            ratio=(D-Lengths(j,1))/(Lengths(j+1,1)-Lengths(j,1));
            InterpContour(i,:)=(1-ratio)*Contour(j,1:2)+ratio*Contour(j+1,1:2);
            break
        end
    end
end
Ltot=0;
LengthsIn=zeros(size(ContourIn,1)+1,1);
for j=1:size(ContourIn,1)-1
    Ltot=Ltot+sqrt((ContourIn(j+1,1)-ContourIn(j,1))^2+(ContourIn(j+1,2)-ContourIn(j,2))^2);
    LengthsIn(j+1,1)=Ltot;
end
Ltot=Ltot+sqrt((ContourIn(j+1,1)-ContourIn(1,1))^2+(ContourIn(j+1,2)-ContourIn(1,2))^2);
LengthsIn(j+2,1)=Ltot;
Ninterp=ceil(Ltot/Average_Nodal_Distance);
DistInterp=Ltot/Ninterp;
InterpContourIn=zeros(Ninterp,2);
InterpContourIn(1,:)=ContourIn(1,1:2);
for i=2:Ninterp
    D=(i-1)*DistInterp;
    for j=1:size(Lengths,1)
        if D>=LengthsIn(j,1) & D<=LengthsIn(j+1,1)
            ratio=(D-LengthsIn(j,1))/(LengthsIn(j+1,1)-LengthsIn(j,1));
            InterpContourIn(i,:)=(1-ratio)*ContourIn(j,1:2)+ratio*ContourIn(j+1,1:2);
            break
        end
    end
end

% Position the field nodes and triangulate
Box=[min(InterpContour(:,1)),min(InterpContour(:,2));max(InterpContour(:,1)),max(InterpContour(:,2))];
%
xc=mean(InterpContour(:,1));
yc=mean(InterpContour(:,2));
rmax=max(((InterpContour(:,1)-xc).^2+(InterpContour(:,2)-yc).^2).^0.5);
pmin=0;
pmax=rmax;
fmin=Average_Nodal_Distance*Mesh_Ratio(1);
fmax=Average_Nodal_Distance*Mesh_Ratio(2);
a=(fmax-fmin)/rmax;
b=fmin;
fh=@(p) b+a*((p(:,1)-xc).^2+(p(:,2)-yc).^2).^0.5;
fd=@(p) ddiff(dpoly(p,InterpContour),dpoly(p,InterpContourIn));
%[INITIAL_POSITIONS,TRIANGULATION]=distmesh2d_plot(Activate_Plot,0.01,@dpoly,fh,Average_Nodal_Distance,Box,InterpContour,InterpContour);
%[INITIAL_POSITIONS,TRIANGULATION]=distmesh2d_plot(Activate_Plot,0.01,@dpoly,@huniform,Average_Nodal_Distance,Box,InterpContour,InterpContour);
[INITIAL_POSITIONS,TRIANGULATION]=distmesh2d_plot(Activate_Plot,0.01,fd,@huniform,Average_Nodal_Distance,Box,[InterpContour;InterpContourIn]);
close

Flag_Correct=0;
for i=2:size(INITIAL_POSITIONS,1)
    if norm(INITIAL_POSITIONS(i,:)-INITIAL_POSITIONS(i-1,:))<1e-3*Average_Nodal_Distance
        Flag_Correct=1;
        break
    end
end
if Flag_Correct==1
    INITIAL_POSITIONS=[INITIAL_POSITIONS(1:i-1,:);INITIAL_POSITIONS(i+1:size(INITIAL_POSITIONS,1),:)];
    for j=1:size(TRIANGULATION,1)
        if TRIANGULATION(j,1)>=i
            TRIANGULATION(j,1)=TRIANGULATION(j,1)-1;
        end
        if TRIANGULATION(j,2)>=i
            TRIANGULATION(j,2)=TRIANGULATION(j,2)-1;
        end
        if TRIANGULATION(j,3)>=i
            TRIANGULATION(j,3)=TRIANGULATION(j,3)-1;
        end
    end
end

% Fix the orientation of the triangular facets
for i=1:size(TRIANGULATION,1)
    P1=INITIAL_POSITIONS(TRIANGULATION(i,1),:);
    P2=INITIAL_POSITIONS(TRIANGULATION(i,2),:);
    P3=INITIAL_POSITIONS(TRIANGULATION(i,3),:);
    u=P2-P1;
    v=(P1+P2+P3)/3-P1;
    if u(1)*v(2)-u(2)*v(1)<0
        TRIANGULATION(i,:)=flipud(TRIANGULATION(i,:)')';
    end
end

% Detect and sort the boundary nodes
Segments=zeros(3*size(TRIANGULATION,1),3);
for i=1:size(TRIANGULATION,1)
    Segments(3*(i-1)+1,1:2)=TRIANGULATION(i,[1,2]);
    if TRIANGULATION(i,1)<TRIANGULATION(i,2)
        Segments(3*(i-1)+1,3)=1;
    else
        Segments(3*(i-1)+1,3)=2;
    end
    Segments(3*(i-1)+2,1:2)=TRIANGULATION(i,[2,3]);
    if TRIANGULATION(i,2)<TRIANGULATION(i,3)
        Segments(3*(i-1)+2,3)=1;
    else
        Segments(3*(i-1)+2,3)=2;
    end
    Segments(3*(i-1)+3,1:2)=TRIANGULATION(i,[3,1]);
    if TRIANGULATION(i,3)<TRIANGULATION(i,1)
        Segments(3*(i-1)+3,3)=1;
    else
        Segments(3*(i-1)+3,3)=2;
    end
end
Segments=sortrows(cat(2,sort(Segments(:,1:2)')',Segments(:,3)),[1,2]);
Toremove=[];
for i=2:size(Segments,1)
    if Segments(i,1)==Segments(i-1,1) & Segments(i,2)==Segments(i-1,2)
        Toremove=cat(1,Toremove,[i-1;i]);
    end
end
Unsorted=zeros(size(Segments,1)-size(Toremove,1),2);
nb=1;
ntr=1;
for i=1:size(Segments,1)
    if i==Toremove(ntr,1)
        if ntr<size(Toremove,1)
            ntr=ntr+1;
        end
    else
        if Segments(i,3)==1
            Unsorted(nb,:)=Segments(i,[1,2]);
        elseif Segments(i,3)==2
            Unsorted(nb,:)=Segments(i,[2,1]);
        end
        nb=nb+1;
    end
end

Sorted=SortSegments(Unsorted);
if polyarea(INITIAL_POSITIONS(Sorted{1,1},1),INITIAL_POSITIONS(Sorted{1,1},2))>polyarea(INITIAL_POSITIONS(Sorted{2,1},1),INITIAL_POSITIONS(Sorted{2,1},2))
    BORDERS={'Closed',size(Sorted{1,1},1),Sorted{1,1},Interpolant;'Closed',size(Sorted{2,1},1),Sorted{2,1},Interpolant};
else
    BORDERS={'Closed',size(Sorted{2,1},1),Sorted{2,1},Interpolant;'Closed',size(Sorted{1,1},1),Sorted{1,1},Interpolant};
end







return
BORDERS=zeros(size(Unsorted,1),2);
BORDERS(1,:)=Unsorted(1,:);
for i=2:size(BORDERS,1)
    for j=2:size(Unsorted,1)
        if BORDERS(i-1,2)==Unsorted(j,1)
            BORDERS(i,:)=Unsorted(j,:);
            break
        end
    end
end
for i=1:size(BORDERS,1)
    Nbound=i;
    Nnode=BORDERS(i,1);
    Vbound=INITIAL_POSITIONS(Nbound,:);
    Vnode=INITIAL_POSITIONS(Nnode,:);
    INITIAL_POSITIONS(Nnode,:)=Vbound;
    INITIAL_POSITIONS(Nbound,:)=Vnode;
    for j=i+1:size(BORDERS,1)
        if BORDERS(j,1)==Nbound
            BORDERS(j,1)=Nnode;
        end
    end
    for j=1:size(TRIANGULATION,1)
        if TRIANGULATION(j,1)==Nnode
            TRIANGULATION(j,1)=Nbound;
        elseif TRIANGULATION(j,1)==Nbound
            TRIANGULATION(j,1)=Nnode;
        end
        if TRIANGULATION(j,2)==Nnode
            TRIANGULATION(j,2)=Nbound;
        elseif TRIANGULATION(j,2)==Nbound
            TRIANGULATION(j,2)=Nnode;
        end
        if TRIANGULATION(j,3)==Nnode
            TRIANGULATION(j,3)=Nbound;
        elseif TRIANGULATION(j,3)==Nbound
            TRIANGULATION(j,3)=Nnode;
        end
    end
end
BORDERS(:,1)=[1:size(BORDERS,1)]';
BORDERS(1:size(BORDERS,1)-1,2)=[2:size(BORDERS,1)]';
BORDERS(size(BORDERS,1),2)=1;
BORDERS={'Closed',size(BORDERS,1),BORDERS(:,1),Interpolant};