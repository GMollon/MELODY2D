function [INITIAL_POSITIONS,BORDERS,TRIANGULATION]=SetRigidNodesForPolygonContour(Contour,Interpolant,Average_Nodal_Distance,Activate_Plot)

% Close the contour if needed
if norm(Contour(1,:)-Contour(size(Contour,1),:))>1e-3*Average_Nodal_Distance
    Contour=[Contour;Contour(1,:)];
else
    Contour(size(Contour,1),:)=Contour(1,:);
end
Contours=cell(size(Contour,1)-1,1);
Ncontours=size(Contours,1);
for i=1:Ncontours
    Contours{i,1}=[Contour(i,1),Contour(i,2);Contour(i+1,1),Contour(i+1,2)];
end

% Compute regular boundary nodes
%Contour=[];
%for n=1:Ncontours
%    Ltot=sqrt((Contours{n,1}(2,1)-Contours{n,1}(1,1))^2+(Contours{n,1}(2,2)-Contours{n,1}(1,2))^2);
%    Ninterp=ceil(Ltot/Average_Nodal_Distance);
%    DistInterp=Ltot/Ninterp;
%    InterpContour=zeros(Ninterp+1,2);
%    InterpContour(1,:)=Contours{n,1}(1,1:2);
%    InterpContour(size(InterpContour,1),:)=Contours{n,1}(size(Contours{n,1},1),1:2);
%    for i=2:Ninterp
%        D=(i-1)*DistInterp;
%        ratio=D/Ltot;
%        InterpContour(i,:)=(1-ratio)*Contours{n,1}(1,1:2)+ratio*Contours{n,1}(2,1:2);
%    end
%    if isempty(Contour)==0
%        InterpContour=InterpContour(2:size(InterpContour,1),:);
%    end
%    Contour=cat(1,Contour,InterpContour);
%end

Contour=[];
for n=1:Ncontours
    Ltot=sqrt((Contours{n,1}(2,1)-Contours{n,1}(1,1))^2+(Contours{n,1}(2,2)-Contours{n,1}(1,2))^2);
    Ninterp=ceil(Ltot/Average_Nodal_Distance);
    DistInterp=Ltot/Ninterp;
    InterpContour=zeros(Ninterp+1,2);
    InterpContour(1,:)=Contours{n,1}(1,1:2);
    InterpContour(end,:)=Contours{n,1}(end,1:2);
    if Ninterp>1
        for i=2:Ninterp
            D=(i-1)*DistInterp;
            ratio=D/Ltot;
            InterpContour(i,:)=(1-ratio)*Contours{n,1}(1,1:2)+ratio*Contours{n,1}(2,1:2);
        end
    end
    if isempty(Contour)==0
        InterpContour=InterpContour(2:size(InterpContour,1),:);
    end
    Contour=cat(1,Contour,InterpContour);
end

% Position the field nodes and triangulate
if size(Contour,1)>4
    Box=[min(Contour(:,1)),min(Contour(:,2));max(Contour(:,1)),max(Contour(:,2))];
    [INITIAL_POSITIONS,TRIANGULATION]=distmesh2d_plot(Activate_Plot,0.001,@dpoly,@huniform,2*Average_Nodal_Distance,Box,Contour,Contour);
else
    INITIAL_POSITIONS=Contour(1:3,:);
    TRIANGULATION=[1,2,3];
end
%close

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
if size(Toremove,1)>0
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
else
    Unsorted=[];
    for i=1:size(Segments,1)
        if Segments(i,3)==1
            Unsorted(i,1:2)=Segments(i,[1,2]);
        elseif Segments(i,3)==2
            Unsorted(i,1:2)=Segments(i,[2,1]);
        end
    end
end
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
%BORDERS(1:size(BORDERS,1)-1,2)=[2:size(BORDERS,1)]';
%BORDERS(size(BORDERS,1),2)=1;
BORDERS(end+1,1)=BORDERS(1,1);
%
BORDERS={'Closed',size(BORDERS,1),BORDERS(:,1),Interpolant,0,0};
close