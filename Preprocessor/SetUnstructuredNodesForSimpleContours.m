function [INITIAL_POSITIONS,BORDERS,TRIANGULATION]=SetUnstructuredNodesForSimpleContours(Contours,Average_Nodal_Distance,Activate_Plot)
for i=1:size(Contours,1)
    if strcmp(Contours{i,1},'Simple')==0
        i=i-1;
        break
    end
end
Ncontours=i;

% Compute regular boundary nodes for all contours
Contour=[];
for n=1:Ncontours
    Ltot=0;
    Lengths=zeros(size(Contours{n,2},1),1);
    for j=1:size(Contours{n,2},1)-1
        Ltot=Ltot+sqrt((Contours{n,2}(j+1,1)-Contours{n,2}(j,1))^2+(Contours{n,2}(j+1,2)-Contours{n,2}(j,2))^2);
        Lengths(j+1,1)=Ltot;
    end
    Ninterp=round(Ltot/Average_Nodal_Distance);
    DistInterp=Ltot/Ninterp;
    InterpContour=zeros(Ninterp,2);
    InterpContour(1,:)=Contours{n,2}(1,1:2);
    InterpContour(size(InterpContour,1),:)=Contours{n,2}(size(Contours{n,2},1),1:2);
    for i=2:Ninterp-1
        D=(i-1)*DistInterp;
        for j=1:size(Lengths,1)
            if D>=Lengths(j,1) & D<=Lengths(j+1,1)
                ratio=(D-Lengths(j,1))/(Lengths(j+1,1)-Lengths(j,1));
                InterpContour(i,:)=(1-ratio)*Contours{n,2}(j,1:2)+ratio*Contours{n,2}(j+1,1:2);
                break
            end
        end
    end
    if isempty(Contour)==0
        InterpContour=InterpContour(2:size(InterpContour,1),:);
    end
    Contour=cat(1,Contour,InterpContour);
end

% Position the field nodes and triangulate
Box=[min(Contour(:,1)),min(Contour(:,2));max(Contour(:,1)),max(Contour(:,2))];
[INITIAL_POSITIONS,TRIANGULATION]=distmesh2d_plot(Activate_Plot,0.001,@dpoly,@huniform,Average_Nodal_Distance,Box,Contour,Contour);

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

% Detect the boundary nodes
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

% Sort the boundary nodes
BORDERS=cell(Ncontours,3);
for i=1:size(Unsorted,1)
    if norm([INITIAL_POSITIONS(Unsorted(i,1),1)-Contours{1,2}(1,1),INITIAL_POSITIONS(Unsorted(i,1),2)-Contours{1,2}(1,2)])<1e-6*Average_Nodal_Distance
        prev=Unsorted(i,1);
        break
    elseif norm([INITIAL_POSITIONS(Unsorted(i,2),1)-Contours{1,2}(1,1),INITIAL_POSITIONS(Unsorted(i,2),2)-Contours{1,2}(1,2)])<1e-6*Average_Nodal_Distance
        prev=Unsorted(i,2);
        break
    end
end
Done=zeros(size(Unsorted,1),1);
Border=prev;
while Border(1)~=Border(size(Border,1)) | size(Border,1)==1
    for i=1:size(Unsorted,1)
        if Unsorted(i,1)==prev
            if Done(i,1)==0
                Done(i,1)=1;
                prev=Unsorted(i,2);
                Border=cat(1,Border,prev);
                break
            end
        elseif Unsorted(i,2)==prev
            if Done(i,1)==0
                Done(i,1)=1;
                prev=Unsorted(i,1);
                Border=cat(1,Border,prev);
                break
            end
        end
    end
end
v1=[INITIAL_POSITIONS(Border(2),1)-INITIAL_POSITIONS(Border(1),1),INITIAL_POSITIONS(Border(2),2)-INITIAL_POSITIONS(Border(1),2)];
v2=[Contours{1,2}(2,1)-Contours{1,2}(1,1),Contours{1,2}(2,2)-Contours{1,2}(1,2)];
if abs(dot(v1/norm(v1),v2/norm(v2))-1)>1e-8
    Border=flipud(Border);
end
deb=1;
for n=1:Ncontours
    for i=1:size(Unsorted,1)
        if norm([INITIAL_POSITIONS(Unsorted(i,1),1)-Contours{n,2}(size(Contours{n,2},1),1),INITIAL_POSITIONS(Unsorted(i,1),2)-Contours{n,2}(size(Contours{n,2},1),2)])<1e-6*Average_Nodal_Distance
            last=Unsorted(i,1);
            break
        elseif norm([INITIAL_POSITIONS(Unsorted(i,2),1)-Contours{n,2}(size(Contours{n,2},1),1),INITIAL_POSITIONS(Unsorted(i,2),2)-Contours{n,2}(size(Contours{n,2},1),2)])<1e-6*Average_Nodal_Distance
            last=Unsorted(i,2);
            break
        end
    end
    for i=deb:size(Border,1)
        if Border(i,1)==last
            BORDERS{n,1}='Simple';
            BORDERS{n,2}=i-deb+1;
            BORDERS{n,3}=Border(deb:i,1);
            BORDERS{n,4}=Contours{n,3};
            BORDERS{n,5}=Contours{n,4}-1;
            BORDERS{n,6}=Contours{n,5}-1;
            deb=i;
            break
        end
    end
end