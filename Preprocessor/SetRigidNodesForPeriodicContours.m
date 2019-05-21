function [INITIAL_POSITIONS,BORDERS,TRIANGULATION]=SetRigidNodesForPeriodicContours(Contour1,Contour2,Interpolant1,Interpolant2,Average_Nodal_Distance,Mesh_Ratio,Periodic_Boundaries,Activate_Plot)

% Compute regular boundary nodes for contour1
HR=(Contour2(1,2)-Contour1(size(Contour1,1),2))/(max(Contour2(:,2))-min(Contour1(:,2)));
Ltot=0;
Lengths=zeros(size(Contour1,1),1);
for j=1:size(Contour1,1)-1
    Ltot=Ltot+sqrt((Contour1(j+1,1)-Contour1(j,1))^2+(Contour1(j+1,2)-Contour1(j,2))^2);
    Lengths(j+1,1)=Ltot;
end
Ninterp=round(Ltot/(HR*Average_Nodal_Distance*Mesh_Ratio(1)));
DistInterp=Ltot/Ninterp;
InterpContour1=zeros(Ninterp+1,2);
InterpContour1(1,:)=Contour1(1,1:2);
InterpContour1(size(InterpContour1,1),:)=Contour1(size(Contour1,1),1:2);
for i=2:Ninterp
    D=(i-1)*DistInterp;
    for j=1:size(Lengths,1)
        if D>=Lengths(j,1) & D<=Lengths(j+1,1)
            ratio=(D-Lengths(j,1))/(Lengths(j+1,1)-Lengths(j,1));
            InterpContour1(i,:)=(1-ratio)*Contour1(j,1:2)+ratio*Contour1(j+1,1:2);
            break
        end
    end
end

% Compute regular boundary nodes for contour2
Ltot=0;
Lengths=zeros(size(Contour2,1),1);
for j=1:size(Contour2,1)-1
    Ltot=Ltot+sqrt((Contour2(j+1,1)-Contour2(j,1))^2+(Contour2(j+1,2)-Contour2(j,2))^2);
    Lengths(j+1,1)=Ltot;
end
Ninterp=round(Ltot/(HR*Average_Nodal_Distance*Mesh_Ratio(2)));
DistInterp=Ltot/Ninterp;
InterpContour2=zeros(Ninterp+1,2);
InterpContour2(1,:)=Contour2(1,1:2);
InterpContour2(size(InterpContour2,1),:)=Contour2(size(Contour2,1),1:2);
for i=2:Ninterp
    D=(i-1)*DistInterp;
    for j=1:size(Lengths,1)
        if D>=Lengths(j,1) & D<=Lengths(j+1,1)
            ratio=(D-Lengths(j,1))/(Lengths(j+1,1)-Lengths(j,1));
            InterpContour2(i,:)=(1-ratio)*Contour2(j,1:2)+ratio*Contour2(j+1,1:2);
            break
        end
    end
end

% Close contour
L12=Contour2(1,2)-Contour1(size(Contour1,1),2);

fmin=Average_Nodal_Distance*Mesh_Ratio(2)*0.6*HR;
fmax=Average_Nodal_Distance*Mesh_Ratio(1)*0.6*HR;
Ninterp=floor(2*L12/(fmin+fmax));
fmin=2*L12/Ninterp/(Mesh_Ratio(1)/Mesh_Ratio(2)+1);
fmax=2*L12/Ninterp*(1-1/(Mesh_Ratio(1)/Mesh_Ratio(2)+1));
Dd=(fmin-fmax)/(Ninterp-1);
Contour12=zeros(Ninterp-1,2);
for i=1:size(Contour12,1)
    D=i*fmax+(i*(i-1)/2)*Dd;
    Contour12(i,1)=Periodic_Boundaries(2);
    Contour12(i,2)=Contour1(size(Contour1,1),2)+D;
end
Contour21=flipud(Contour12);
Contour21(:,1)=ones(size(Contour21,1),1)*Periodic_Boundaries(1);
Contour=cat(1,InterpContour1,Contour12,InterpContour2,Contour21);
Contour=[Contour;Contour(1,:)];

%%%%
%figure;plot(InterpContour1(:,1),InterpContour1(:,2),'.b');hold on
%plot(InterpContour2(:,1),InterpContour2(:,2),'.r')
%plot(Contour12(:,1),Contour12(:,2),'.g')
%plot(Contour21(:,1),Contour21(:,2),'.k')
%axis equal;pause
%%%%

% Position the field nodes and triangulate
xdist=max(Contour(:,1))-min(Contour(:,1));
ydist=max(Contour(:,2))-min(Contour(:,2));
Box=[min(Contour(:,1))-xdist/100,min(Contour(:,2))-ydist/100;max(Contour(:,1))+xdist/100,max(Contour(:,2))+ydist/100];
pmin=min(Contour(:,2));
pmax=max(Contour(:,2));
fmin=Average_Nodal_Distance*Mesh_Ratio(1);
fmax=Average_Nodal_Distance*Mesh_Ratio(2);
a=(fmax-fmin)/(pmax-pmin);
b=fmin-pmin*a;
fh=@(p) b+a*p(:,2);
[INITIAL_POSITIONS,TRIANGULATION]=distmesh2d_plot(Activate_Plot,0.01,@dpoly,fh,4*0.99*Average_Nodal_Distance,Box,Contour,Contour);

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

% Introduce periodicity in nodes and triangulation
ListLeft=[];
ListRight=[];
for i=1:size(INITIAL_POSITIONS,1)
    if abs(INITIAL_POSITIONS(i,1)-Periodic_Boundaries(1))<1e-2*Average_Nodal_Distance
        ListLeft=cat(1,ListLeft,[i,INITIAL_POSITIONS(i,2)]);
    elseif abs(INITIAL_POSITIONS(i,1)-Periodic_Boundaries(2))<1e-2*Average_Nodal_Distance
        ListRight=cat(1,ListRight,[i,INITIAL_POSITIONS(i,2)]);
    end
end
ListLeft=sortrows(ListLeft,2);
ListLeft=ListLeft(:,1);
UpperLeft=ListLeft(size(ListLeft,1));
LowerLeft=ListLeft(1);
ListRight=sortrows(ListRight,2);
ListRight=ListRight(:,1);
for i=1:size(TRIANGULATION,1)
    if isempty(find(ListRight==TRIANGULATION(i,1)))==0
        TRIANGULATION(i,1)=ListLeft(find(ListRight==TRIANGULATION(i,1)));
    end
    if isempty(find(ListRight==TRIANGULATION(i,2)))==0
        TRIANGULATION(i,2)=ListLeft(find(ListRight==TRIANGULATION(i,2)));
    end
    if isempty(find(ListRight==TRIANGULATION(i,3)))==0
        TRIANGULATION(i,3)=ListLeft(find(ListRight==TRIANGULATION(i,3)));
    end
end
ListKept=[];
for i=1:size(INITIAL_POSITIONS,1)
    if isempty(find(ListRight==i))
        ListKept=cat(1,ListKept,i);
    end
end
INITIAL_POSITIONS=INITIAL_POSITIONS(ListKept,:);
for i=1:size(TRIANGULATION,1)
    TRIANGULATION(i,1)=find(ListKept==TRIANGULATION(i,1));
    TRIANGULATION(i,2)=find(ListKept==TRIANGULATION(i,2));
    TRIANGULATION(i,3)=find(ListKept==TRIANGULATION(i,3));
end
UpperLeft=find(ListKept==UpperLeft);
LowerLeft=find(ListKept==LowerLeft);

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
Border1=LowerLeft;
Done=zeros(size(Unsorted,1),1);
prev=LowerLeft;
while Border1(1)~=Border1(size(Border1,1)) | size(Border1,1)==1
    for i=1:size(Unsorted,1)
        if Unsorted(i,1)==prev
            if Done(i,1)==0
                Done(i,1)=1;
                prev=Unsorted(i,2);
                Border1=cat(1,Border1,prev);
                break
            end
        elseif Unsorted(i,2)==prev
            if Done(i,1)==0
                Done(i,1)=1;
                prev=Unsorted(i,1);
                Border1=cat(1,Border1,prev);
                break
            end
        end
    end
end
Border1=Border1(1:size(Border1,1)-1);
if Border1(2)>Border1(size(Border1,1))
    Border1=flipud(Border1);
end
Border2=UpperLeft;
Done=zeros(size(Unsorted,1),1);
prev=UpperLeft;
while Border2(1)~=Border2(size(Border2,1)) | size(Border2,1)==1
    for i=1:size(Unsorted,1)
        if Unsorted(i,1)==prev
            if Done(i,1)==0
                Done(i,1)=1;
                prev=Unsorted(i,2);
                Border2=cat(1,Border2,prev);
                break
            end
        elseif Unsorted(i,2)==prev
            if Done(i,1)==0
                Done(i,1)=1;
                prev=Unsorted(i,1);
                Border2=cat(1,Border2,prev);
                break
            end
        end
    end
end
Border2=Border2(1:size(Border2,1)-1);
if Border2(2)<Border2(size(Border2,1))
    Border2=flipud(Border2);
end
%
Border1=[Border1;Border1(1)];
Border2=[Border2;Border2(1)];
%
BORDERS={'Periodic',size(Border1,1),Border1,Interpolant1,0,0;'Periodic',size(Border2,1),Border2,Interpolant2,1,1};