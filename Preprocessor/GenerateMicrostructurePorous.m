function [C,Cup,Cdown]=GenerateMicrostructurePorous(Type,LowerBoundary,UpperBoundary,SideBoundaries,Npoints,xmin,xmax,ymin,ymax,noise,critLength,critAngle,voidFrac,cumulativeFreq)

boundaries= [xmin xmax ymin ymax];
%Data%
%Type='Halton';
%LowerBoundary='Embedded';
%UpperBoundary='Flat';
%Npoints=100;
%xmin=0;
%xmax=5;
%ymin=0;
%ymax=1;


%Boundaries%
if strcmp(Type,'Halton')
    Points=Halton(Npoints,2);
    Points(:,1)=xmin+(xmax-xmin)*Points(:,1);
    Points(:,2)=ymin+(ymax-ymin)*Points(:,2);
elseif strcmp(Type,'Rand')
    Points=rand(Npoints,2);
    Points(:,1)=xmin+(xmax-xmin)*Points(:,1);
    Points(:,2)=ymin+(ymax-ymin)*Points(:,2);
elseif strcmp(Type,'Lattice')
    Points=Lattice(Npoints,xmin,xmax,ymin,ymax,noise);
    Npoints=size(Points,1);
end
if strcmp(SideBoundaries,'Periodic')
    Pointsxp=[Points(:,1)+(xmax-xmin),Points(:,2)];
    Pointsxm=[Points(:,1)-(xmax-xmin),Points(:,2)];
elseif strcmp(SideBoundaries,'Flat')
    Pointsxp=[xmax-(Points(:,1)-xmax),Points(:,2)];
    Pointsxm=[xmin-(Points(:,1)-xmin),Points(:,2)];
end
Points=[Points;Pointsxp;Pointsxm];
if strcmp(LowerBoundary,'Flat')
    Pointsym=[Points(:,1),ymin-(Points(:,2)-ymin)];
elseif strcmp(LowerBoundary,'Embedded')
    Pointsym=[Points(:,1),Points(:,2)-(ymax-ymin)];
end
if strcmp(UpperBoundary,'Flat')
    Pointsyp=[Points(:,1),ymax+(ymax-Points(:,2))];
elseif strcmp(UpperBoundary,'Embedded')
    Pointsyp=[Points(:,1),Points(:,2)+(ymax-ymin)];
end
Points=[Points;Pointsyp;Pointsym];

%plot(Points(:,1),Points(:,2),'r*'); axis equal; hold on;
%axis([0 1 0 1]);
%Tesselation%
[vertices,cells]=voronoin(Points); %fct matlab

%________________________________________
%plot structure before refinement

if false
    figure();
    hold on; axis equal;
    axis(boundaries);
    title('Voronoi Cells')
    
    for i=1:size(cells,1)
        %plot cells:
        plotCell(cells(i),vertices,'-b');
    end
    
end



%________________________________________
%refinement of structure
%critLength=0.05;
%critAngle=pi/2*1.01;
plotDetails=false;
plotSteps=false;
plotFinal=false;
maxIterations=500;
[vertices,cells]=removeSmallEdges(vertices,cells,critLength,plotDetails,plotFinal,boundaries,maxIterations);
%Simon Massa 16.11.2016
%________________________________________




%Detecting Boundary Segments%
Segments=[];
SegmentsUp=[];
SegmentsDown=[];
for i=1:3*Npoints
    for j=1:size(cells{i,1},2)-1
        Segments=cat(1,Segments,[cells{i,1}(1,j),cells{i,1}(1,j+1),i]);
    end
    Segments=cat(1,Segments,[cells{i,1}(1,j+1),cells{i,1}(1,1),i]);
end
Segments(:,1:2)=sort(Segments(:,1:2)')';
Segments=sortrows(Segments,[1,2,3]);
for i=2:size(Segments,1)-1
    if Segments(i,1)~=Segments(i-1,1) | Segments(i,2)~=Segments(i-1,2)
        if Segments(i,1)~=Segments(i+1,1) | Segments(i,2)~=Segments(i+1,2)
            if vertices(Segments(i,1),2)<(ymin+ymax)/2 & Segments(i,3)<=Npoints
                SegmentsDown=cat(1,SegmentsDown,Segments(i,:));
            elseif Segments(i,3)<=Npoints
                SegmentsUp=cat(1,SegmentsUp,Segments(i,:));
            end
        end
    end
end
i=1;
if Segments(i,1)~=Segments(i+1,1) & Segments(i,2)~=Segments(i+1,2)
    if vertices(Segments(i,1),2)<(ymin+ymax)/2 & Segments(i,3)<=Npoints
        SegmentsDown=cat(1,SegmentsDown,Segments(i,:));
    elseif Segments(i,3)<=Npoints
        SegmentsUp=cat(1,SegmentsUp,Segments(i,:));
    end
end
i=size(Segments,1);
if Segments(i,1)~=Segments(i-1,1) & Segments(i,2)~=Segments(i-1,2)
    if vertices(Segments(i,1),2)<(ymin+ymax)/2 & Segments(i,3)<=Npoints
        SegmentsDown=cat(1,SegmentsDown,Segments(i,:));
    elseif Segments(i,3)<=Npoints
        SegmentsUp=cat(1,SegmentsUp,Segments(i,:));
    end
end

%Sorting Segments%
VerticesDown=sort(reshape(SegmentsDown(:,1:2),[],1));
VerticesDown(:,2)=vertices(VerticesDown,1);
VerticesDown=sortrows(VerticesDown,2);
done=0;
last=VerticesDown(size(VerticesDown,1),1);
ContourDown=last;
while done<size(SegmentsDown,1)
    for i=1:size(SegmentsDown,1)
        if SegmentsDown(i,1)==last & SegmentsDown(i,3)>0
            last=SegmentsDown(i,2);
            SegmentsDown(i,3)=0;
            ContourDown=[ContourDown;last];
            done=done+1;
            break
        elseif SegmentsDown(i,2)==last & SegmentsDown(i,3)>0
            last=SegmentsDown(i,1);
            SegmentsDown(i,3)=0;
            ContourDown=[ContourDown;last];
            done=done+1;
            break
        end
    end
end
cdown=vertices(ContourDown,:);
if cdown(1,1)>=xmax
    for i=2:size(cdown,1)
        if cdown(i,1)<xmax
            newpoint=[xmax,cdown(i-1,2)+(xmax-cdown(i-1,1))*(cdown(i,2)-cdown(i-1,2))/(cdown(i,1)-cdown(i-1,1))];
            Cdown=[newpoint;cdown(i:size(cdown,1)-1,:);[cdown(1:i-1,1)-(xmax-xmin),cdown(1:i-1,2)];[xmin,newpoint(2)]];
            break
        end
    end
elseif cdown(size(cdown,1),1)<=xmin
    cdown=flipud([(xmin-cdown(:,1)+xmax),cdown(:,2)]);
    for i=2:size(cdown,1)
        if cdown(i,1)<xmax
            newpoint=[xmax,cdown(i-1,2)+(xmax-cdown(i-1,1))*(cdown(i,2)-cdown(i-1,2))/(cdown(i,1)-cdown(i-1,1))];
            Cdown=[newpoint;cdown(i:size(cdown,1)-1,:);[cdown(1:i-1,1)-(xmax-xmin),cdown(1:i-1,2)];[xmin,newpoint(2)]];
            break
        end
    end
    Cdown=flipud([(xmin-Cdown(:,1)+xmax),Cdown(:,2)]);
end

VerticesUp=sort(reshape(SegmentsUp(:,1:2),[],1));
VerticesUp(:,2)=vertices(VerticesUp,1);
VerticesUp=sortrows(VerticesUp,2);
done=0;
last=VerticesUp(1,1);
ContourUp=last;
while done<size(SegmentsUp,1)
    for i=1:size(SegmentsUp,1)
        if SegmentsUp(i,1)==last & SegmentsUp(i,3)>0
            last=SegmentsUp(i,2);
            SegmentsUp(i,3)=0;
            ContourUp=[ContourUp;last];
            done=done+1;
            break
        elseif SegmentsUp(i,2)==last & SegmentsUp(i,3)>0
            last=SegmentsUp(i,1);
            SegmentsUp(i,3)=0;
            ContourUp=[ContourUp;last];
            done=done+1;
            break
        end
    end
end
cup=flipud(vertices(ContourUp,:));
if cup(1,1)>=xmax
    for i=2:size(cup,1)
        if cup(i,1)<xmax
            newpoint=[xmax,cup(i-1,2)+(xmax-cup(i-1,1))*(cup(i,2)-cup(i-1,2))/(cup(i,1)-cup(i-1,1))];
            Cup=[newpoint;cup(i:size(cup,1)-1,:);[cup(1:i-1,1)-(xmax-xmin),cup(1:i-1,2)];[xmin,newpoint(2)]];
            break
        end
    end
elseif cup(size(cup,1),1)<=xmin
    cup=flipud([(xmin-cup(:,1)+xmax),cup(:,2)]);
    for i=2:size(cup,1)
        if cup(i,1)<xmax
            newpoint=[xmax,cup(i-1,2)+(xmax-cup(i-1,1))*(cup(i,2)-cup(i-1,2))/(cup(i,1)-cup(i-1,1))];
            Cup=[newpoint;cup(i:size(cup,1)-1,:);[cup(1:i-1,1)-(xmax-xmin),cup(1:i-1,2)];[xmin,newpoint(2)]];
            break
        end
    end
    Cup=flipud([(xmin-Cup(:,1)+xmax),Cup(:,2)]);
end
Cup=flipud(Cup);

%________________________________________
%refinement of structure
[vertices,cells]=removeSmallAngles(vertices,cells,critLength,critAngle*2*pi/360,plotDetails,plotSteps,plotFinal,boundaries,maxIterations,false);
%Simon Massa 16.11.2016
%________________________________________


%________________________________________
%treat Cdown
%figure();
%plot(Cdown(:,1),Cdown(:,2))
%axis equal
plotDetails=false;plotSteps=false;plotFinal=false;
width=max(Cdown(:,1))-min(Cdown(:,1));
CdownPre=[Cdown(:,1)-width Cdown(:,2)];
CdownPost=[Cdown(:,1)+width Cdown(:,2)];
CDownVertices=[CdownPost;Cdown;CdownPre];
CDownCells{1}=1:size(CDownVertices,1);
disp('treat Cdown');
[CDownVertices,CDownCells]=removeSmallAngles(CDownVertices,CDownCells,critLength,critAngle*2*pi/360,plotDetails,plotSteps,plotFinal,[-10 10 -10 10],maxIterations,true);
Cdown=CDownVertices(CDownCells{1},:);



%find vertices closest to x=zero and x=width
%Cdown
numVert=size(Cdown,1);
dist0=width;
dist1=width;
ind0=-1;
ind1=-1;
for i=1:numVert
    currDist0=abs(Cdown(i,1)-0);
    currDist1=abs(Cdown(i,1)-width);
    if currDist0<dist0
        dist0=currDist0;
        ind0=i;
    end
    if currDist1<dist1
        dist1=currDist1;
        ind1=i;
    end
end
Cdown=Cdown(ind1:ind0,:);

Contours=cell(Npoints,1);
for i=1:Npoints
    Contours{i,1}=[vertices(cells{i,1}',:);vertices(cells{i,1}(1),:)];
end
if false
    figure
    hold on
    %title('before porosity');
    
    for i=1:Npoints
        %fill(Contours{i,1}(:,1),Contours{i,1}(:,2),'k')
        plot(Contours{i,1}(:,1),Contours{i,1}(:,2),'-b')
        hold on
    end
end
%introduce Porosity
%simple method
%calculate area of Polygons:
surfContours=surfaceContours(Contours);
[surfContoursAsc,inds]=sort(surfContours,'ascend');
%deleteCells untill reaching aimed void space
totalSurface=sum(surfContours);
aimVoidSpace=voidFrac*totalSurface;
deletedSurface=0;
i=1;
indsDel=[];
while deletedSurface<aimVoidSpace
    indsDel=[indsDel inds(i)];
    deletedSurface=deletedSurface+surfContoursAsc(i);
    i=i+1;
    
end
Contours(indsDel)=[];

disp(['num Contours deleted for porosity: ' num2str(i-1)]);
voidFracOut=deletedSurface/totalSurface;
Npoints=Npoints-i+1;

%Simon Massa 26.11.2016
%________________________________________




C=Contours;
%Plotting%
if false
    figure
    hold on
    title('Output Structure');
    axis equal;
    for i=1:Npoints
        if i==1
            disp('');
        end
        fill(C{i,1}(:,1),C{i,1}(:,2),'k')
        plot(C{i,1}(:,1),C{i,1}(:,2),'-b')
        hold on
    end
    axis equal
    plot(Cdown(:,1),Cdown(:,2),'-r')
    plot(Cup(:,1),Cup(:,2),'-g')
end
end


function surfContours=surfaceContours(cont)
numcont=size(cont,1);
A=zeros(numcont,1);
for i=1:numcont
    
    tempVert=cont{i};
    A(i)=polyarea(tempVert(:,1),tempVert(:,2));
    
end
surfContours=A;
end


function plotCell(cell,vertices,linestyle);
plot(vertices([cell{1,1}';cell{1,1}(1)'],1),vertices([cell{1,1}';cell{1,1}(1)'],2),linestyle);
end

