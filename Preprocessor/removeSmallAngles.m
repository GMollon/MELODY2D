function [ outVertices,outCells ] = removeSmallAngles(vertices,cells,critLength,critAngle,flagPlotDetails,flagPlotSteps,flagPlotFinal,boundaries,maxIterations,removeOnlyAntiClockwiseConvex)%vertices,cells,plotSteps)
%refineVoronoi Modifies a Voronoi-Cell microstructure as to eliminate
%angles of problematic acuteness.
%   Input arguments 'vertices' and 'cells' represent a 2D cellular
%   structure as created e.g by the function voronoin. 'critAngle' [radian] is the target value for minimum inner
%   vertex angle (of any cell). 'critLength' Is the minimum lenght of new Edges introduced to the structure. Plot Options (logicals) allow
%   to illustrate the progress of the refinement. The Refinement is done by
%   an iterative algorithm, set 'maxIterations' to limit calculation time.
%   Please mind that in the case of interruption by 'maxIterations' there
%   are still critical angles present (during runtime, see command window
%   messages for more information about the results)
%
%
%
%   last modified: 15.12.2016
%   by: Simon Massa
axesLim=boundaries;%[xmin xmax ymin ymax];
%
if flagPlotSteps==true && maxIterations>5
    warning(['using plotSteps = true with given number of maximum iterations will result in creation of up to ' num2str(maxIterations+6) ' figures'])
    str = input('continue execution? [y/n]','s');
    if str=='y'
        
    else
        error('terminated by user');
    end
end

%====================================================
%Begin refinement
disp('Refine Voronoi Structure: Remove small angles');
disp(['min angle in cell corner: ' num2str(critAngle*360/(2*pi)) '°']);
disp(['max iterations allowed: ' num2str(maxIterations)]);

%internal parameters:
safety=1.05;
histo=false;
maximise=true;

%initializations
itCount=1;
numAnglesReplaced=1; %initialization to allow entrance to loop
fact=0.9;
warnTriangle=false;



[vertices,cells]=cleanStructure(vertices,cells);

if flagPlotDetails
    figure();
    fig1=gcf;
    title('original microstructure');
    axis equal; hold on;
    axis(axesLim);
    if maximise set(gcf, 'Position', fact*get(0,'Screensize')); end % Maximize figure.
    
    figure();
    fig5=gcf;
    title('Structure with critical angles highlited');
    axis equal; hold on;
    axis(axesLim);
    if maximise set(gcf, 'Position', fact*get(0,'Screensize')); end % Maximize figure.
    
end


mySize=size(cells);
numCells=mySize(1);


if flagPlotDetails
    for i=1:numCells
        %plot cells:
        figure(fig1);
        plotCell(cells(i),vertices,'-b');
        a=vertices(cells{i},:);
        plot(a(:,1),a(:,2),'*b');
        %plotCell(cells(i),vertices,'-b');
        
        
    end
    
end

%%

%Angles
%====================================
%
%find critical angles
angleItCount=0;

while (angleItCount<=maxIterations && numAnglesReplaced>0)
    
    cellsWithCritAng=[];
    numAnglesReplaced=0;
    newVertices=vertices;
    newCells=cells;
    for i=1:numCells
        tempArr=cells{i};
        mySize=size(tempArr);
        Nvert=mySize(2);
        jBis=1;
        
        %only if at least one vertice is in boundarys
        cellIsInArea=false;
        for j=1:Nvert
            if isInArea(vertices(tempArr(j),:),boundaries);
                cellIsInArea=true;
                break;
            end
        end
        
        if cellIsInArea
            if flagPlotDetails
                figure(fig5);
                plotCell(newCells(i),newVertices,'-b');
            end
            
            
            
            for j=1:Nvert
                %mesure angle at each vertice
                
                if j==1
                    vert1=vertices(tempArr(Nvert),:);
                    vert2=vertices(tempArr(1),:);
                    vert3=vertices(tempArr(2),:);
                    
                elseif j==Nvert
                    vert1=vertices(tempArr(Nvert-1),:);
                    vert2=vertices(tempArr(Nvert),:);
                    vert3=vertices(tempArr(1),:);
                else
                    vert1=vertices(tempArr(j-1),:);
                    vert2=vertices(tempArr(j),:);
                    vert3=vertices(tempArr(j+1),:);
                end
                a=vert1-vert2;
                b=vert3-vert2;
                anglea=atan2(a(2),a(1));
                if anglea<0
                    anglea=2*pi+anglea;
                end
                angleb=atan2(b(2),b(1));
                if angleb<0
                    angleb=2*pi+angleb;
                end
                
                angle=(anglea-angleb);
                
                if angle<0
                   angle=2*pi+angle; 
                end
                
                %disp(num2str(angle*360/(2*pi)));
                
                leftturn=angle>0 && angle<pi;
                %rightturn=angle>pi && angle<2*pi;
                
                
                if removeOnlyAntiClockwiseConvex==true
                    isCritical=angle <= critAngle && leftturn;
                else
                    isCritical=angle <= critAngle || angle >= 2*pi-critAngle;
                end
                if isCritical
                    
                    cellsWithCritAng=[cellsWithCritAng ['(' num2str(i) ' ' num2str(j) ') ']];
                    %=========================================
                    %actions to perform if angle is < crit angle:
                    numAnglesReplaced=numAnglesReplaced+1;
                    
                    d=(critLength/2)*safety / (sin(angle/2));
                    relVec=vert1-vert2;
                    relVecNewVert=d*(relVec/norm(relVec));
                    if norm(relVec)-d<critLength
                        %relVecNewVert=relVec./2;
                        P1=NaN;%vert1;
                    else
                        P1=vert2+relVecNewVert;
                    end
                    
                    relVec=vert3-vert2;
                    relVecNewVert=d*(relVec/norm(relVec));
                    if norm(relVec)-d<critLength
                        %relVecNewVert=relVec./2;
                        P2=NaN;%vert3;
                    else
                        P2=vert2+relVecNewVert;
                    end
                    
                    currCellVert=newCells{i};
                    
                    if i==100
                        disp(' ')
                    end
                    if ~(sum(isnan(P1))>0) && ~(sum(isnan(P2))>0)
                        newVertices=[newVertices;P1;P2];
                        mySize=size(newVertices);
                        numVert=mySize(1);
                        currCellVert=[currCellVert(1:jBis-1) numVert-1 numVert currCellVert(jBis+1:end)];
                        plotTri=[P1;vert2;P2];
                        jBis=jBis+1;
                        
                    elseif (sum(isnan(P1))>0) && ~(sum(isnan(P2))>0)
                        newVertices=[newVertices;P2];
                        mySize=size(newVertices);
                        numVert=mySize(1);
                        currCellVert=[currCellVert(1:jBis-1) numVert currCellVert(jBis+1:end)];
                        plotTri=[vert1;vert2;P2];
                        
                        
                    elseif ~(sum(isnan(P1))>0) && (sum(isnan(P2))>0)
                        newVertices=[newVertices;P1];
                        mySize=size(newVertices);
                        numVert=mySize(1);
                        currCellVert=[currCellVert(1:jBis-1) numVert currCellVert(jBis+1:end)];
                        plotTri=[P1;vert2;vert3];
                    elseif (sum(isnan(P1))>0) && (sum(isnan(P2))>0)
                        if size(currCellVert,2)>3 %check if cell is not a triangle
                            currCellVert=[currCellVert(1:jBis-1) currCellVert(jBis+1:end)];
                            plotTri=[vert1;vert2;vert3];
                            jBis=jBis-1;
                        else
                            warnTriangle=true;
                            numAnglesReplaced=numAnglesReplaced-1;
                        end
                    end
                    
                    %                 if size(currCellVert,2)<3
                    %                     disp('I did bullshit')
                    %                 end
                    
                    newCells{i}=currCellVert;
                    %=========================================
                    if flagPlotDetails
                        figure(fig5);
                        fill(plotTri(:,1),plotTri(:,2),'r');
                    end
                    break;% !!!!!
                end
                jBis=jBis+1;
            end
        end
        
        
        
        
        
        
    end
    
    cells=newCells;
    vertices=newVertices;
    
    
    [vertices,cells]=cleanStructure(vertices,cells);
    
    
    disp(['End of refinement iteration ' num2str(itCount)]);
    disp(['Number of critical angles replaced by edge: ' num2str(numAnglesReplaced)]);
    disp(['Was working on (cells vertice): ' num2str(cellsWithCritAng)]);
    itCount=itCount+1;
    
    angleItCount=angleItCount+1;
    
    
    if flagPlotSteps
        figure();
        fig6=gcf;
        title(['Structure after Iteration ' num2str(itCount)]);
        axis equal; hold on;
        axis(axesLim);
        if maximise set(gcf, 'Position', fact*get(0,'Screensize')); end % Maximize figure.
        for i=1:numCells
            
            figure(fig6);
            plotCell(cells(i),vertices,'-b');
        end
    end
end




if numAnglesReplaced==0
    disp('Algorithm converged successfully. Small angles removed');
else
    disp(['Output -structure contains critical angles and/or small edges, consider allowing more Iterations. Stopped after ' num2str(itCount-1) ' iteration(s).']);
end

if warnTriangle
    warning('triangular cells could not be treated, might contain critical angles/edges, consider allowing larger critical Length')
end



outVertices=vertices;
outCells=cells;

if flagPlotFinal
    figure();
    title('Output Structure (small angles removed)');
    axis equal; hold on;
    axis(axesLim);
    numCells=size(cells,1);
    for i=1:numCells
        plotCell(outCells(i),outVertices,'-b');
    end
    if maximise set(gcf, 'Position', fact*get(0,'Screensize')); end % Maximize figure.
    
end

end


function plotCell(cell,vertices,linestyle);
plot(vertices([cell{1,1}';cell{1,1}(1)'],1),vertices([cell{1,1}';cell{1,1}(1)'],2),linestyle);
end

function [ diffMat ] = differenceMatrix( mat )

diffMat=mat(2:end,:)-mat(1:end-1,:);

end
function [flag]=isInArea(point,bounds);
flag=false;
x=point(1,1);
y=point(1,2);
xmin=bounds(1,1);
xmax=bounds(1,2);
ymin=bounds(1,3);
ymax=bounds(1,4);

if x<=xmax && x>=xmin && y<=ymax && y>=ymin
    flag=true;
end
end
function [vertices,cells]=cleanStructure(vertices,cells)
%============================
%clean up unused vertices
%explanation: during the process of cutting out small edges, some vertices
%disappear as they are no longer used by any cell. These vertices are
%cleared out from the vertices matrix, also the indexes in the cells are
%transformed to fit the new length of the vertices matrix
mySize=size(cells);
numCells=mySize(1);
vertInd=[];
for i=1:numCells
    vertInd=[vertInd; cells{i}'];
end
vertInd=unique(vertInd);
mySize=size(vertInd);
numVertClean=mySize(1);
inds=1:numVertClean;
inds=inds';
vertices=vertices(vertInd,:);
indexTransformation=[vertInd inds];


for i=1:numCells
    tempArr=cells{i};
    mySize=size(tempArr);
    Nvert=mySize(2);
    for j=1:Nvert
        tempArr(j)=indexTransformation(indexTransformation(:,1)==tempArr(j),2);
    end
    cells{i}=tempArr;
    
    
end
% end clean up unused vertices
end