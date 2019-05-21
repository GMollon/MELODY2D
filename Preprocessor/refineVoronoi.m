function [ outVertices,outCells ] = refineVoronoi(vertices,cells,critLength,critAngle,flagPlotDetails,flagPlotSteps,flagPlotFinal,boundaries,maxIterations)%vertices,cells,plotSteps)
%refineVoronoi Modifies a Voronoi-Cell microstructure as to eliminate
%angles of problematic acuteness and edges that are to small.
%   Input arguments 'vertices' and 'cells' represent a 2D cellular
%   structure as created e.g by the function voronoin. 'critLength' and
%   'critAngle' [radian] are the target values for minimum edge length and inner
%   vertex angle (of any cell) respectively. Plot Options (logicals) allow
%   to illustrate the progress of the refinement. The Refinement is done by
%   an iterative algorithm, set 'maxIterations' to limit calculation time.
%   Please mind that in the case of interruption by 'maxIterations' there
%   are still critical angles and or edges present (during runtime, see command window
%   messages for more information about the results)
%   refineVoronoi also sorts out cells that contain vertices with 'inf'
%   coordinates.
%
%
%   last modified: 13.11.2016
%   by: Simon Massa
axesLim=boundaries;%[0.44 0.6 0.7 0.86];
%


if flagPlotDetails==true && maxIterations>2
    warning(['using plotSteps = true with given number of maximum iterations will result in creation of up to ' num2str(maxIterations*6) ' figures'])
    str = input('continue execution? [y/n]','s');
    if str=='y'
        
    else
        error('terminated by user');
    end
end

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
disp('Refine Voronoi Structure');
disp(['min Length of Edge: ' num2str(critLength)]);
disp(['min angle in cell corner: ' num2str(critAngle*360/(2*pi)) '°']);
disp(['max iterations allowed: ' num2str(maxIterations)]);

%internal parameters:
safety=1.05;
histo=false;
maximise=true;

%initializations
itCount=1;
numEdgesReplaced=1;
numAnglesReplaced=1; %initialization to allow entrance to loop
fact=0.9;
warnTriangle=false;



[vertices,cells]=cleanStructure(vertices,cells);

if flagPlotDetails
    figure();
    fig1=gcf;
    title('original microstructure, cleaned for inf values');
    axis equal; hold on;
    axis(axesLim);
    if maximise set(gcf, 'Position', fact*get(0,'Screensize')); end % Maximize figure.
    
    
    figure();
    fig2=gcf;
    title('critical edges')
    axis equal; hold on;
    axis(axesLim);
    if maximise set(gcf, 'Position', fact*get(0,'Screensize')); end % Maximize figure.
    
    if histo
        figure();
        fig3=gcf;
        title('length distribution of cell segments (a segment that is associated to two cells is counted twice)');
        hold on;
        if maximise set(gcf, 'Position',fact*get(0,'Screensize')); end % Maximize figure.
    end
    
    figure();
    fig4=gcf;
    title('Structure after removal of short segments');
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
if flagPlotSteps
    figure();
    fig6=gcf;
    title(['Structure after Iteration ' num2str(itCount)]);
    axis equal; hold on;
    axis(axesLim);
    if maximise set(gcf, 'Position', fact*get(0,'Screensize')); end % Maximize figure.
end


segmentNorms=[];
vertX_vertY=[];

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

%handle short segments
%========================
replacements=true;
segmentsItCount=0;
numEdgesReplacedTot=0;

while replacements && segmentsItCount<maxIterations
    if flagPlotDetails
        figure(fig2);
        cla;
        for i=1:numCells
            plotCell(cells(i),vertices,'-b');
            %figure 2 is for highlighting critical edges
        end
    end
    
    
    numEdgesReplaced=0;
    for i=1:numCells
        
        
        currCellVert=cells{i};
        newCurrCellVert=cells{i};
        if size(currCellVert,1)>3
            vertX_vertY__temp=[vertices([cells{i,1}';cells{i,1}(1)'],1),vertices([cells{i,1}';cells{i,1}(1)'],2)];
            %a matrix that represents the segments of the i-th Cell as relative vectors
            %(on the rows of the matrix):
            segmX_segmY__temp=differenceMatrix(vertX_vertY__temp);
            %get lengths of segments:
            normTemp=(segmX_segmY__temp(:,1).^2+segmX_segmY__temp(:,2).^2).^0.5;
            %segmX_segmY=[segmX_segmY;segmX_segmY__temp];
            vertX_vertY=[vertX_vertY;vertX_vertY__temp];
            segmentNorms=[segmentNorms;normTemp];
            
            %highlight critical segments and cell
            mySize=size(normTemp);
            numEdge=mySize(1);
            j=1;
            
            replacement=false;
            %remove maximum one segment per cell
            while j<=numEdge && ~replacement
                %disp(['cellVert= ' num2str(currCellVert)]);
                %disp(['j (edge or vertice)= ' num2str(j)]);
                %disp(['vertice index= ' num2str(currCellVert(j))]);
                %disp(['nomVert= ' num2str(nomVert)]);
                
                if normTemp(j)<critLength
                    replacement=true;
                    if j==numEdge
                        
                        vertInd1=cells{i}(j);
                        vertInd2=cells{i}(1);
                    else
                        vertInd1=cells{i}(j) ;
                        vertInd2=cells{i}(j+1);
                    end
                    replacementVertice=(vertices(vertInd1,:)+vertices(vertInd2,:))/2;
                    vertices=[vertices;replacementVertice];
                    indRepVert=size(vertices,1);
                    
                    if flagPlotDetails
                        figure(fig2);
                        plot(vertices([cells{i,1}';cells{i,1}(1)'],1),vertices([cells{i,1}';cells{i,1}(1)'],2),'-c')
                        plot(vertices([vertInd1;vertInd2],1),vertices([vertInd1;vertInd2],2),'lineStyle','-','color','r','lineWidth',4)
                        plot(replacementVertice(1,1),replacementVertice(1,2),'db');
                    end
                    
                    %disp('replacing!');
                    if j==numEdge
                        newCurrCellVert(j)=indRepVert;
                        newCurrCellVert(1)=-1;
                        %end loop
                    else
                        newCurrCellVert(j)=indRepVert;
                        newCurrCellVert(j+1)=-1;
                        %end loop
                    end
                    newCurrCellVert=newCurrCellVert((newCurrCellVert==-1)==false);
                    cells{i}=newCurrCellVert;
                    
                    for l=1:numCells
                        currCellVert=cells{l};
                        newCurrCellVert=cells{l};
                        numVert=size(currCellVert,2);
                        %replaced2=false;
                        m=1;
                        while m<=numVert %&& ~replaced2
                            if m==numVert
                                if sum([currCellVert(m) currCellVert(1)]==[vertInd1 vertInd2])==2 ||sum([currCellVert(m) currCellVert(1)]==[vertInd2 vertInd1])==2
                                    newCurrCellVert(m)=indRepVert;
                                    newCurrCellVert(1)=-1;
                                elseif currCellVert(m)==vertInd1 || currCellVert(m)==vertInd2
                                    newCurrCellVert(m)=indRepVert;
                                end
                            else
                                if sum(currCellVert(m:m+1)==[vertInd1 vertInd2])==2 || sum(currCellVert(m:m+1)==[vertInd2 vertInd1])==2
                                    newCurrCellVert(m)=indRepVert;
                                    newCurrCellVert(m+1)=-1;
                                elseif currCellVert(m)==vertInd1 || currCellVert(m)==vertInd2
                                    newCurrCellVert(m)=indRepVert;
                                end
                            end
                            m=m+1;
                        end
                        newCurrCellVert=newCurrCellVert((newCurrCellVert==-1)==false);
                        cells{l}=newCurrCellVert;
                    end
                    
                    
                    
                    
                    
                    numEdgesReplaced=numEdgesReplaced+1;
                    numEdgesReplacedTot=numEdgesReplacedTot+numEdgesReplaced;
                end
                
                
                
                
                
                j=j+1;
            end
            
        else
            warnTriangle=true;
        end
        
        
        
        
    end
    
    if numEdgesReplaced== 0, replacements=false; end
    segmentsItCount=segmentsItCount+1;
end
disp(['removal of short segments took ' num2str(segmentsItCount) ' iterations']);

%plot length distribution histogram
if flagPlotDetails && histo
    %forget about ridiculously long segments:
    segmentNorms=segmentNorms((segmentNorms==inf)==false);
    stdNorms=std(segmentNorms);
    meanNorms=mean(segmentNorms);
    segmentNorms=segmentNorms(segmentNorms<meanNorms+0.5*stdNorms);
    
    %plot length distribution histogram:
    figure(fig3);
    h=histogram(segmentNorms,40);
    hold on;
    lim=axis;
    %plot vertical red line to show critical length value
    plot([critLength,critLength],[0,lim(4)],'LineWidth',2,'lineStyle','-','color','r');
end
[vertices,cells]=cleanStructure(vertices,cells);


%==============
%plot structure after cleanout of small segments
if flagPlotDetails
    for i=1:numCells
        
        
        figure(fig4);
        plotCell(cells(i),vertices,'-b');
        %plot(repVertices(:,1),repVertices(:,2),'db');
    end
    %disp('next cell');
    
end



%%

%Angles
%====================================
%
%find critical angles
angleItCount=0;

while (angleItCount<=maxIterations && (numEdgesReplaced>0 || numAnglesReplaced>0))
    
    
    numAnglesReplaced=0;
    newVertices=vertices;
    newCells=cells;
    for i=1:numCells
        if flagPlotDetails
            figure(fig5);
            plotCell(newCells(i),newVertices,'-b');
        end
        
        tempArr=cells{i};
        mySize=size(tempArr);
        Nvert=mySize(2);
        jBis=1;
        
        if Nvert>3 %check if cell is not a triangle
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
                a=[a 0];
                b=[b 0];
                angle=atan2(norm(cross(a,b)),dot(a,b));
                if angle <= critAngle
                    
                    
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
                    
                    if ~(sum(isnan(P1))>0) && ~(sum(isnan(P2))>0)
                        newVertices=[newVertices;P1;P2];
                        mySize=size(newVertices);
                        numVert=mySize(1);
                        currCellVert=newCells{i};
                        currCellVert=[currCellVert(1:jBis-1) numVert-1 numVert currCellVert(jBis+1:end)];
                        plotTri=[P1;vert2;P2];
                        jBis=jBis+1;
                        
                    elseif (sum(isnan(P1))>0) && ~(sum(isnan(P2))>0)
                        newVertices=[newVertices;P2];
                        mySize=size(newVertices);
                        numVert=mySize(1);
                        currCellVert=newCells{i};
                        currCellVert=[currCellVert(1:jBis-1) numVert currCellVert(jBis+1:end)];
                        plotTri=[vert1;vert2;P2];
                        
                        
                    elseif ~(sum(isnan(P1))>0) && (sum(isnan(P2))>0)
                        newVertices=[newVertices;P1];
                        mySize=size(newVertices);
                        numVert=mySize(1);
                        currCellVert=newCells{i};
                        currCellVert=[currCellVert(1:jBis-1) numVert currCellVert(jBis+1:end)];
                        plotTri=[P1;vert2;vert3];
                    elseif (sum(isnan(P1))>0) && (sum(isnan(P2))>0)
                        currCellVert=newCells{i};
                        currCellVert=[currCellVert(1:jBis-1) currCellVert(jBis+1:end)];
                        plotTri=[vert1;vert2;vert3];
                        jBis=jBis-1;
                    end
                    
                    newCells{i}=currCellVert;
                    %=========================================
                    if flagPlotDetails
                        figure(fig5);
                        fill(plotTri(:,1),plotTri(:,2),'r');
                    end
                    
                end
                jBis=jBis+1;
            end
        else
            warnTriangle=true; %cell is a triangle
        end
        
        
        
        
        
        
        
    end
    
    cells=newCells;
    vertices=newVertices;
    
    
    [vertices,cells]=cleanStructure(vertices,cells);
    
    for i=1:numCells
        
        if flagPlotSteps
            figure(fig6);
            plotCell(cells(i),vertices,'-b');
        end
    end
    
    disp(['End of refinement iteration ' num2str(itCount)]);
    disp(['Number of short edges replaced by angle: ' num2str(numEdgesReplacedTot)]);
    disp(['Number of critical angles replaced by edge: ' num2str(numAnglesReplaced)]);
    
    itCount=itCount+1;
    
    angleItCount=angleItCount+1;
end




if (numEdgesReplaced==0 && numAnglesReplaced==0)
    disp('refinement algorithm converged successfully.');
else
    disp(['Output -structure contains critical angles and/or small edges, consider allowing more Iterations. Stopped after ' num2str(itCount-1) ' iteration(s).']);
end

if warnTriangle
    warning('triangular cells could not be treated, might contain critical angles/edges, considerr allowing larger critical Length')
end



outVertices=vertices;
outCells=cells;

if flagPlotFinal
    figure();
    title('Output Structure');
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