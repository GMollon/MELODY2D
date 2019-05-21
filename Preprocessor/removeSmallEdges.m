function [ outVertices,outCells ] = removeSmallEdges(vertices,cells,critLength,flagPlotDetails,flagPlotFinal,boundaries,maxIterations)
%removeSmallEdges Modifies a Voronoi-Cell microstructure as to eliminate
%  edges that are to small.
%   Input arguments 'vertices' and 'cells' represent a 2D cellular
%   structure as created e.g by the function voronoin. 'critLength' is the
%   target value for minimum edge length (of any cell). Plot Options (logicals) allow
%   to illustrate the progress of the refinement. The Refinement is done by
%   an iterative algorithm, set 'maxIterations' to limit calculation time.
%   Please mind that in the case of interruption by 'maxIterations' there
%   are still critical angles and or edges present (during runtime, see command window
%   messages for more information about the results)
%
%
%
%   last modified: 15.11.2016
%   by: Simon Massa
axesLim=boundaries;%[0.44 0.6 0.7 0.86];
%



%====================================================
%Begin refinement
disp('Refine Voronoi Structure: Remove small edges');
disp(['min Length of Edge: ' num2str(critLength)]);
%disp(['max iterations allowed: ' num2str(maxIterations)]);

%internal parameters:
safety=1.05;
histo=false;
maximise=true;

%initializations
itCount=1;
numEdgesReplaced=1;
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

if flagPlotDetails
    figure(fig2);
    cla;
    for i=1:numCells
        plotCell(cells(i),vertices,'-b');
        %figure 2 is for highlighting critical edges
    end
end

while replacements && segmentsItCount<maxIterations
    
    
    numEdgesReplaced=0;
    for i=1:numCells
        %only if at least one vertice is in boundarys
        currCellVert=cells{i};
        newCurrCellVert=cells{i};

        cellIsInArea=false;
        Nvert=size(currCellVert,1);
        for j=1:Nvert
            if isInArea(vertices(currCellVert(j),:),boundaries);
                cellIsInArea=true;
            end
        end
        
        if true%cellIsInArea
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
                if size(currCellVert,2)>3
                    replacement=true;
                    replacements=true;
                    
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
                    
                else
                    warnTriangle=true;
                end
            end
            
            
            
            
            
            j=j+1;
        end
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










if numEdgesReplaced==0
    disp('Algorithm converged successfully. Small Edges removed');
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
    title('Output Structure (small Edges removed)');
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