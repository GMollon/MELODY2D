function [INITIAL_POSITIONS,BORDERS,TRIANGULATION]=SetStructuredNodesForClosedContour(Contour,Interpolant,Average_Nodal_Distance)
if norm(Contour(1,:)-Contour(size(Contour,1),:))<1e-8*Average_Nodal_Distance
    Contour=Contour(1:size(Contour,1)-1,:);
end
INITIAL_POSITIONS=Contour;
BORDERS=[];

% Boundary Field nodes
for i=1:size(Contour,1)
    if i==size(Contour,1)
        Nini=i;
        Nfin=1;
    else
        Nini=i;
        Nfin=i+1;
    end
    Xini=Contour(Nini,1);
    Yini=Contour(Nini,2);
    Xfin=Contour(Nfin,1);
    Yfin=Contour(Nfin,2);
    Nnewnodes=floor(norm([Xfin-Xini,Yfin-Yini])/Average_Nodal_Distance)-1;
    if Nnewnodes>0
        if abs(Xini-Xfin)<1e-8*Average_Nodal_Distance
            XNewnodes=Xini*ones(Nnewnodes,1);
        else
            XNewnodes=[Xini+0.999999*(Xfin-Xini)/(Nnewnodes+1):0.999999*(Xfin-Xini)/(Nnewnodes+1):Xfin-0.999999*(Xfin-Xini)/(Nnewnodes+1)]';
        end
        if abs(Yini-Yfin)<1e-8*Average_Nodal_Distance
            YNewnodes=Yini*ones(Nnewnodes,1);
        else
            YNewnodes=[Yini+0.999999*(Yfin-Yini)/(Nnewnodes+1):0.999999*(Yfin-Yini)/(Nnewnodes+1):Yfin-0.999999*(Yfin-Yini)/(Nnewnodes+1)]';
        end
        NewNodes=[size(INITIAL_POSITIONS,1)+1:size(INITIAL_POSITIONS,1)+Nnewnodes]';
        NewSegments=cat(2,NewNodes,NewNodes+1);
        NewSegments=cat(1,[Nini,NewSegments(1,1)],NewSegments);
        NewSegments(size(NewSegments,1),2)=Nfin;
        INITIAL_POSITIONS=cat(1,INITIAL_POSITIONS,[XNewnodes,YNewnodes]);
        BORDERS=cat(1,BORDERS,NewSegments);
    else
        NewSegments=[Nini,Nfin];
        BORDERS=cat(1,BORDERS,NewSegments);
    end
end

% Fill with nodes and triangulate
Xmin=min(Contour(:,1));Xmax=max(Contour(:,1));
Ymin=min(Contour(:,2));Ymax=max(Contour(:,2));
Nxinter=round((Xmax-Xmin)/Average_Nodal_Distance);
Nyinter=round((Ymax-Ymin)/Average_Nodal_Distance);
Xinter=(Xmax-Xmin)/Nxinter;
Yinter=(Ymax-Ymin)/Nyinter;
for i=Xmin+Xinter:Xinter:Xmax-Xinter
    for j=Ymin+Yinter:Yinter:Ymax-Yinter
        Xcandidate=i;
        Ycandidate=j;
        if inpolygon(Xcandidate,Ycandidate,Contour(:,1),Contour(:,2))
            add=1;
            for j=1:size(BORDERS,1)
                xp=INITIAL_POSITIONS(BORDERS(j,1),1);
                yp=INITIAL_POSITIONS(BORDERS(j,1),2);
                if sqrt((Xcandidate-xp)^2+(Ycandidate-yp)^2)<Average_Nodal_Distance/2
                    add=0;
                    break
                end
            end
            if add==1
                INITIAL_POSITIONS=cat(1,INITIAL_POSITIONS,[Xcandidate,Ycandidate]);
            end
        end
    end
end
TRIANGULATION=delaunay(INITIAL_POSITIONS(:,1),INITIAL_POSITIONS(:,2));
BORDERS={'Closed',size(BORDERS,1),BORDERS(:,1),Interpolant};