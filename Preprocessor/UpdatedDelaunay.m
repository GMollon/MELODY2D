function Integration_Cells=UpdatedDelaunay(Tri,X,Y)

Nodes_Cells=cell(size(X,1),1);
for i=1:size(Tri,1)
    Nodes_Cells{Tri(i,1),1}=cat(1,Nodes_Cells{Tri(i,1),1},i);
    Nodes_Cells{Tri(i,2),1}=cat(1,Nodes_Cells{Tri(i,2),1},i);
    Nodes_Cells{Tri(i,3),1}=cat(1,Nodes_Cells{Tri(i,3),1},i);
end

Nodes_Segments=cell(size(X,1),1);
Integration_Cells=cell(size(X,1),1);
for i=1:size(Nodes_Segments,1)
    for j=1:size(Nodes_Cells{i,1},1)
        N1=Tri(Nodes_Cells{i,1}(j,1),1);
        N2=Tri(Nodes_Cells{i,1}(j,1),2);
        N3=Tri(Nodes_Cells{i,1}(j,1),3);
        if N1==i
            Nodes_Segments{i,1}=cat(1,Nodes_Segments{i,1},[N2,N3]);
        elseif N2==i
            Nodes_Segments{i,1}=cat(1,Nodes_Segments{i,1},[N1,N3]);
        elseif N3==i
            Nodes_Segments{i,1}=cat(1,Nodes_Segments{i,1},[N1,N2]);
        end
    end
    
    List=sort(cat(1,Nodes_Segments{i,1}(:,1),Nodes_Segments{i,1}(:,2)));
    Flag_Boundary=0;
    for j=1:2:size(List,1)
        if List(j)~=List(j+1)
            Flag_Boundary=1;
            First_Node=List(j);
            for k=size(List,1):-2:1
                if List(k)~=List(k-1)
                    Last_Node=List(k);
                    break
                end
            end
            break
        end
    end
    if Flag_Boundary==0
        Cell=zeros(size(List,1)+1,2);
        First_Node=List(1);
        Starting_Node=First_Node;
        Previous_Node=0;
        Filled=0;
        Flag_Exit=0;
        while Flag_Exit==0
            for j=1:size(Nodes_Segments{i,1},1)
                if Nodes_Segments{i,1}(j,1)==Starting_Node & Nodes_Segments{i,1}(j,2)~=Previous_Node
                    Ending_Node=Nodes_Segments{i,1}(j,2);
                    break
                elseif Nodes_Segments{i,1}(j,2)==Starting_Node & Nodes_Segments{i,1}(j,1)~=Previous_Node
                    Ending_Node=Nodes_Segments{i,1}(j,1);
                    break
                end
            end
            Cell(Filled+1,:)=[(X(Starting_Node)+X(i))/2,(Y(Starting_Node)+Y(i))/2];
            Cell(Filled+2,:)=[(X(Starting_Node)+X(i)+X(Ending_Node))/3,(Y(Starting_Node)+Y(i)+Y(Ending_Node))/3];
            Filled=Filled+2;
            Previous_Node=Starting_Node;
            Starting_Node=Ending_Node;
            if Starting_Node==First_Node
                Cell(Filled+1,:)=[(X(Starting_Node)+X(i))/2,(Y(Starting_Node)+Y(i))/2];
                Flag_Exit=1;
            end
        end
    elseif Flag_Boundary==1
        Cell=zeros(size(List,1)+3,2);
        Cell(1,:)=[X(i),Y(i)];
        Cell(size(Cell,1),:)=[X(i),Y(i)];
        Starting_Node=First_Node;
        Previous_Node=0;
        Filled=1;
        Flag_Exit=0;
        while Flag_Exit==0
            for j=1:size(Nodes_Segments{i,1},1)
                if Nodes_Segments{i,1}(j,1)==Starting_Node & Nodes_Segments{i,1}(j,2)~=Previous_Node
                    Ending_Node=Nodes_Segments{i,1}(j,2);
                    break
                elseif Nodes_Segments{i,1}(j,2)==Starting_Node & Nodes_Segments{i,1}(j,1)~=Previous_Node
                    Ending_Node=Nodes_Segments{i,1}(j,1);
                    break
                end
            end
            Cell(Filled+1,:)=[(X(Starting_Node)+X(i))/2,(Y(Starting_Node)+Y(i))/2];
            Cell(Filled+2,:)=[(X(Starting_Node)+X(i)+X(Ending_Node))/3,(Y(Starting_Node)+Y(i)+Y(Ending_Node))/3];
            Filled=Filled+2;
            Previous_Node=Starting_Node;
            Starting_Node=Ending_Node;
            if Starting_Node==Last_Node
                Cell(Filled+1,:)=[(X(Starting_Node)+X(i))/2,(Y(Starting_Node)+Y(i))/2];
                Flag_Exit=1;
            end
        end
        %%%%
        %Cell=Cell(2:size(Cell,1)-1,:);
        %%%%
    end
    
    %[Theta,Rho] = cart2pol(Cell(:,1)-X(i),Cell(:,2)-Y(i));
    [Theta,Rho] = cart2pol(Cell(:,1)-mean(Cell(:,1)),Cell(:,2)-mean(Cell(:,2)));
    Pos_xmax=find(Cell(:,1)==max(Cell(:,1)));
    Pos_xmax=Pos_xmax(1);
    Flag_Reverse=0;
    if Pos_xmax==size(Cell,1)
        if Theta(Pos_xmax,1)<Theta(Pos_xmax-1,1)
            Flag_Reverse=1;
        end
    else
        if Theta(Pos_xmax+1,1)<Theta(Pos_xmax,1)
            Flag_Reverse=1;
        end
    end
    if Flag_Reverse==1
        Cell=flipud(Cell);
    end
    Integration_Cells{i,1}=Cell;
end





return
figure;
for i=1:size(Nodes_Segments,1)
    if size(Integration_Cells{i,1},1)==0
        continue
    end
    for j=1:size(Nodes_Segments{i,1},1)
        n1=Nodes_Segments{i,1}(j,1);
        n2=Nodes_Segments{i,1}(j,2);
        x1=X(n1,1);
        y1=Y(n1,1);
        x2=X(n2,1);
        y2=Y(n2,1);
        %plot([x1,x2],[y1,y2],'-b')
        hold on
    end
    plot(Integration_Cells{i,1}(:,1),Integration_Cells{i,1}(:,2),'-r')
end
plot(X,Y,'.r','markersize',8)
axis equal