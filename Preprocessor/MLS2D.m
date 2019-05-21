function Shape_Functions=MLS2D(Poi_Position,Fnodes_Positions,List_Domain_Fnodes,Domains_Dimensions)
Number_Fnodes=size(Fnodes_Positions,1);
Number_Fnodes_Domain=size(List_Domain_Fnodes,1);

% Compute polynomial basis
GP=zeros(3,3);
GP(1,1)=1;
GP(1,2)=Poi_Position(1);
GP(1,3)=Poi_Position(2);
GP(2,2)=1;
GP(3,3)=1;

% Compute matrices A and B
Number_Fnodes_Domain=size(List_Domain_Fnodes,1);
Poi_Relative_Position=zeros(Number_Fnodes_Domain,2);
Moments=zeros(3,Number_Fnodes_Domain);
A=zeros(3,3);
B=zeros(3,Number_Fnodes_Domain);
Moments(1,:)=ones(1,Number_Fnodes_Domain);
Moments(2,:)=Fnodes_Positions(List_Domain_Fnodes,1)';
Moments(3,:)=Fnodes_Positions(List_Domain_Fnodes,2)';
Poi_Relative_Position(:,1)=Poi_Position(1)-Fnodes_Positions(List_Domain_Fnodes,1);
Poi_Relative_Position(:,2)=Poi_Position(2)-Fnodes_Positions(List_Domain_Fnodes,2);
W=WeightW1(Poi_Relative_Position,List_Domain_Fnodes,Domains_Dimensions);
for ii=1:3
    for jj=1:Number_Fnodes_Domain
        B(ii,jj)=Moments(ii,jj)*W(jj);
    end
end
for iii=1:Number_Fnodes_Domain
    for ik=1:3
        for jk=1:3
            Moments2D(ik,jk)=Moments(ik,iii)*Moments(jk,iii);
        end
    end
    for ikk=1:3
        for jkk=1:3
            A(ikk,jkk)=A(ikk,jkk)+W(iii)*Moments2D(ikk,jkk);
        end
    end
end
c=GP(1,1:3)';
Am=A(:,:,1)^-1;

% Compute gamma
gamma(:,1)=Am*c;

% Compute Phi and derivatives
Shape_Functions=zeros(Number_Fnodes_Domain,1);
for iph=1:Number_Fnodes_Domain
    for jph=1:3
        Shape_Functions(iph,1)=Shape_Functions(iph,1)+gamma(jph)*B(jph,iph,1);
    end
end