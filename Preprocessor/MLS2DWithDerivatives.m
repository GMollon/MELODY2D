function [Shape_Functions,Weight_Functions]=MLS2DWithDerivatives(Poi_Position,Fnodes_Positions,List_Domain_Fnodes,Domains_Dimensions,Number_Monomials)
%x=Poi_Position(1);
%y=Poi_Position(2);
%xi=Fnodes_Positions(List_Domain_Fnodes,1);
%yi=Fnodes_Positions(List_Domain_Fnodes,2);
%c=Domains_Dimensions(List_Domain_Fnodes,2);
%d=Domains_Dimensions(List_Domain_Fnodes,1);
%
%W=(exp(-((x-xi).^2+(y-yi).^2)./c.^2)-exp(-(d./c).^2))./(1-exp(-(d./c).^2));
%dW=exp(-((x-xi).^2+(y-yi).^2)./c.^2)./(1-exp(-(d./c).^2));
%dWdx=-2*(x-xi)./(c.^2).*dW;
%dWdy=-2*(y-yi)./(c.^2).*dW;
%Weight_Functions=[W,dWdx,dWdy];
%W=diag(W);
%dWdx=diag(dWdx);
%dWdy=diag(dWdy);
%
%P=[ones(1,size(xi,1));xi';yi'];%;(xi.^2)';(xi.*yi)';(yi.^2)';(xi.^3)';(xi.^2.*yi)';(xi.*yi.^2)';(yi.^3)'];
%dPdx=[zeros(1,size(xi,1));ones(1,size(xi,1));zeros(1,size(xi,1));2*xi';yi';zeros(1,size(xi,1));3*(xi.^2)';(2*xi.*yi)';(yi.^2)';zeros(1,size(xi,1))];
%dPdy=[zeros(1,size(xi,1));zeros(1,size(xi,1));ones(1,size(xi,1));zeros(1,size(xi,1));xi';2*yi';zeros(1,size(xi,1));(xi.^2)';(2*xi.*yi)';3*(yi.^2)'];
%p=[1;x;y];%;x^2;x*y;y^2;x^3;x^2*y;x*y^2;y^3];
%dpdx=[0;1;0];%;2*x;y;0;3*x^2;2*x*y;y^2;0];
%dpdy=[0;0;1];%;0;x;2*y;0;x^2;2*x*y;3*y^2];
%
%P=P(1:Number_Monomials,:);
%dPdx=dPdx(1:Number_Monomials,:);
%dPdy=dPdy(1:Number_Monomials,:);
%p=p(1:Number_Monomials,:);
%dpdx=dpdx(1:Number_Monomials,:);
%dpdy=dpdy(1:Number_Monomials,:);
%
%B=P*W;
%A=B*P';
%invA=inv(A);
%
%dBdx=dPdx*W+P*dWdx;
%dBdy=dPdy*W+P*dWdy;
%dAdx=dBdx*P'+B*dPdx';
%dAdy=dBdy*P'+B*dPdy';
%dinvAdx=-invA*dAdx*invA;
%dinvAdy=-invA*dAdy*invA;
%
%Phi=p'*invA*B;
%dPhidx=dpdx'*invA*B+p'*dinvAdx*B+p'*invA*dBdx;
%dPhidy=dpdy'*invA*B+p'*dinvAdy*B+p'*invA*dBdy;
%Shape_Functions=[Phi',dPhidx',dPhidy'];
%return
%
%Compute MLS shape function and their derivatives
%
%Input variables :
%    Poi_Position (gpos) : coordinates x and y of the point of interest
%    Fnodes_Positions (x) : x and y coordinates for all field nodes
%    List_Domain_Fnodes (Nv) : Field nodes used in the support domain
%    Domains_Dimensions (Ds) : x and y dimensions of the support domain of each field node
%    Number_Monomials (mm) : number of monomials used in the basis
%
%Output variables :
%    Shape_Functions, with the following lines :
%        phi
%        dphix
%        dphiy
%

Support_Type=1;

Number_Fnodes=size(Fnodes_Positions,1);
Number_Fnodes_Domain=size(List_Domain_Fnodes,1);

% Compute polynomial basis
GP=zeros(6,6);
GP(1,1)=1;
GP(1,2)=Poi_Position(1);
GP(1,3)=Poi_Position(2);
GP(1,4)=Poi_Position(1)^2;
GP(1,5)=Poi_Position(1)*Poi_Position(2);
GP(1,6)=Poi_Position(2)^2;
GP(2,2)=1;
GP(2,4)=2*Poi_Position(1);
GP(2,5)=Poi_Position(2);
GP(3,3)=1;
GP(3,5)=Poi_Position(1);
GP(3,6)=2*Poi_Position(2);
GP(4,4)=2;
GP(5,5)=1;
GP(6,6)=2;

% Compute matrices A and B
Number_Fnodes_Domain=size(List_Domain_Fnodes,1);
Poi_Relative_Position=zeros(Number_Fnodes_Domain,2);
Moments=zeros(6,Number_Fnodes_Domain);
A=zeros(Number_Monomials,Number_Monomials,3);
B=zeros(Number_Monomials,Number_Fnodes_Domain,3);
Moments(1,:)=ones(1,Number_Fnodes_Domain);
Moments(2,:)=Fnodes_Positions(List_Domain_Fnodes,1)';
Moments(3,:)=Fnodes_Positions(List_Domain_Fnodes,2)';
Moments(4,:)=Fnodes_Positions(List_Domain_Fnodes,1)'.^2;
Moments(5,:)=Fnodes_Positions(List_Domain_Fnodes,1)'.*Fnodes_Positions(List_Domain_Fnodes,2)';
Moments(6,:)=Fnodes_Positions(List_Domain_Fnodes,2)'.^2;
Poi_Relative_Position(:,1)=Poi_Position(1)-Fnodes_Positions(List_Domain_Fnodes,1);
Poi_Relative_Position(:,2)=Poi_Position(2)-Fnodes_Positions(List_Domain_Fnodes,2);
W=WeightW1WithDerivatives(Poi_Relative_Position,List_Domain_Fnodes,Domains_Dimensions);
for ii=1:Number_Monomials
    for jj=1:Number_Fnodes_Domain
        for kk=1:3
            B(ii,jj,kk)=Moments(ii,jj)*W(jj,kk);
        end
    end
end
pp=zeros(Number_Monomials,Number_Monomials);
for iii=1:Number_Fnodes_Domain
    for ik=1:Number_Monomials
        for jk=1:Number_Monomials
            Moments2D(ik,jk)=Moments(ik,iii)*Moments(jk,iii);
        end
    end
    for ikk=1:Number_Monomials
        for jkk=1:Number_Monomials
            for kkk=1:3
                A(ikk,jkk,kkk)=A(ikk,jkk,kkk)+W(iii,kkk)*Moments2D(ikk,jkk);
            end
        end
    end
end
c=GP(1,1:Number_Monomials)';
Am=A(:,:,1)^-1;

gamma=zeros(Number_Monomials,3);

% Compute gamma
gamma(:,1)=Am*c;

% Compute dgamma/dx
c=zeros(Number_Monomials,1);
for in=1:Number_Monomials
    for jn=1:Number_Monomials
        c(in)=c(in)+A(in,jn,2)*gamma(jn,1);
    end
end
for kn=1:Number_Monomials
    c(kn)=GP(2,kn)-c(kn);
end
gamma(:,2)=Am*c;

% Compute dgamma/dy
c=zeros(Number_Monomials,1);
for in=1:Number_Monomials
    for jn=1:Number_Monomials
        c(in)=c(in)+A(in,jn,3)*gamma(jn,1);
    end
end
for kn=1:Number_Monomials
    c(kn)=GP(3,kn)-c(kn);
end
gamma(:,3)=Am*c;

% Compute Phi and derivatives
Shape_Functions=zeros(Number_Fnodes_Domain,3);
for iph=1:Number_Fnodes_Domain
    for jph=1:Number_Monomials
        Shape_Functions(iph,1)=Shape_Functions(iph,1)+gamma(jph,1)*B(jph,iph,1);
        Shape_Functions(iph,2)=Shape_Functions(iph,2)+gamma(jph,2)*B(jph,iph,1)+gamma(jph,1)*B(jph,iph,2);
        Shape_Functions(iph,3)=Shape_Functions(iph,3)+gamma(jph,3)*B(jph,iph,1)+gamma(jph,1)*B(jph,iph,3);
    end
end
Weight_Functions=W;