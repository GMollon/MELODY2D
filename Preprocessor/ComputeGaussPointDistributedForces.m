function Local_Distributed_Forces=ComputeGaussPointDistributedForces(Shape_Functions,Weight,Jacobian,Gravity,Rho)
%
%Provides the local distributed forces corresponding to a given Gauss point
%
%Input variables :
%    Shape_Functions : Values and derivatives (at the Gauss point) of the
%        shape functions of the relevant field nodes
%    Weight : weight of the Gauss points
%    Jacobian : Jacobian of the affine transform at the Gauss point
%    Gravity : vector of the gravity field [gx,gy]
%    Rho : local unit mass of the material at the Gauss point
%
%Output variables :
%    Local_Distributed_Forces : Distributed forces vector for the Gauss point
%
%Nfnodes=size(Shape_Functions,1);
%Phimat=zeros(2,2*Nfnodes);
%for in=1:Nfnodes
%    Phimat(1,2*in-1)=Shape_Functions(in,1);
%    Phimat(2,2*in)=Shape_Functions(in,1);
%end
%
%Local_Distributed_Forces=zeros(2*Nfnodes,1);
%for ii=1:2*Nfnodes
%    for kk=1:2
%        Local_Distributed_Forces(ii,1)=Local_Distributed_Forces(ii,1)+Weight*Jacobian*Phimat(kk,ii)*Rho*Gravity(kk);
%    end
%end

Nfnodes=size(Shape_Functions,1);
Local_Distributed_Forces=zeros(2*Nfnodes,1);
for ii=1:Nfnodes
    Local_Distributed_Forces(2*ii-1,1)=Weight*Jacobian*Shape_Functions(ii,1)*Rho*Gravity(1);
    Local_Distributed_Forces(2*ii,1)=Weight*Jacobian*Shape_Functions(ii,1)*Rho*Gravity(2);
end
