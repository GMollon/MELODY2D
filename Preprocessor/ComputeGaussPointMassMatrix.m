function Local_Mass_Matrix=ComputeGaussPointMassMatrix(Shape_Functions,Weight,Jacobian,Rho)
Nfnodes=size(Shape_Functions,1);
Local_Mass_Matrix=zeros(2*Nfnodes,2*Nfnodes);
% for ii=1:Nfnodes
%     %for jj=1:Nfnodes
%     %    Local_Mass_Matrix(2*ii-1,2*jj-1)=Weight*Jacobian*Rho*Shape_Functions(ii,1)*Shape_Functions(jj,1);
%     %    Local_Mass_Matrix(2*ii,2*jj)=Weight*Jacobian*Rho*Shape_Functions(ii,1)*Shape_Functions(jj,1);
%     %end
%     %
%     jj=[1:Nfnodes];
%     Local_Mass_Matrix(2*ii-1,2*jj-1)=Weight*Jacobian*Rho*Shape_Functions(ii,1)*Shape_Functions(jj,1);
%     Local_Mass_Matrix(2*ii,2*jj)=Weight*Jacobian*Rho*Shape_Functions(ii,1)*Shape_Functions(jj,1);
%     %
% end
%
ii=[1:Nfnodes]';
jj=[1:Nfnodes];
Local_Mass_Matrix(2*ii-1,2*jj-1)=Weight*Jacobian*Rho*Shape_Functions(ii,1)'.*Shape_Functions(jj,1);
Local_Mass_Matrix(2*ii,2*jj)=Weight*Jacobian*Rho*Shape_Functions(ii,1)'.*Shape_Functions(jj,1);