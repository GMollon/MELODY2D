function Local_Stiffness_Matrix=ComputeGaussPointStiffnessMatrix(Shape_Functions,Weight,Jacobian,Elasticity_Matrix)
%
%Provides the local stiffness matrix corresponding to a given Gauss point
%
%Input variables :
%    Shape_Functions : Values and derivatives (at the Gauss point) of the
%        shape functions of the relevant field nodes
%    Weight : weight of the Gauss points
%    Jacobian : Jacobian of the affine transform at the Gauss point
%    Elasticity_Matrix : Matrix of the local elasticity tensor
%
%Output variables :
%    Local_Matrix : Stiffness matrix for the Gauss point
%
Nfnodes=size(Shape_Functions,1);
Local_Stiffness_Matrix=zeros(2*Nfnodes,2*Nfnodes);
% for ii=1:Nfnodes
%     i=2*ii-1;
%     for jj=1:Nfnodes
%         j=2*jj-1;
%         Local_Stiffness_Matrix(i,j)=Local_Stiffness_Matrix(i,j)+Shape_Functions(ii,2)*Elasticity_Matrix(1,1)*Shape_Functions(jj,2);
%         Local_Stiffness_Matrix(i,j)=Local_Stiffness_Matrix(i,j)+Shape_Functions(ii,2)*Elasticity_Matrix(1,3)*Shape_Functions(jj,3);
%         Local_Stiffness_Matrix(i,j)=Local_Stiffness_Matrix(i,j)+Shape_Functions(ii,3)*Elasticity_Matrix(3,1)*Shape_Functions(jj,2);
%         Local_Stiffness_Matrix(i,j)=Local_Stiffness_Matrix(i,j)+Shape_Functions(ii,3)*Elasticity_Matrix(3,3)*Shape_Functions(jj,3);
%     end
%     for jj=1:Nfnodes
%         j=2*jj;
%         Local_Stiffness_Matrix(i,j)=Local_Stiffness_Matrix(i,j)+Shape_Functions(ii,2)*Elasticity_Matrix(1,2)*Shape_Functions(jj,3);
%         Local_Stiffness_Matrix(i,j)=Local_Stiffness_Matrix(i,j)+Shape_Functions(ii,2)*Elasticity_Matrix(1,3)*Shape_Functions(jj,2);
%         Local_Stiffness_Matrix(i,j)=Local_Stiffness_Matrix(i,j)+Shape_Functions(ii,3)*Elasticity_Matrix(3,2)*Shape_Functions(jj,3);
%         Local_Stiffness_Matrix(i,j)=Local_Stiffness_Matrix(i,j)+Shape_Functions(ii,3)*Elasticity_Matrix(3,3)*Shape_Functions(jj,2);
%     end
% end
% for ii=1:Nfnodes
%     i=2*ii;
%     for jj=1:Nfnodes
%         j=2*jj-1;
%         Local_Stiffness_Matrix(i,j)=Local_Stiffness_Matrix(i,j)+Shape_Functions(ii,3)*Elasticity_Matrix(2,1)*Shape_Functions(jj,2);
%         Local_Stiffness_Matrix(i,j)=Local_Stiffness_Matrix(i,j)+Shape_Functions(ii,3)*Elasticity_Matrix(2,3)*Shape_Functions(jj,3);
%         Local_Stiffness_Matrix(i,j)=Local_Stiffness_Matrix(i,j)+Shape_Functions(ii,2)*Elasticity_Matrix(3,1)*Shape_Functions(jj,2);
%         Local_Stiffness_Matrix(i,j)=Local_Stiffness_Matrix(i,j)+Shape_Functions(ii,2)*Elasticity_Matrix(3,3)*Shape_Functions(jj,3);
%     end
%     for jj=1:Nfnodes
%         j=2*jj;
%         Local_Stiffness_Matrix(i,j)=Local_Stiffness_Matrix(i,j)+Shape_Functions(ii,3)*Elasticity_Matrix(2,2)*Shape_Functions(jj,3);
%         Local_Stiffness_Matrix(i,j)=Local_Stiffness_Matrix(i,j)+Shape_Functions(ii,3)*Elasticity_Matrix(2,3)*Shape_Functions(jj,2);
%         Local_Stiffness_Matrix(i,j)=Local_Stiffness_Matrix(i,j)+Shape_Functions(ii,2)*Elasticity_Matrix(3,2)*Shape_Functions(jj,3);
%         Local_Stiffness_Matrix(i,j)=Local_Stiffness_Matrix(i,j)+Shape_Functions(ii,2)*Elasticity_Matrix(3,3)*Shape_Functions(jj,2);
%     end
% end
%
ii=[1:Nfnodes]';
i=2*ii-1;
jj=[1:Nfnodes];
j=2*jj-1;
Local_Stiffness_Matrix(i,j)=Local_Stiffness_Matrix(i,j)+Shape_Functions(ii,2)*Elasticity_Matrix(1,1)*Shape_Functions(jj,2)';
Local_Stiffness_Matrix(i,j)=Local_Stiffness_Matrix(i,j)+Shape_Functions(ii,2)*Elasticity_Matrix(1,3)*Shape_Functions(jj,3)';
Local_Stiffness_Matrix(i,j)=Local_Stiffness_Matrix(i,j)+Shape_Functions(ii,3)*Elasticity_Matrix(3,1)*Shape_Functions(jj,2)';
Local_Stiffness_Matrix(i,j)=Local_Stiffness_Matrix(i,j)+Shape_Functions(ii,3)*Elasticity_Matrix(3,3)*Shape_Functions(jj,3)';
j=2*jj;
Local_Stiffness_Matrix(i,j)=Local_Stiffness_Matrix(i,j)+Shape_Functions(ii,2)*Elasticity_Matrix(1,2)*Shape_Functions(jj,3)';
Local_Stiffness_Matrix(i,j)=Local_Stiffness_Matrix(i,j)+Shape_Functions(ii,2)*Elasticity_Matrix(1,3)*Shape_Functions(jj,2)';
Local_Stiffness_Matrix(i,j)=Local_Stiffness_Matrix(i,j)+Shape_Functions(ii,3)*Elasticity_Matrix(3,2)*Shape_Functions(jj,3)';
Local_Stiffness_Matrix(i,j)=Local_Stiffness_Matrix(i,j)+Shape_Functions(ii,3)*Elasticity_Matrix(3,3)*Shape_Functions(jj,2)';
i=2*ii;
j=2*jj-1;
Local_Stiffness_Matrix(i,j)=Local_Stiffness_Matrix(i,j)+Shape_Functions(ii,3)*Elasticity_Matrix(2,1)*Shape_Functions(jj,2)';
Local_Stiffness_Matrix(i,j)=Local_Stiffness_Matrix(i,j)+Shape_Functions(ii,3)*Elasticity_Matrix(2,3)*Shape_Functions(jj,3)';
Local_Stiffness_Matrix(i,j)=Local_Stiffness_Matrix(i,j)+Shape_Functions(ii,2)*Elasticity_Matrix(3,1)*Shape_Functions(jj,2)';
Local_Stiffness_Matrix(i,j)=Local_Stiffness_Matrix(i,j)+Shape_Functions(ii,2)*Elasticity_Matrix(3,3)*Shape_Functions(jj,3)';
j=2*jj;
Local_Stiffness_Matrix(i,j)=Local_Stiffness_Matrix(i,j)+Shape_Functions(ii,3)*Elasticity_Matrix(2,2)*Shape_Functions(jj,3)';
Local_Stiffness_Matrix(i,j)=Local_Stiffness_Matrix(i,j)+Shape_Functions(ii,3)*Elasticity_Matrix(2,3)*Shape_Functions(jj,2)';
Local_Stiffness_Matrix(i,j)=Local_Stiffness_Matrix(i,j)+Shape_Functions(ii,2)*Elasticity_Matrix(3,2)*Shape_Functions(jj,3)';
Local_Stiffness_Matrix(i,j)=Local_Stiffness_Matrix(i,j)+Shape_Functions(ii,2)*Elasticity_Matrix(3,3)*Shape_Functions(jj,2)';
%
Local_Stiffness_Matrix=Local_Stiffness_Matrix*Weight*Jacobian;