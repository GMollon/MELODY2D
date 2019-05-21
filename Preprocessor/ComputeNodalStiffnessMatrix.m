function Local_Stiffness_Matrix=ComputeNodalStiffnessMatrix(WeightX,WeightY,Surface,Elasticity_Matrix)

Nfnodes=size(WeightX,1);
Local_Stiffness_Matrix=zeros(2*Nfnodes,2*Nfnodes);
for ii=1:Nfnodes
    i=2*ii-1;
    for jj=1:Nfnodes
        j=2*jj-1;
        Local_Stiffness_Matrix(i,j)=Local_Stiffness_Matrix(i,j)+WeightX(ii,1)*Elasticity_Matrix(1,1)*WeightX(jj,1);
        Local_Stiffness_Matrix(i,j)=Local_Stiffness_Matrix(i,j)+WeightX(ii,1)*Elasticity_Matrix(1,3)*WeightY(jj,1);
        Local_Stiffness_Matrix(i,j)=Local_Stiffness_Matrix(i,j)+WeightY(ii,1)*Elasticity_Matrix(3,1)*WeightX(jj,1);
        Local_Stiffness_Matrix(i,j)=Local_Stiffness_Matrix(i,j)+WeightY(ii,1)*Elasticity_Matrix(3,3)*WeightY(jj,1);
    end
    for jj=1:Nfnodes
        j=2*jj;
        Local_Stiffness_Matrix(i,j)=Local_Stiffness_Matrix(i,j)+WeightX(ii,1)*Elasticity_Matrix(1,2)*WeightY(jj,1);
        Local_Stiffness_Matrix(i,j)=Local_Stiffness_Matrix(i,j)+WeightX(ii,1)*Elasticity_Matrix(1,3)*WeightX(jj,1);
        Local_Stiffness_Matrix(i,j)=Local_Stiffness_Matrix(i,j)+WeightY(ii,1)*Elasticity_Matrix(3,2)*WeightY(jj,1);
        Local_Stiffness_Matrix(i,j)=Local_Stiffness_Matrix(i,j)+WeightY(ii,1)*Elasticity_Matrix(3,3)*WeightX(jj,1);
    end
end
for ii=1:Nfnodes
    i=2*ii;
    for jj=1:Nfnodes
        j=2*jj-1;
        Local_Stiffness_Matrix(i,j)=Local_Stiffness_Matrix(i,j)+WeightY(ii,1)*Elasticity_Matrix(2,1)*WeightX(jj,1);
        Local_Stiffness_Matrix(i,j)=Local_Stiffness_Matrix(i,j)+WeightY(ii,1)*Elasticity_Matrix(2,3)*WeightY(jj,1);
        Local_Stiffness_Matrix(i,j)=Local_Stiffness_Matrix(i,j)+WeightX(ii,1)*Elasticity_Matrix(3,1)*WeightX(jj,1);
        Local_Stiffness_Matrix(i,j)=Local_Stiffness_Matrix(i,j)+WeightX(ii,1)*Elasticity_Matrix(3,3)*WeightY(jj,1);
    end
    for jj=1:Nfnodes
        j=2*jj;
        Local_Stiffness_Matrix(i,j)=Local_Stiffness_Matrix(i,j)+WeightY(ii,1)*Elasticity_Matrix(2,2)*WeightY(jj,1);
        Local_Stiffness_Matrix(i,j)=Local_Stiffness_Matrix(i,j)+WeightY(ii,1)*Elasticity_Matrix(2,3)*WeightX(jj,1);
        Local_Stiffness_Matrix(i,j)=Local_Stiffness_Matrix(i,j)+WeightX(ii,1)*Elasticity_Matrix(3,2)*WeightY(jj,1);
        Local_Stiffness_Matrix(i,j)=Local_Stiffness_Matrix(i,j)+WeightX(ii,1)*Elasticity_Matrix(3,3)*WeightX(jj,1);
    end
end
Local_Stiffness_Matrix=Local_Stiffness_Matrix/Surface;