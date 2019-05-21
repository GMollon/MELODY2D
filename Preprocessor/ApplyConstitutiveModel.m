function [PK,Cauchy,Elasticity_Matrix,PK3,Cauchy3]=ApplyConstitutiveModel(F,Material_Type,Material_Parameters)
if strcmp(Material_Type,'ElasticLinear')
    GL=0.5*(F'*F-[1,0;0,1]);
    GL=[GL(1,1);GL(2,2);2*GL(1,2)];
    J=det(F);
    Young_Modulus=Material_Parameters(1,4);
    Poisson_Coef=Material_Parameters(1,5);
    if Poisson_Coef==0.5
        Poisson_Coef=0.49;
    end
    Elasticity_Matrix=zeros(3,3);
    Elasticity_Matrix(1,1)=1-Poisson_Coef;
    Elasticity_Matrix(1,2)=Poisson_Coef;
    Elasticity_Matrix(2,1)=Poisson_Coef;
    Elasticity_Matrix(2,2)=1-Poisson_Coef;
    Elasticity_Matrix(3,3)=0.5*(1-2*Poisson_Coef);
    Elasticity_Matrix=Elasticity_Matrix*Young_Modulus/((1+Poisson_Coef)*(1-2*Poisson_Coef));
    PK=Elasticity_Matrix*GL;
    PK3=Poisson_Coef*(GL(1)+GL(2))*Young_Modulus/((1+Poisson_Coef)*(1-2*Poisson_Coef));
elseif strcmp(Material_Type,'LargeElasticLinear')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% NB : doit être linéarisé !!! %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Young_Modulus=Material_Parameters(1,4);
    Poisson_Coef=Material_Parameters(1,5);
    if Poisson_Coef==0.5
        Poisson_Coef=0.49;
    end
    %Elasticity_Matrix=zeros(3,3);
    %Elasticity_Matrix(1,1)=1-Poisson_Coef;
    %Elasticity_Matrix(1,2)=Poisson_Coef;
    %Elasticity_Matrix(2,1)=Poisson_Coef;
    %Elasticity_Matrix(2,2)=1-Poisson_Coef;
    %Elasticity_Matrix(3,3)=0.5*(1-2*Poisson_Coef);
    %Elasticity_Matrix=Elasticity_Matrix*Young_Modulus/((1+Poisson_Coef)*(1-2*Poisson_Coef));
    %
    %A=Elasticity_Tensor(Young_Modulus,Poisson_Coef);
    %F=[F,[0;0];[0,0,1]];
    %J=det(F);
    %GL=0.5*(F'*F-diag(ones(3,1)));
    %invF=inv(F);
    %B=invF'*GL*invF;
    %B=Double_dot(A,B);
    %PK=J*invF*B*invF';
    %
    %PK3=PK(3,3);
    %PK=[PK(1,1);PK(2,2);PK(1,2)];
    %F=F(1:2,1:2);
    %
    J=det(F);
    GL=0.5*(F'*F-[1,0;0,1]);
    GL=[GL(1,1);GL(2,2);2*GL(1,2)];
    %A=Elasticity_Tensor(Young_Modulus,Poisson_Coef);
    Lambda=Young_Modulus*Poisson_Coef/((1+Poisson_Coef)*(1-2*Poisson_Coef));
    Mu=Young_Modulus/(2*(1+Poisson_Coef));
    %
    F=[F,[0;0];[0,0,1]];
    G=inv(F);
    Elasticity_Matrix=zeros(3,3);
    for i=1:3
        for j=1:3
            for k=1:3
                for l=1:3
                    I=(Kronecker(i,k)*Kronecker(j,l)+Kronecker(j,k)*Kronecker(i,l))/2;
                    J=Kronecker(i,j)*Kronecker(k,l)/3;
                    A=3*Lambda*J+2*Mu*I;
                    Elasticity_Matrix(1,1)=Elasticity_Matrix(1,1)+G(1,i)*G(1,j)*A*G(1,k)*G(1,l);
                    Elasticity_Matrix(1,2)=Elasticity_Matrix(1,2)+G(1,i)*G(1,j)*A*G(2,k)*G(2,l);
                    Elasticity_Matrix(1,3)=Elasticity_Matrix(1,3)+G(1,i)*G(1,j)*A*G(1,k)*G(2,l);
                    Elasticity_Matrix(2,1)=Elasticity_Matrix(2,1)+G(2,i)*G(2,j)*A*G(1,k)*G(1,l);
                    Elasticity_Matrix(2,2)=Elasticity_Matrix(2,2)+G(2,i)*G(2,j)*A*G(2,k)*G(2,l);
                    Elasticity_Matrix(2,3)=Elasticity_Matrix(2,3)+G(2,i)*G(2,j)*A*G(1,k)*G(2,l);
                    Elasticity_Matrix(3,1)=Elasticity_Matrix(3,1)+G(1,i)*G(2,j)*A*G(1,k)*G(1,l);
                    Elasticity_Matrix(3,2)=Elasticity_Matrix(3,2)+G(1,i)*G(2,j)*A*G(2,k)*G(2,l);
                    Elasticity_Matrix(3,3)=Elasticity_Matrix(3,3)+G(1,i)*G(2,j)*A*G(1,k)*G(2,l);
                end
            end
        end
    end
    Elasticity_Matrix=J*Elasticity_Matrix;
    PK=Elasticity_Matrix*GL;
    
    PK3=Poisson_Coef*(GL(1)+GL(2))*Young_Modulus/((1+Poisson_Coef)*(1-2*Poisson_Coef));
    F=F(1:2,1:2);
elseif strcmp(Material_Type,'NeoHookean')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% NB : doit être linéarisé !!! %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    GL=0.5*(F'*F-[1,0;0,1]);
    GL=[GL(1,1);GL(2,2);2*GL(1,2)];
    J=det(F);
    Young_Modulus=Material_Parameters(1,4);
    Poisson_Coef=Material_Parameters(1,5);
    if Poisson_Coef==0.5
        Poisson_Coef=0.49;
    end
    Lambda=Young_Modulus*Poisson_Coef/((1+Poisson_Coef)*(1-2*Poisson_Coef));
    Mu=Young_Modulus/(2*(1+Poisson_Coef));
    LambdaP=Lambda/J;
    MuP=(Mu-Lambda*log(J))/J;
    Young_ModulusP=MuP*(3*LambdaP+2*MuP)/(LambdaP+MuP);
    Poisson_CoefP=LambdaP/(2*(LambdaP+MuP));
    Elasticity_Matrix=zeros(3,3);
    Elasticity_Matrix(1,1)=1-Poisson_CoefP;
    Elasticity_Matrix(1,2)=Poisson_CoefP;
    Elasticity_Matrix(2,1)=Poisson_CoefP;
    Elasticity_Matrix(2,2)=1-Poisson_CoefP;
    Elasticity_Matrix(3,3)=0.5*(1-2*Poisson_CoefP);
    Elasticity_Matrix=Elasticity_Matrix*Young_ModulusP/((1+Poisson_CoefP)*(1-2*Poisson_CoefP));
    PK=Elasticity_Matrix*GL;
    PK3=Poisson_CoefP*(GL(1)+GL(2))*Young_ModulusP/((1+Poisson_CoefP)*(1-2*Poisson_CoefP));
end
Cauchy=F*[PK(1),PK(3);PK(3),PK(2)]*F'/J;
Cauchy=[Cauchy(1,1);Cauchy(2,2);Cauchy(1,2)];
Cauchy3=PK3;