function [GaussPoints,GaussWeights,Jacobian]=TriangleGaussPoints(X1,Y1,X2,Y2,X3,Y3,N_Gauss_Triangle)
%
%Provides the Gauss points and the corresponding weights for a triangular
%cell of the backgroud mesh, as well as the corresponding jacobian
%
%Input variables :
%    Xi : x-coordinate of the background triangle
%    Yi : x-coordinate of the background triangle
%    N_Gauss_Triangle : number of Gauss points in the triangle (3, 6, or 7)
%
%Output variables :
%    GaussPoints : coordinates of the gauss points of the cell
%    GaussWeights : corresponding weights
%    Jacobian : jacobian of the affine transform
%
GaussPW=GaussCoefficientsTriangle(N_Gauss_Triangle);
GaussWeights=GaussPW(:,3);
GaussPoints=zeros(N_Gauss_Triangle,2);
for i=1:N_Gauss_Triangle
    a1=GaussPW(i,1);
    a2=GaussPW(i,2);
    a3=1-a1-a2;
    x=a3*X1+a1*X2+a2*X3;
    y=a3*Y1+a1*Y2+a2*Y3;
    GaussPoints(i,:)=[x,y];
end
Jacobian=abs(X1*(Y2-Y3)+X2*(Y3-Y1)+X3*(Y1-Y2))/2;