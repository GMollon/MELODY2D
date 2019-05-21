function List_DOF=CreateListDOF(List_Domain_Fnodes)
List_DOF=zeros(size(List_Domain_Fnodes,1)*2,1);
for i=1:size(List_Domain_Fnodes,1)
    List_DOF(2*i-1,1)=2*List_Domain_Fnodes(i,1)-1;
    List_DOF(2*i,1)=2*List_Domain_Fnodes(i,1);
end