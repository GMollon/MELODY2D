function Listen=SqueezeList(Liste)
Liste=sortrows(Liste,1);
Listen=zeros(size(Liste,1),size(Liste,2));
Listen(1,:)=Liste(1,:);
num=1;
for i=2:size(Liste,1)
    if Liste(i,1)==Liste(i-1,1)
        continue
    else
        num=num+1;
        Listen(num,:)=Liste(i,:);
    end
end
Listen=Listen(1:num,:);
    