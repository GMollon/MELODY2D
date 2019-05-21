function Sorted=SortSegments(Unsorted)
%NB : only coded for a body with one hole, otherwise an additional loop is needed
Sorted=cell(1,1);

Sorted{1,1}=Unsorted(1,:)';
Last=Unsorted(1,2);
FlagOver=0;
while FlagOver==0
    FlagOver=1;
    for i=2:size(Unsorted,1)
        if Unsorted(i,1)==Last
            Last=Unsorted(i,2);
            if Last~=Sorted{1,1}(1)
                Sorted{1,1}=[Sorted{1,1};Last];
                FlagOver=0;
            end
        end
    end
end

for i=2:size(Unsorted,1);
    if isempty(find(Sorted{1,1}==Unsorted(i,1)))
        break
    end
end

Sorted{2,1}=Unsorted(i,:)';
Last=Unsorted(i,2);
FlagOver=0;
while FlagOver==0
    FlagOver=1;
    for i=2:size(Unsorted,1)
        if Unsorted(i,1)==Last
            Last=Unsorted(i,2);
            if Last~=Sorted{2,1}(1)
                Sorted{2,1}=[Sorted{2,1};Last];
                FlagOver=0;
            end
        end
    end
end
