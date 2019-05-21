function List_Domain_Fnodes=SupportDomain(Poi_Position,Fnodes_Positions,Domains_Dimensions,Support_Neighbours,Connections,Domain)
% FlagOut=0;
% FlagNext=0;
% while FlagOut==0
%     NewDomain=[];
%     for j=Domain'
%         NewDomain=[NewDomain;Connections{j,1}];
%     end
%     Domain=unique(NewDomain);
%     if FlagNext==1
%         FlagOut=1;
%     elseif size(Domain,1)>Support_Neighbours
%         FlagNext=1;
%     end
% end
Domain=[1:size(Fnodes_Positions,1)]';
Distances=[(Fnodes_Positions(Domain,1)-Poi_Position(1)).^2+(Fnodes_Positions(Domain,2)-Poi_Position(2)).^2-Domains_Dimensions(Domain).^2,Domain];
List_Domain_Fnodes=Distances(find(Distances(:,1)<0),2);





% List_Domain_Fnodes=[];
% for k=1:size(Fnodes_Positions)
%     if sqrt((Fnodes_Positions(k,1)-Poi_Position(1))^2+(Fnodes_Positions(k,2)-Poi_Position(2))^2)<Domains_Dimensions(k,1)
%         List_Domain_Fnodes=cat(1,List_Domain_Fnodes,k);
%     end
% end