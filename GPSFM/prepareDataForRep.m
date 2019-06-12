function [ Ps ,visible,MRep] = prepareDataForRep(PsList,Points3Dold,M)

Ps = zeros(3*length(PsList),4);
for i = 1:length(PsList)
    Ps(3*i-2:3*i,1:4) = PsList{i};
end

visible = M~=0;
visible = visible(1:2:end,:);

MRep = [];
for i = 1:2:size(M,1)
    MRep = [MRep;[M(i:i+1,:);ones(1,size(M,2))]];
end

end