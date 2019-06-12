function [triplpets]=extractTripletsFromViewingGraph(tripletGraph)
L=laplacian(tripletGraph);

triplpets=zeros(size(tripletGraph.Nodes,1)^2,3);
counter = 0;
for i=1:size(L,1)
    rel=find(L(i,i+1:end)<0);
    rel=rel+i;
    if length(rel)<2
        continue;
    end
    possiblePairs=nchoosek(rel,2);
    for j=1:size(possiblePairs,1)
        if L(possiblePairs(j,1),possiblePairs(j,2))

            counter = counter +1;
            triplpets(counter,:) =[i possiblePairs(j,1) possiblePairs(j,2)];
            
        end
    end
end
triplpets = triplpets(1:counter,:);
end