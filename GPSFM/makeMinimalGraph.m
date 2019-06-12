function [ firstGroup,G]=makeMinimalGraph(selectedTriplets,firstGroup,G,nodesNum)


remainedNodes=zeros(nodesNum*2,length(firstGroup));
remainedNodesGlobalIndex=zeros(nodesNum*2,1);
indexremainedNodes=1;
maskTriplets=ones(size(selectedTriplets,1),1);
[~,indsss]=sort(sum(selectedTriplets(firstGroup,4:6),2));
indsssInGraph=indsss;
for i=1:length(firstGroup)   
    maskTripletstemp=maskTriplets;
    maskTripletstemp(indsss(i))=0;
    reduced=selectedTriplets(logical(maskTripletstemp),1:3);
    H= rmnode(G,(indsssInGraph(i)));
    indss=reduced(:);  
    if isempty(setdiff(1:nodesNum,unique(indss)))==0
        continue;
    end
    ccCurent=conncomp(H);
    if length(unique(ccCurent))~=1
        remainedNodes(indexremainedNodes,logical(maskTripletstemp))=ccCurent;
        remainedNodesGlobalIndex(indexremainedNodes)=indsss(i);
        indexremainedNodes=indexremainedNodes+1;
    else
        G=H;
        indsssInGraph(indsssInGraph>indsssInGraph(i))=indsssInGraph(indsssInGraph>indsssInGraph(i))-1;
        foundSomething=true;
        maskTriplets=maskTripletstemp;
        remainedNodes=remainedNodes.*repmat(maskTripletstemp',size(remainedNodes,1),1); 
        for ii=1:indexremainedNodes-1
            if length(unique(remainedNodes(ii,:)))==2  
                maskTripletstemp=maskTriplets;
                curInd=remainedNodesGlobalIndex(ii);
                maskTripletstemp(curInd)=0;
                reduced=selectedTriplets(logical(maskTripletstemp),1:3);
                ttt=find(find(maskTriplets)==curInd);
                indss=reduced(:);
                if isempty((setdiff(1:nodesNum,unique(indss))))==0
                    continue;
                end
                H= rmnode(G,ttt);
                ccCurent=conncomp(H);
                if length(unique(ccCurent))==1
                    G=H;
                    indsssInGraph(indsssInGraph>ttt)=indsssInGraph(indsssInGraph>ttt)-1;
                    remainedNodes(ii,:)=0;
                    maskTriplets=maskTripletstemp;
                    remainedNodes=remainedNodes.*repmat(maskTripletstemp',size(remainedNodes,1),1);
                end            
            end 
        end
    end
end


 
 for ii=1:indexremainedNodes-1
            if length(unique(remainedNodes(ii,:)))==2  
                maskTripletstemp=maskTriplets;
                curInd=remainedNodesGlobalIndex(ii);
                maskTripletstemp(curInd)=0;
                reduced=selectedTriplets(logical(maskTripletstemp),1:3);
                ttt=find(find(maskTriplets)==curInd);
                indss=reduced(:);
                if isempty((setdiff(1:nodesNum,unique(indss))))==0
                    continue;
                end
                H= rmnode(G,ttt);
                ccCurent=conncomp(H);
                if length(unique(ccCurent))==1
                    G=H;
                    indsssInGraph(indsssInGraph>ttt)=indsssInGraph(indsssInGraph>ttt)-1;
                    remainedNodes(ii,:)=0;
                    maskTriplets=maskTripletstemp;
                    remainedNodes=remainedNodes.*repmat(maskTripletstemp',size(remainedNodes,1),1);
                end            
            end 
        end


firstGroup=find(maskTriplets);
end

function  [G,maskTriplets,indsssInGraph,foundSomething]=removePossibleLeaves(maskTriplets,G,indsssInGraph,indsss,i,selectedTriplets,nodesNum)
originalIndicesInGraph=find(maskTriplets);
AA=adjacency(G);

leavesGraphIndices=find(sum(AA)==1);

leavesOriginalIndices=originalIndicesInGraph(leavesGraphIndices);

%We need to sort them according to indsss
[~,inds]=sort(sum(selectedTriplets(leavesOriginalIndices,4:6),2));

leavesOriginalIndices=leavesOriginalIndices(inds);
leavesGraphIndices=leavesGraphIndices(inds);
foundSomething=false;
for jj=1:length(leavesOriginalIndices)
    j=leavesOriginalIndices(jj);
    
    if find(indsss==j)>i
        %we are not there yets
        
        dd=0;
        break;
    end
    maskTripletsTemp=maskTriplets;
    maskTripletsTemp(j)=0;
    if length(unique(selectedTriplets(logical(maskTripletsTemp),1:3)))==nodesNum
%         j
        foundSomething=true;
        maskTriplets=maskTripletsTemp;
        G=rmnode(G,leavesGraphIndices(jj));
          indsssInGraph(indsssInGraph>leavesGraphIndices(jj))=indsssInGraph(indsssInGraph>leavesGraphIndices(jj))-1;
          leavesGraphIndices(leavesGraphIndices>leavesGraphIndices(jj))=leavesGraphIndices(leavesGraphIndices>leavesGraphIndices(jj))-1;
    end
    
end



end

