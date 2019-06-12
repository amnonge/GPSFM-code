function [ FN ] = normalizeForbineusNorm( FN )
%NORMALIZRFORBINEUSNORM Summary of this function goes here
%   Detailed explanation goes here
nodesNum=size(FN,1)/3;
for ii=1:(nodesNum-1)
    for jj=ii+1:nodesNum
        FN(jj*3-2:jj*3,ii*3-2:ii*3)=FN(jj*3-2:jj*3,ii*3-2:ii*3)/norm(FN(jj*3-2:jj*3,ii*3-2:ii*3),'fro');
        FN(ii*3-2:ii*3,jj*3-2:jj*3)=FN(ii*3-2:ii*3,jj*3-2:jj*3)/norm(FN(ii*3-2:ii*3,jj*3-2:jj*3),'fro');
    end
end
end

