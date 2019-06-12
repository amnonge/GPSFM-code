function [ pointsNorm ] = pnorm( points )
%PNORM Summary of this function goes here
%   Detailed explanation goes here
pointsNorm=points./repmat(points(end,:),size(points,1),1);

end

