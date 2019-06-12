function [ res ] = getCrossM( v )
%GETCROSSM Summary of this function goes here
%   Detailed explanation goes here
t1=cross(v,[1;0;0]);
t2=cross(v,[0;1;0]);
t3=cross(v,[0;0;1]);

res=[t1 t2 t3];
end

