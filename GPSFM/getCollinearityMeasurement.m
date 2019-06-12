function [ output_args ] = getCollinearityMeasurement( e12,e13,e21,e23,e31,e32,center )
%GETCOLLINEARITYMEASUREMENT Summary of this function goes here
%   Detailed explanation goes here

norm1=norm((e12+e13)/2-center);
m1=norm(e12/norm1-e13/norm1);


norm2=norm((e21+e23)/2-center);
m2=norm(e21/norm2-e23/norm2);

norm3=norm((e31+e32)/2-center);
m3=norm(e31/norm3-e32/norm3);

output_args=min([m1,m2,m3]);
end

