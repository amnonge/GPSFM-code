function [ error ] = TripletError( F )
%SWEENYTRIPLETERROR Summary of this function goes here
%   Detailed explanation goes here


[Xs] = optimizeGivenTriplets(F,[1 2 3],false);

error=norm(Xs{1}-F,'fro');

end

