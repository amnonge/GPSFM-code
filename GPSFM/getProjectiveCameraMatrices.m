function Ps=getProjectiveCameraMatrices(Xs,v,finalTriplets)

%GETPROJECTIVECAMERAMATRICES Summary of this function goes here
%   Detailed explanation goes here


Ps={};
[ ~,~,~,Psc,~ ] = getCameraMatrices(  projectF( Xs{1} ) );%getCameraMatrices( projectF( Xs{1} ));
Ps{1}=Psc;



for i=2:length(Xs)
    [ ~,~,~,Psc,~ ] = getCameraMatrices( projectF( Xs{i} ) );
    Ps{i}=Psc;
end
for i=1:size(v,1)
    curEdge=v(i,:);
    [C,ia,ib]=intersect(finalTriplets(curEdge(1),1:3),finalTriplets(curEdge(2),1:3));
    [H3]=  findHomography({Ps{curEdge(2)}{ib(1)},Ps{curEdge(2)}{ib(2)}},{Ps{curEdge(1)}{ia(1)},Ps{curEdge(1)}{ia(2)}});
    H3=double(H3);
    H3 =  H3/norm(H3,'fro');
    svd(Ps{curEdge(2)}{1});
    Ps{curEdge(2)}{1}=Ps{curEdge(2)}{1}*(H3);
    Ps{curEdge(2)}{2}=Ps{curEdge(2)}{2}*(H3);
    Ps{curEdge(2)}{3}=Ps{curEdge(2)}{3}*(H3); 
end
end

