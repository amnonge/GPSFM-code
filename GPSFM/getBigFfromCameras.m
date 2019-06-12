function [Fnew,Anew]=getBigFfromCameras(Ps,finalTriplets,nodesNum)
Us=[];
Vs=[];
for i=1:nodesNum
    [is,js]=find(finalTriplets==i);
    curP=Ps{is(1)}{js(1)};
     V4f=inv(curP(1:3,1:3))';
T4f=null(curP);
T4f=T4f(1:3)/T4f(4);
T4f=getCrossM(T4f);
U4f=V4f*T4f;
Us=[Us;U4f];
Vs=[Vs;V4f];
end

for i=1:size(Vs,1)/3
    normalization=norm(Vs(3*i-2:3*i),'fro');
    Vs(3*i-2:3*i,:)=Vs(3*i-2:3*i,:)/normalization;
    Us(3*i-2:3*i,:)=Us(3*i-2:3*i,:)/normalization;
end
Fnew=Vs*Us.'+Us*Vs.';
Anew=Vs*Us.';