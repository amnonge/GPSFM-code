function [ Fc,A ] = projectF( F )
[V,D]=eig(F);

dd=diag(D);
[~,inds]=sort(abs(dd),'descend');
d=dd(inds(1:6));
pd=d(d>0);
nd=d(d<0);

DD=diag(d);
VV=V(:,inds(1:6));
VP=VV(:,d>0);
VN=VV(:,d<0);
DDN=diag(nd);
DDP=diag(pd);


F-VV*DD*VV';
Fc=VP*DDP*VP'+VN*DDN*VN';
X=VP*DDP.^0.5;
Y=VN*abs(DDN).^0.5;
Fc2=X*X'-Y*Y';

Fc3= [[X+Y]*[X-Y]'+[X-Y]*[X+Y]']/2;

A=[X+Y]*[X-Y]'/2;

norm(F-Fc);
end

