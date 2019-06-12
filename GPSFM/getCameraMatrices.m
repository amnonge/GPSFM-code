function [ Us,Vs,Ts,Ps,Ffound ] = getCameraMatrices( F )
%GETCAMERAMATRICES Summary of this function goes here
%   Detailed explanation goes here

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



P1 = perms(1:3);

VPg=VP;
DDPg=DDP;
VNg=VN;
DDNg=DDN;
errorsss=zeros(6,6);
FFFs=cell(6,6);
for iii=1:size(P1,1)
    for jjj=1:size(P1,1)
        cur111=P1(iii,:);
        cur222=P1(jjj,:);
        
        
        VP=VPg(:,cur111);
        DDP=diag(DDPg);
        DDP=DDP(cur111);
        DDP=diag(DDP);
        
        X=VP*DDP.^0.5;
        VN=VNg(:,cur222);
        DDN=diag(DDNg);
        DDN=DDN(cur222);
        DDN=diag(DDN);
        
        Y=VN*abs(DDN).^0.5;
        
        if cond(X(1:3,1:3)+Y(1:3,1:3))<cond(X(1:3,1:3)-Y(1:3,1:3))
            Vfound=[X+Y];
            Ufound=[X-Y]/2;
        else
            Vfound=[X-Y];
            Ufound=[X+Y]/2;
        end
        Tfound = zeros(size(Vfound));
        Us = zeros(size(Ufound));
        
        for i = 1:3:size(Vfound,1)
            tempT = inv(Vfound(i:i+2,:))*Ufound(i:i+2,:);
            tempT = 0.5*(tempT-tempT');
            tempTvec = [tempT(3,2) tempT(1,3) tempT(2,1)]';
            Tfound(i:i+2,:) = getCrossM(tempTvec);
            Us(i:i+2,:) = Vfound(i:i+2,:)*Tfound(i:i+2,:);
        end
        
        Vs = Vfound;
        Ts = Tfound;
        Ps=cell(size(Vs,1)/3,1);
        for i=1:(size(Vs,1)/3)
            t=zeros(3,1);curT=Ts(i*3-2:i*3,:);curV=Vs(i*3-2:i*3,:);
            t(1)=curT(3,2);
            t(2)=curT(1,3);
            
            t(3)=curT(2,1);
            Ps{i}=inv(curV)'*[eye(3) -t];
        end
        
        Ffound = Us*Vs'+Vs*Us';
        errorsss(iii,jjj)=norm(Ffound-F,'fro');
        FFFs{iii,jjj}=Ffound;
    end
end
[~,ind]=min(errorsss(:));
[iii,jjj]=ind2sub([6 6] ,ind);
cur111=P1(iii,:);
cur222=P1(jjj,:);


VP=VPg(:,cur111);
DDP=diag(DDPg);
DDP=DDP(cur111);
DDP=diag(DDP);

X=VP*DDP.^0.5;
VN=VNg(:,cur222);
DDN=diag(DDNg);
DDN=DDN(cur222);
DDN=diag(DDN);


Y=VN*abs(DDN).^0.5;


if cond(X(1:3,1:3)+Y(1:3,1:3))<cond(X(1:3,1:3)-Y(1:3,1:3))
    Vfound=[X+Y];
    Ufound=[X-Y]/2;
else
    Vfound=[X-Y];
    Ufound=[X+Y]/2;
end
Tfound = zeros(size(Vfound));
Us = zeros(size(Ufound));

for i = 1:3:size(Vfound,1)
    tempT = inv(Vfound(i:i+2,:))*Ufound(i:i+2,:);
    tempT = 0.5*(tempT-tempT');
    tempTvec = [tempT(3,2) tempT(1,3) tempT(2,1)]';
    Tfound(i:i+2,:) = getCrossM(tempTvec);
    Us(i:i+2,:) = Vfound(i:i+2,:)*Tfound(i:i+2,:);
end

Vs = Vfound;
Ts = Tfound;
Ps=cell(size(Vs,1)/3,1);
for i=1:(size(Vs,1)/3)
    t=zeros(3,1);curT=Ts(i*3-2:i*3,:);curV=Vs(i*3-2:i*3,:);
    t(1)=curT(3,2);
    t(2)=curT(1,3);
    
    t(3)=curT(2,1);
    Ps{i}=inv(curV)'*[eye(3) -t];
end

Ffound = Us*Vs'+Vs*Us';

