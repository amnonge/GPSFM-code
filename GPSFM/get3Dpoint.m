function [Points3D, visible,MRep] = get3Dpoint(Ps,Mignored,Mfull)
W=Mignored~=0;
Points3D=zeros(4,size(Mignored,2));

for i=1:size(Mignored,2)
    xsss=[];
    Pscc={};
    eqs=[];
    for j=1:size(Mignored,1)/2
        if(W((j-1)*2+1,i))
            x1=Mignored((j-1)*2+1:(j-1)*2+2,i);
            curEqs=[x1(1)*Ps{j}(3,:)-Ps{j}(1,:);x1(2)*Ps{j}(3,:)-Ps{j}(2,:)];
            eqs=[eqs;curEqs];
            xsss=[xsss x1];
            Pscc{length(Pscc)+1}=Ps{j};
        end
    end
    [u,d,v]=svd(eqs);
    X=v(:,end);
    Points3D(:,i)=X(1:4)/X(4);

    X = vgg_X_from_xP_nonlin(xsss,Pscc,repmat([576   720 ]',length(Pscc)));

    Points3D(:,i)=X(1:4)/X(4);
end
    visible = Mfull~=0;
    visible = visible(1:2:end,:);
    
    MRep = [];
for i = 1:2:size(Mfull,1)
    MRep = [MRep;[Mfull(i:i+1,:);ones(1,size(Mfull,2))]];
end
    
end