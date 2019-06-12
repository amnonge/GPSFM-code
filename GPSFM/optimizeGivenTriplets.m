function [Xs,Y] = optimizeGivenTriplets( FN,Cf,show)

if nargin<3
    show=true;
end
nn=size(FN,1);
n=size(Cf,1);

GAMMAS=cell(n,1);
FNS=cell(n,1);
YS=cell(n,1);
Xs=cell(n,1);
for i=1:n
    GAMMAS{i}=zeros(9,9);
    FNS{i}=getTrippleInds(FN, Cf(i,:) );
    
    YS{i}=FNS{i};
end



errors=zeros(n,1000);
sings=zeros(n,1000);
sings5=zeros(n,1000);
tic

wss=ones(n,1);
count=1;


for i=1:25
    
    for j=1:n
        
        Xs{j}= proxSVD( YS{j}-GAMMAS{j}, FNS{j});
        
    end
    
    
    
    Y=solve_Y4(GAMMAS,Xs,nn,Cf,ones(size(wss)));

    for j=1:n
        YS{j}=getTrippleInds(Y, Cf(j,:) );
        GAMMAS{j}=GAMMAS{j}+Xs{j}-YS{j};
        errors(j,count)=norm(FNS{j}-YS{j},'fro');
        ttt=svd(YS{j});
        sings(j,count)=ttt(7);
        
        sings5(j,i)=ttt(6);
        
    end
    count=count+1;
end

if show
    figure('DefaultAxesFontSize',15), plot(log10(errors'),'LineWidth',2),ylim([-3.5 0]),xlabel('iterations'),ylabel('log_{10} ||F_i^-F_i||'),legend('Tr_1','Tr_2','Tr_3','Tr_4','Tr_5','Tr_6','Tr_7','Tr_8','Tr_9')
    figure('DefaultAxesFontSize',15), plot(log10(sings'),'LineWidth',2),xlabel('iterations'),ylabel('log_{10} (\sigma_i(7))'),legend('Tr_1','Tr_2','Tr_3','Tr_4','Tr_5','Tr_6','Tr_7','Tr_8','Tr_9')
end
Xs=YS;

objectives=zeros(nn/3,nn/3);
for i=1:nn/3-1
    for j=i:nn/3
        objectives(i,j)=norm(Y(3*i-2:3*i,3*j-2:3*j)-FN(3*i-2:3*i,3*j-2:3*j),'fro');
    end
end


end


function Y=solve_Y4(GAMMAS,Xs,nn,Cf,wss)
Y=zeros(nn);
W=zeros(nn);

for i=1:length(GAMMAS)
    [~,tr]=getTrippleInds(zeros(nn), Cf(i,:) );
    Y(tr,tr)=Y(tr,tr)+Xs{i};
    W(tr,tr)= W(tr,tr)+ones(9)*wss(i);
end
Y=Y./W;
Y(W==0)=0;
Y=(Y+Y')/2;
for i=1:nn/3
    Y(3*i-2:3*i,3*i-2:3*i)=0;
end


end



function F=proxSVD( Y1, Y2 )
res=10*Y1/11+Y2/11;

[u,d,v]=svd(res);
temp=diag(d);
temp(7:end)=0;
d=diag(temp);
F=u*d*v';
end



function [ma,tr]=getTrippleInds(BigMatrix, tripleInds )
tr=[tripleInds(1)*3-2:tripleInds(1)*3 tripleInds(2)*3-2:tripleInds(2)*3 tripleInds(3)*3-2:tripleInds(3)*3];
ma=BigMatrix(tr,tr);
end
