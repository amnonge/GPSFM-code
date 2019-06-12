function [Xs,Y] = optimizeTripletsIRLS( FN,Cf,show,doIRLS)
%OPTIMZEOLSSON Summary of this function goes here
%   Detailed explanation goes here

%Olson formulation is  min_X R_r0(X) + |X-X0|^2 + rho * |X-P|^2
%in our case x0 in FN, P comes from the ADMM

%so our formulation is 
%X_{t+1}=argmin R_r0(X)+|X-X0|^2+row|X-Y_{t}+GAMMA_{t}||^2
%Y_{t+1}=argmin row||X_{t+1}-Y+GAMMA_t||^2 s.t Y symmetric, and Y_ii=0
%GAMMA_{t+1}=GAMMA_{t}+X_{t+1}-Y_{t+1}


%GAMMA_{0}=zeros
%Y_{0}=FN
% [C]=nchoosek(1:size(FN,2)/3,3);
% Cf=[];
% maxDistCameras=3;
% 
% minimalInds=[];
% for i=1:size(C,1)
%     cur=C(i,:);
%     minC=min(cur);
%     maxC=max(cur);
%     
%     if (maxC-minC)<maxDistCameras
%        Cf=[Cf;cur];
%         if (maxC-minC)<3
%          minimalInds=[minimalInds;size(Cf,1)];
%          end
%     end
%    
% end
if nargin<3
    show=true;
end
% Cf=[Cf;18 19 20;19 20 21];
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


% 
errors=zeros(n,30000);
sings=zeros(n,30000);
sings5=zeros(n,30000);
% tic

   wss=ones(n,1);
   count=1;
   
   if doIRLS
       iterssIrlss=20;
   else
       iterssIrlss=1;
   end
for irls=1:iterssIrlss
% irls
for i=1:1000%1000
    
    for j=1:n

 Xs{j}= proxSVD( YS{j}-GAMMAS{j}, YS{j}-GAMMAS{j});
        
    end
    
   
    
    Y=solve_Y4(GAMMAS,Xs,nn,Cf,wss,FNS);
  
for j=1:n
     YS{j}=getTrippleInds(Y, Cf(j,:) );
     GAMMAS{j}=GAMMAS{j}+Xs{j}-YS{j};
     errors(j,count)=norm(FNS{j}-YS{j},'fro');
     ttt=svd(YS{j});
     sings(j,count)= ttt(7)/ttt(6);
     
     sings5(j,i)=ttt(6);
     
end
count=count+1;
end
for j=1:n
    wss(j)=norm(getTrippleInds(Y, Cf(j,:))-getTrippleInds(FN, Cf(j,:) ),'fro')^-1;
end
wss=wss/mean(wss);
%  wss=wss.^2;

end

% rrr=0;
if show
% figure('DefaultAxesFontSize',15), plot(log10(errors'),'LineWidth',2),ylim([-3.5 0]),xlabel('iterations'),ylabel('log_{10} ||F_i^-F_i||'),legend('Tr_1','Tr_2','Tr_3','Tr_4','Tr_5','Tr_6','Tr_7','Tr_8','Tr_9')
%  figure('DefaultAxesFontSize',15), plot(log10(sings'),'LineWidth',2),xlabel('iterations'),ylabel('log_{10} (\sigma_i(7))'),legend('Tr_1','Tr_2','Tr_3','Tr_4','Tr_5','Tr_6','Tr_7','Tr_8','Tr_9')
a = log10(mean(sings,1));
a = a(1:1000);
noiseVec = 1:length(a);
p  =regress(a',[1./noiseVec;noiseVec;ones(size(noiseVec))]');%polyfit(noiseVec,a,2);%regress(a',[1./noiseVec;noiseVec;ones(size(noiseVec))]'); %polyfit(noiseVec,a,2);regress(a',[1./noiseVec;noiseVec;ones(size(noiseVec))]'); %
newa = [];

for i = 1:length(noiseVec)
    newa(i) = p(1)/noiseVec(i)+p(2)*noiseVec(i)+p(3);%p(1)*noiseVec(i)^2+p(2)*noiseVec(i)+p(3);% p(1)/noiseVec(i)+p(2)*noiseVec(i)+p(3);%
end
newa(800:end) = newa(800);
figure('DefaultAxesFontSize',19), plot((  newa),'LineWidth',3),xlabel('Iterations'),ylabel('$\displaystyle log_{10}(\frac{\sigma_7}{\sigma_6})$','interpreter','latex')
% figure('DefaultAxesFontSize',19), plot((  a),'LineWidth',3),xlabel('Iterations'),ylabel('$\displaystyle log_{10}(\frac{\sigma_7}{\sigma_6})$','interpreter','latex')
end
Xs=YS;

objectives=0;%zeros(nn/3,nn/3);
% for i=1:nn/3-1
%        for j=i:nn/3
%         objectives(i,j)=norm(Y(3*i-2:3*i,3*j-2:3*j)-FN(3*i-2:3*i,3*j-2:3*j),'fro');
%        end
% end


end


function Y=solve_Y4(GAMMAS,Xs,nn,Cf,wss,FNS)
    Y=zeros(nn);
    W=zeros(nn);
    param = 0.001;
    for i=1:length(GAMMAS)
        [~,tr]=getTrippleInds(zeros(nn), Cf(i,:) );
        Y(tr,tr)=Y(tr,tr)+Xs{i}+GAMMAS{i}+param*wss(i)*FNS{i};
         W(tr,tr)= W(tr,tr)+(1+param*wss(i))*ones(9);
    end
    Y=Y./W;
    Y(W==0)=0;
    Y=(Y+Y')/2;
    for i=1:nn/3
        Y(3*i-2:3*i,3*i-2:3*i)=0;
    end
    
end

% function Y=solve_Y3(GAMMAS,XS,n,Cf)
%     [xx,yy]=meshgrid(4:6,1:3);
% 
% [xx2,yy2]=meshgrid(7:9,1:3);
% [xx3,yy3]=meshgrid(7:9,4:6);
% idsX=[xx(:);xx2(:);xx3(:)];
% idsy=[yy(:);yy2(:);yy3(:)];
% inds=sub2ind([9 9], idsy,idsX);
% 
% allInds=[];
% allValues=[];
% for j=1:length(GAMMAS)
%  
%          [~,inss]=getTrippleInds(zeros(n), Cf(j,:) );
%          indsCur=sub2ind([n n], inss(idsy),inss(idsX) )';
%             curMat=GAMMAS{j}+XS{j};
%             curMat_t=curMat';
%             allInds=[allInds;indsCur;indsCur];
%             allValues=[allValues;curMat(inds);curMat_t(inds)];
%             
%             
%             
%     
%     
% end
% AA=sparse(1:size(allValues,1),allInds,ones(size(allValues,1),1),size(allValues,1),n*n);
%  Y=reshape(AA\allValues,n,n);
%  Y=Y+Y';
% 
% end
% 
% function Y=solve_Y(GAMMAS,Xs,n,Cf)
% % n=18;
% 
% for j=1:length(GAMMAS)
%     if j==1
%         [~,inss]=getTrippleInds(zeros(n), Cf(j,:) );
%         inss=mat2str(inss);
%         strr=sprintf(' sum(sum_square(GAMMAS{%d}+Xs{%d}-Y(%s,%s)))', j,j,inss,inss);
%   
%     else
%          [~,inss]=getTrippleInds(zeros(n), Cf(j,:) );
%         inss=mat2str(inss);
%       strr=sprintf('%s + sum(sum_square(GAMMAS{%d}+Xs{%d}-Y(%s,%s)))',strr, j,j,inss,inss);
%     end
%     
%     
% end
% cvx_begin sdp quiet
%     variables   Y(n,n)
%     
%     
% %     minimize( sum(sum_square(Fhat-A+A.')) + 0.5*beta*sum(sum_square(A-B)) +trace(Y'*(A-B)) + tau*(norm_nuc(A+A.') - trace(Al*(A+A.')*Bl')  )   )% norm(A-W,'fro')   0.5*beta*sum(sum_square(B-A+gamma))    + tau*(norm_nuc(Fhat) - trace(Al*(Fhat)*Bl')  ) 
%    % minimize( sum(sum_square(GAMMA1+X1-Y(1:9,1:9)))+sum(sum_square(GAMMA2+X2-Y(4:12,4:12)))+ sum(sum_square(GAMMA3+X3-Y(7:15,7:15)))+sum(sum_square(GAMMA4+X4-Y(10:18,10:18)))) %{+ 0.5*beta*sum(sum_square(B-A+gamma))%}  )%tau*(norm_nuc(Fhat) - trace(Al*(Fhat)*Bl') )  + tau*(norm_nuc(Fhat) - trace(Al*(Fhat)*Bl'))
%   minimize(eval(strr))
%    
%     subject to
%         Y-Y.'==0;
% %         lambda(1,n) == 1
%         for i = 3:3:n
%             Y(i-2:i,i-2:i) == 0;
%         end
%         
%     
% 
%         
%     
% cvx_end
% end


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
