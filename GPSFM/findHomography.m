function [ Hf ,Psa] = findHomography( Ps,Pso )
Ps_s1=Ps{1};
Ps_s2=Ps{2};
Ps_d1=Pso{1};
Ps_d2=Pso{2};


csnosym=zeros(24,18);
csnosym(1:3,1)=-Ps_s1(:,1);
csnosym(13:15,1)=-Ps_s2(:,1);


csnosym(4:6,2)=-Ps_s1(:,1);
csnosym(16:18,2)=-Ps_s2(:,1);

csnosym(7:9,3)=-Ps_s1(:,1);
csnosym(19:21,3)=-Ps_s2(:,1);

csnosym(10:12,4)=-Ps_s1(:,1);
csnosym(22:24,4)=-Ps_s2(:,1);




csnosym(1:3,5)=-Ps_s1(:,2);
csnosym(13:15,5)=-Ps_s2(:,2);

csnosym(4:6,6)=-Ps_s1(:,2);
csnosym(16:18,6)=-Ps_s2(:,2);

csnosym(7:9,7)=-Ps_s1(:,2);
csnosym(19:21,7)=-Ps_s2(:,2);

csnosym(10:12,8)=-Ps_s1(:,2);
csnosym(22:24,8)=-Ps_s2(:,2);


csnosym(1:3,9)=-Ps_s1(:,3);
csnosym(13:15,9)=-Ps_s2(:,3);

csnosym(4:6,10)=-Ps_s1(:,3);
csnosym(16:18,10)=-Ps_s2(:,3);

csnosym(7:9,11)=-Ps_s1(:,3);
csnosym(19:21,11)=-Ps_s2(:,3);

csnosym(10:12,12)=-Ps_s1(:,3);
csnosym(22:24,12)=-Ps_s2(:,3);


csnosym(1:3,13)=-Ps_s1(:,4);
csnosym(13:15,13)=-Ps_s2(:,4);

csnosym(4:6,14)=-Ps_s1(:,4);
csnosym(16:18,14)=-Ps_s2(:,4);

csnosym(7:9,15)=-Ps_s1(:,4);
csnosym(19:21,15)=-Ps_s2(:,4);

csnosym(10:12,16)=-Ps_s1(:,4);
csnosym(22:24,16)=-Ps_s2(:,4);

csnosym(1:12,17)=Ps_d1(:);
csnosym(13:24,18)=Ps_d2(:);

cs=csnosym;

[u,d,v]=svd(cs);
size(v);
Hf=[v(1,end) v(2,end) v(3,end) v(4,end);v(5,end) v(6,end) v(7,end) v(8,end);
    v(9,end) v(10,end) v(11,end) v(12,end);v(13,end) v(14,end) v(15,end) v(16,end)];
Psa=cell(length(Ps),1);
for i=1:length(Ps)
    Psa{i}=Ps{i}*Hf;
    Pso{i}*v(16+i,end);
end



end

