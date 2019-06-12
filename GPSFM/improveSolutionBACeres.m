function [ Ps2 ,Points3D , pointsInCameras] = improveSolutionBACeres( Ps,M)

W=M~=0;

pointsInCameras=cell(size(M,2),size(M,1)/2);
Points3D=zeros(4,size(M,2));

goodInds=ones(size(M,2),1);

for i=1:size(M,2)
    try
        dim = sum(  W(1:2:end,i));
        eqs = zeros(2*dim,4);
        Pscc = cell(1,dim);
        xsss = zeros(2,dim);
        counter=1;
        for j=1:size(M,1)/2
            %for each point go over its occurance
            
            if(W((j-1)*2+1,i))
                %points j
                pointsInCameras{i,j}=[M((j-1)*2+1:(j-1)*2+2,i);1];
                x1=M((j-1)*2+1:(j-1)*2+2,i);
                xsss(:,counter) = x1;
                Pscc{counter}=Ps{j};
                curEqs=[x1(1)*Ps{j}(3,:)-Ps{j}(1,:);x1(2)*Ps{j}(3,:)-Ps{j}(2,:)];
                eqs(2*(counter-1)+1:2*(counter-1)+2,:) = curEqs;
                counter=counter+1;
            end
        end
        [u,d,v]=svd(eqs);
        X=v(:,end);
        Points3D(:,i)=X(1:4)/X(4);
        
        Points3D(:,i)=X(1:4)/X(4);   
    catch
        goodInds(i)=0;
    end   
end
goodInds=logical(goodInds);
reffs=find(goodInds);
Points3D=Points3D(:,goodInds);


pointsInCamerasO=pointsInCameras;
pointsInCameras=cell(size(Points3D,2),size(M,1)/2);
for i=1:size(Points3D,2)
    for j=1:(size(M,1)/2)
        pointsInCameras{i,j}=pointsInCamerasO{reffs(i),j};
    end
end


num2d= sum(sum(~cellfun(@isempty,pointsInCameras),2));
[js,is]=find((~cellfun(@isempty,pointsInCameras))');
camPointmap=[js,is];
twoDPoints=zeros(num2d,2);

for i=1:num2d
    temp=pnorm(pointsInCameras{is(i),js(i)});
    twoDPoints(i,:)=temp(1:2)';
end

Psrep=zeros(length(Ps),12);
for i=1:length(Ps)
    Psrep(i,:)=[Ps{i}(:)]';
end


[Xsu,Psu]=BAceres(1,Points3D(1:3,:),camPointmap,Psrep',twoDPoints);
newXs=Points3D(1:3,:)'+Xsu';
newPs=Psrep+Psu';


Ps2=cell(size(newPs,1),1);
for i=1:size(newPs,1)
    Ps2{i}=reshape(newPs(i,:)',3,4);
end


Points3D = [newXs ones(size(newXs,1),1)]';



end

