function preProcess(dataset)
load(dataset);
if strcmp('Pumpkin 2',dataset)
    dataset
    indices=[1:184 186:196];
    indicesx=indices*2-1;
    indicesy=indices*2;
    Mc=zeros(length(indices)*2,size(M,2));
    Mc(1:2:end,:)=M(indicesx,:);
    Mc(2:2:end,:)=M(indicesy,:);
    M=Mc;
    
    M=M(:,(sum(abs(M)>10^-5,1)>=4));
    FNc=zeros(3*length(indices),3*length(indices));
    Fc=zeros(3*length(indices),3*length(indices));
    for i=1:length(indices)-1
        for j=i+1:length(indices)
            FNc(3*i-2:3*i,3*j-2:3*j)=FN(3*indices(i)-2:3*indices(i),3*indices(j)-2:3*indices(j));
            FNc(3*j-2:3*j,3*i-2:3*i)=FN(3*indices(j)-2:3*indices(j),3*indices(i)-2:3*indices(i));
            Fc(3*i-2:3*i,3*j-2:3*j)=F(3*indices(i)-2:3*indices(i),3*indices(j)-2:3*indices(j));
            Fc(3*j-2:3*j,3*i-2:3*i)=F(3*indices(j)-2:3*indices(j),3*indices(i)-2:3*indices(i));
        end
    end
    FN=FNc;
    F=Fc;
    pointMatchesGroundc=cell(length(indices),length(indices),2);
    pointMatchesInliersc=cell(length(indices),length(indices),2);
    
    for i=1:length(indices)-1
        for j=i+1:length(indices)
            pointMatchesGroundc{i,j,1}=pointMatchesGround{indices(i),indices(j),1};
            pointMatchesGroundc{i,j,2}=pointMatchesGround{indices(i),indices(j),2};
            
            pointMatchesInliersc{i,j,1}=pointMatchesInliers{indices(i),indices(j),1};
            pointMatchesInliersc{i,j,2}=pointMatchesInliers{indices(i),indices(j),2};
        end
    end
    pointMatchesGround=pointMatchesGroundc;
    pointMatchesInliers=pointMatchesInliersc;
end

[ FN ] = normalizeForbineusNorm( FN );
FN(isnan(FN)) = 0;

pointMatchesInliers = converPointInlier(pointMatchesInliers);
M=M(:,(sum(abs(M)>10^-5,1)>=4));
save([dataset ' new'],'pointMatchesInliers','FN','M')
end
