clear
addpath(genpath('.'));
dataset = 'Example';
load([dataset '.mat'])


width = 1296;
hight = 1936;
image_size = repmat([width,hight]',1,length(P));
centers =repmat([width/2,hight/2]',1,length(P));
points = u_uncalib.points;
index = u_uncalib.index;
measurements = zeros(2*length(points),u_uncalib.pointnr);
for i = 1:length(points)
    for j = 1:size(points{i},2)
        measurements(2*i-1:2*i,index{i}(j)) = points{i}(1:2,j);
    end
end

measurements(isnan(measurements))=0;




camsNum=size(measurements,1)/2;

FN=zeros(camsNum*3,camsNum*3);
pointMatchesInliers=cell(3,3,2);
pointMatchesGround=cell(3,3,2);
M=measurements;
for i=1:size(measurements,1)/2-1
    for j=i+1:size(measurements,1)/2
        ttemp=[i j]*2;
        ttemp1=[i j]*2-1;
        ttt=[[i j] [i j]];
        ttt(1:2:end)=ttemp1;
        ttt(2:2:end)=ttemp;
        indsTriplets= find(sum(M(ttt,:)~=0)>3);
             if length(indsTriplets)<8
                 continue;
             end
         xa=M(ttt(1:2),indsTriplets);
         xa=[xa ; ones(1,size(xa,2))];
         xb=M(ttt(3:4),indsTriplets);
         xb=[xb ; ones(1,size(xb,2))];
        %  
         [normMatA,normMatB] = getNormalizationMatBuilding(xa,xb);
         xa = normMatA * xa;
         xb = normMatB * xb;

            [FijN,inliersIndex] = estimateFundamentalMatrix(xa(1:2,:)',...
            xb(1:2,:)','Method','MSAC',...
            'NumTrials',2000,'InlierPercentage',50);

    if size(xa(:,inliersIndex),2)> 50
      FijN = normMatB' * FijN * normMatA;
      FijN=FijN'/FijN(3,3);
      xa = inv(normMatA) * xa;
      xb = inv(normMatB) * xb;
    else
        FijN = 0;
    end
        FN(3*i-2:3*i,3*j-2:3*j)=FijN;
        FN(3*j-2:3*j,3*i-2:3*i)=FijN';

        FijG = FijN;
        pointMatchesInliers{i,j,1}=xa;
        pointMatchesInliers{i,j,2}=xb;
        xij_a = xa;
        xij_b = xb;
        pointMatchesGround{i,j,1}=xij_a;
        pointMatchesGround{i,j,2}=xij_b;

    end
end

if strcmp('Pumpkin 2',dataset)

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
save([dataset ' pro'],'pointMatchesInliers','FN','M');

