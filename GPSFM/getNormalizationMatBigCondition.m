function normMat = getNormalizationMatBigCondition(stdParam,M)
newpoints = cell(size(M,1)/2,1);
for i = 1:(size(M,1)/2)
    
    relcur=M(i*2-1:i*2,:);
    relcur=relcur(:,sum(abs(relcur))>0);
    newpoints{i}=relcur;
end
normMatTotal = zeros(size(M,1)/2);
ttts=zeros(length(newpoints),1);
for i = 1:length(newpoints)
    allPoints = newpoints{i};
    Points_mean = mean(allPoints,2);
    translated_points = allPoints - repmat(Points_mean,1,size(allPoints, 2)); 
    tt=mean(translated_points.^2,2);
    ttts(i) =(min(tt(1)/tt(2),tt(2)/tt(1)));
end
if mean(ttts)<0.4
    anisotropic=true;
else
    anisotropic=false;
end
for i = 1:length(newpoints)
    allPoints = newpoints{i};
    Points_mean = mean(allPoints,2);
    c_u = Points_mean(1);
    c_v = Points_mean(2);
    translated_points = allPoints - repmat(Points_mean,1,size(allPoints, 2));
    if anisotropic==true
        [trans, ~] = normtrans(allPoints, false);
    else
        temp_mean_from_origin = (mean(sum((translated_points.^2),1)))/stdParam;
        scale = sqrt(1./(temp_mean_from_origin));
        trans=[scale 0 0 ;0 scale 0; 0 0 1] * [1 0 -c_u; 0 1 -c_v; 0 0 1];
    end

    normMat=inv(trans);
    close all;
    normMat = normMat/normMat(3,3);
   
    normMatTotal(3*i-2:3*i,3*i-2:3*i) = normMat;
end
normMat = normMatTotal;
end