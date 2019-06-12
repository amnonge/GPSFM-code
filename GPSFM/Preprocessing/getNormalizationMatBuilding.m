function [normMatA,normMatB] = getNormalizationMatBuilding(xa,xb)

Points_mean = mean(xa,2);
c_u = Points_mean(1);
c_v = Points_mean(2);
translated_points = xa - repmat(Points_mean,1,size(xa, 2));
temp_mean_from_origin = (mean(sum((translated_points.^2),1)))/1;
%compute scale
scale = sqrt(1/(temp_mean_from_origin));
%compute Transformation matrix as described in the project page
normMatA = [scale 0 0 ;0 scale 0; 0 0 1] * [1 0 -c_u; 0 1 -c_v; 0 0 1];

normMatA = normMatA/normMatA(3,3);

Points_mean = mean(xb,2);
c_u = Points_mean(1);
c_v = Points_mean(2);
translated_points = xb - repmat(Points_mean,1,size(xb, 2));
temp_mean_from_origin = (mean(sum((translated_points.^2),1)))/1;
%compute scale
scale = sqrt(1/(temp_mean_from_origin));
%compute Transformation matrix as described in the project page
normMatB = [scale 0 0 ;0 scale 0; 0 0 1] * [1 0 -c_u; 0 1 -c_v; 0 0 1];

normMatB = normMatB/normMatB(3,3);


end