function newMat = converPointInlier(pointMatchesInliers)
newMat = zeros(size(pointMatchesInliers));
for i = 1:size(pointMatchesInliers,1)
    for j = 1:size(pointMatchesInliers,2)
        newMat(i,j)= size(pointMatchesInliers{i,j,1},2);
    end
end
end