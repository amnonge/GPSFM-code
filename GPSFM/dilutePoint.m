function [newM,keep_pts,rm_pts] = dilutePoint(M)

if size(M,2)> 20000
    param = 4;
else
    param = 3;
end
visible = M ~= 0;
visible = visible(1:2:end,:) & visible(2:2:end,:);

rm_pts = sum(visible, 1) < param;
keep_pts = ~rm_pts;
newM = M(:,keep_pts);

end