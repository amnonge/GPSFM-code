function [ FN ] = normalizeF( FN,normalization )

for i = 1:3:size(FN,1)-2
    for j = 1:3:size(FN,2)-2
        FN(i:i+2,j:j+2) = normalization'*FN(i:i+2,j:j+2)*normalization;
    end
end

end

