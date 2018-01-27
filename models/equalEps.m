function [isEqualStr, isEqual] = equalEps(a, b, epsilon)
if nargin == 2
    epsilon = eps('single');
end
isEqual = sum(sum(abs((a - b) ./ max(a, b)))) < epsilon;
if(isEqual)
    isEqualStr = 'Yes';
else
    isEqualStr = 'No';
end
end

