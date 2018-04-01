function [isEqualStr, isEqual] = equalEps(a, b, epsilon)
isEqualStr = 'cannot be underfined';
if nargin == 2
    epsilon = 1e-7;
end
absErr = abs(a - b);
maxAB = max(abs(a), abs(b));
maxAB(maxAB == 0) = epsilon;
relErr = absErr ./ maxAB;
isEqualM = absErr.*relErr < epsilon;
isEqual = isempty(find(isEqualM == 0, 1));

if(isEqual)
    isEqualStr = 'Yes';
else
    isEqualStr = 'No';
end

end

