function [isEqualStr, isEqual] = equalEps(a, b, epsilon)
isEqualStr = 'cannot be underfined';
if nargin == 2
    epsilon = 1e-7;
end
absErr = abs(a - b);
relErr = absErr./max(abs(a), abs(b));
isEqualM = absErr.*relErr < epsilon;
isEqual = isempty(find(isEqualM == 0, 1));
if(isEqual)
    isEqualStr = 'Yes';
else
    isEqualStr = 'No';
end

end

