function [answerStr, isPass] = equalEpsN(varargin)
for i=1:nargin-1
    [answerStr, isPass] = equalEps(varargin{i}, varargin{i+1});
    if ~isPass
        break;
    end
end
end