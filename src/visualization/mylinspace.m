function [ix, iy] = mylinspace(range, N)
%MYLINSPACE returns linspace(range(i,1), range(i,2), N).
%
%
%   See Also: LINSPACE, DETERMINE_BORDER

%   $ Hyunwoo J. Kim $  $ 2015/04/01 16:06:33 (CDT) $

    ix = linspace(range(1,1), range(1, 2), N);
    iy = linspace(range(2,1),  range(2, 2), N);
   
    
end