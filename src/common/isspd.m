function T = isspd(mx,varargin)
%ISSPD checks if mx (2d or 3d) is a set of spd matrices.
%   
%   MX is dxd or dxdxn matrix. 
%   T is 1 if spd which means eigen values are greater than tolerance c.
%   Otherwise T is 0. For dxdxn matrix, T is a vector.
%   
%   T = isspd(mx) returns 1 or n boolean values with default tolerance 0.
%   T = isspd(mx, 1e-10) returns boolean values with tolerance c=1e-10.
%
%   See Also: MYISREAL

%   $ Hyunwoo J. Kim $  $ 2015/01/28 20:16:18 (CST) $

% Check matrices are symmetric positive definite.

    if nargin  == 2
        c = varargin{1};
    else
        c = 0;
    end
   
    T = zeros(size(mx,3),1);
    for i=1:size(mx,3)
        if myisreal(mx(:,:,i))
            T(i) = sum(eig(mx(:,:,i)) <= 0 +c) ==0;
        else
            T(i) = false;
        end
    end
    T = logical(T);
end