function R = getRotationfromAffine(A)
%GETROTATIONFROMAFFINE returns rotation matrices from affine transform
%matrix A.
%
%
%   See Also:

%   $ Hyunwoo J. Kim $  $ 2015/04/13 14:37:29 (CDT) $
    R = sqrtm(A*A')\A;
end