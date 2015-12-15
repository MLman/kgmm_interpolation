function mxstacknew = diag3D(mxstack)
%DIAG3D is a 3D version diag function.
%
%   See Also:

%   $ Hyunwoo J. Kim $  $ 2015/03/25 00:34:54 (CDT) $
    dims = size(mxstack);
    assert(length(dims) == 3);
    if dims(1) == 1 && dims(1) < dims(2)
        mxstacknew = zeros(dims(2),dims(2),dims(3));
        for i =1:size(mxstack,3)
            mxstacknew(:,:,i) = diag(mxstack(1,:,i));
        end
    elseif dims(1) > 1 && dims(1) == dims(2)
        mxstacknew = zeros(1,dims(2),dims(3));
        for i =1:size(mxstack,3)
            mxstacknew(1,:,i) = diag(mxstack(:,:,i));
        end
    else
        mxstacknew = mxstack;
    end
end
