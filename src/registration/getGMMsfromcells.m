function gmms = getGMMsfromcells(gmmcells, coordinate)
%GETGMMSFROMCELLS returns a subset of GMMS from GMMCELLS AT COORDINATES.
%
%   COORDINATE is a set of row vectors.   
%
%   See Also: EAP_AFFINE_DEMO

%   $ Hyunwoo J. Kim $  $ 2015/04/13 17:26:52 (CDT) $
    ngmms = size(coordinate,1);
    dim  = size(coordinate,2);
    gmms = cell(1,size(coordinate,1));
    if dim == 2
        for i = 1:ngmms
            gmms{i} = gmmcells{coordinate(i,1),coordinate(i,2)};
        end
    else
        for i = 1:ngmms
            gmms{i} = gmmcells{coordinate(i,1),coordinate(i,2),coordinate(i,3)};
        end
    end

end