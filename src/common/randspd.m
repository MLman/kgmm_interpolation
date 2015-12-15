function P = randspd(n, varargin)
%RANDSPD generates a random spd matrix.
%
%
%   See Also: WISHRND

%   $ Hyunwoo J. Kim $  $ 2015/09/15 13:08:44 (CDT) $

    while true
        if nargin == 2
            c = varargin{1};
        else
            c = 3;
        end
        P = c*(rand(n)-0.5);
        P = P*P';

        if isspd(P,1e-1)
            break
        end
    end
end

