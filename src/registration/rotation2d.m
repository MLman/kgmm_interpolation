function R = rotation2d(theta)
%ROTATION2D returns a rotation matrix in 2d. Theta is the angle to rotate
%in countclockwise.
%
%
%   See Also: 

%   $ Hyunwoo J. Kim $  $ 2015/04/01 12:27:48 (CDT) $

    R = [cos(theta), -sin(theta); 
        sin(theta),  cos(theta);];
end