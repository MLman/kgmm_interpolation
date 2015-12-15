function [A,b] = getaffinetransGauss(g1, data1, g2, data2, varargin)
%GETAFFINETRANSGAUSS finds affine transform between two Gaussian
%distribution based Cholesky decomposition and rotation in the isotrophic
%intermediate space.
%
%   Affine transform from g1 to g2.
%
%   getaffinetransGauss(g1, data1, g2, data2)
%   getaffinetransGauss(g1, data1, g2, data2, 'norotation')
%
%   G1, G2 are structures of Gaussian distributions.
%   G1.mu = a row vector in Rd
%   G1.Sigma = full matrix in R^dxd
%   data1.x = data points. Row vectors
%   data1.w = weight of data points.
%
%   See Also: GETAFFINETRANSGMM, AFFINETRANSGMM

%   $ Hyunwoo J. Kim $  $ 2015/04/01 21:10:23 (CDT) $

	L1 = chol(g1.Sigma,'lower');
    L2 = chol(g2.Sigma,'lower');
    
    x1c = data1.x-repmat(g1.mu, [size(data1.x,1),1]);
    x2c = data2.x-repmat(g2.mu, [size(data2.x,1),1]);

    invL1x1c = (x1c/L1');
    invL2x2c = (x2c/L2');
    
    if ~(nargin == 5 && strcmpi(varargin{1},'norotation'))
        [theta1, lx1c] = getthetafromxy(invL1x1c);
        [theta2, lx2c] = getthetafromxy(invL2x2c);
    end
    
    if nargin == 5 && strcmpi(varargin{1},'norotation')
        R = eye(2);
    elseif nargin == 4 || isempty(varargin{1}) ||strcmpi(varargin{1},'l2') %% Square of length
        % Default option.
        w1 = (lx1c.^2).*data1.w(:);
        w2 = (lx2c.^2).*data2.w(:);
    elseif nargin == 5 && strcmpi(varargin{1},'l1')
        w1 = lx1c.*data1.w(:);
        w2 = lx2c.*data2.w(:);
    elseif nargin == 5 && strcmpi(varargin{1},'theta')
        w1 = data1.w(:);
        w2 = data2.w(:);
    elseif nargin == 5 && strcmpi(varargin{1},'principalD') %% Principal direction
        PD1 = getPrincipalDirection(data1.x);
        PD2 = getPrincipalDirection(data2.x);
        invL1PD1 = L1\PD1;
        invL2PD2 = L2\PD2;
        [thetaPD1inI, ~] = getthetafromxy(invL1PD1');
        [thetaPD2inI, ~] = getthetafromxy(invL2PD2');
        R = rotation2d(thetaPD2inI-thetaPD1inI);
    end
    
    if nargin <=4 || ~( strcmpi(varargin{1},'norotation') || strcmpi(varargin{1},'principalD'))
        mtheta1 = sum(theta1.*w1)/sum(w1);
        mtheta2 = sum(theta2.*w2)/sum(w2);
        R = rotation2d(mtheta2-mtheta1);
    end
    A = L2*R/L1; % Double check
    fprintf('getaffinetransGauss Rotation needs to be check.')
    b = -A*g1.mu'+g2.mu';

end