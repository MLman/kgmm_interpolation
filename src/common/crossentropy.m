function H = crossentropy(m_p, S_p, m_q, S_q)
%CROSSENTROPY calculates between two gaussians.
%
%   m_p, m_q are means. Column vectors 
%   S_p, S_q are covariance matrices. Full matrix.
%   H(p,q) = E_p[-log q]
%   H(p,q) = 1/2*(k log 2pi + log |Sigma_q| +tr(Simga_q^{-1}, Sigma_p)
%              +(m_p - m_q)'S_q(m_p - m_q)
%
%
%   See Also: EM_GMM_CLOSEST_FULL

%   $ Hyunwoo J. Kim $  $ 2015/03/27 00:42:56 (CDT) $

    k = length(m_p); % Dimension of m_p

%    H = 1/2*( k*log(2*pi) + log(det(S_q)) + trace(inv(S_q)*S_p) ...
%        + (m_p-m_q)'*inv(S_q)*(m_p-m_q));

    [V,D] = eig(S_q);
    invS_q = V*diag(1./diag(D))*V';
    H = 1/2*( k*log(2*pi) + sum(log(diag(D))) + trace(invS_q*S_p) ...
        + (m_p-m_q)'*invS_q*(m_p-m_q));
end