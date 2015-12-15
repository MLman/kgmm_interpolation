function dS = dfds(gmms)
%DFDS returns the first derivative w.r.t S.
%
%
%   See Also: DFDPI, DFDM, L2DISTGMM, L2DISTGMMS

%   $ Hyunwoo J. Kim $  $ 2014/11/02 00:40:10 (CDT) $

    N = length(gmms);
    K = gmms{1}.NComponents;
    d = gmms{1}.NDimensions;
    dS = zeros(d, d, K, N-2);  % K Low vectors and N-2 free GMMs

    % PComponents is pi
    % Confusing indices. The first gmm and the last gmm are fixed.
    % dpi has two less vectors than #gmms.
    for i = 2:N-1
        for j =1:K
            for jj = 1:K
%                     if j == jj
%                         c = 2;
%                     else
%                         c = 4;
%                     end
                    dS(:,:,j,i-1) = dS(:,:,j,i-1) + 4*gmms{i}.PComponents(j)...
                    *gmms{i}.PComponents(jj)*dCdS(gmms, i, i, j, jj)...
                    -2*gmms{i}.PComponents(j)...
                    *gmms{i+1}.PComponents(jj)*dCdS(gmms, i, i+1, j, jj)...
                    -2*gmms{i}.PComponents(j)...
                    *gmms{i-1}.PComponents(jj)*dCdS(gmms, i, i-1, j, jj);
            end
        end
    end
end
    
function c = dCdS(gmms, i, ii, j, jj)
   % Equivalent and simpler implmenetation.
   % If i == ii &&  j == jj, c needs to be multiplied by 2.
   % However, we multiply 4 instead of 2 in dfds.
   m_i_j = gmms{i}.mu(j,:)';
   m_ii_jj = gmms{ii}.mu(jj,:)';
   S_i_j = gmms{i}.Sigma(:,:,j);
   S_ii_jj = gmms{ii}.Sigma(:,:,jj);
   c = -1/2*mvnpdf(m_i_j,m_ii_jj,S_i_j+S_ii_jj)*inv(S_i_j+S_ii_jj) ...
       +1/2*mvnpdf(m_i_j,m_ii_jj,S_i_j+S_ii_jj)*inv(S_i_j + S_ii_jj)*(m_i_j-m_ii_jj)*(m_i_j-m_ii_jj)'*inv(S_i_j + S_ii_jj)';
end