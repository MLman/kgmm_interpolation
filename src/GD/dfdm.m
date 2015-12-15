function dm = dfdm(gmms)
%DFDM returns the first derivative w.r.t m.
%
%
%   See Also: DFDPI, L2DISTGMM, L2DISTGMMS

%   $ Hyunwoo J. Kim $  $ 2014/10/30 02:12:53 (CDT) $

    N = length(gmms);
    K = gmms{1}.NComponents;
    d = gmms{1}.NDimensions;
    dm = zeros(K, d, N-2);  % K Low vectors and N-2 free GMMs

    % PComponents is pi
    % Confusing indices. The first gmm and the last gmm are fixed.
    % dpi has two less vectors than #gmms.
    for i = 2:N-1
        for j =1:K
            for jj = 1:K
                    dm(j,:,i-1) = dm(j,:,i-1) + 4*gmms{i}.PComponents(j)...
                    *gmms{i}.PComponents(jj)*dCdm(gmms, i, i, j, jj)...
                    -2*gmms{i}.PComponents(j)...
                    *gmms{i+1}.PComponents(jj)*dCdm(gmms, i, i+1, j, jj)...
                    -2*gmms{i}.PComponents(j)...
                    *gmms{i-1}.PComponents(jj)*dCdm(gmms, i, i-1, j, jj);
            end
        end
    end
end
    
function c = dCdm(gmms, i, ii, j, jj)
   % A little speed up?
   if i ==ii && j == jj 
       c = 0;
       return 
   end
   m_i_j = gmms{i}.mu(j,:)';
   m_ii_jj = gmms{ii}.mu(jj,:)';
   S_i_j = gmms{i}.Sigma(:,:,j);
   S_ii_jj = gmms{ii}.Sigma(:,:,jj);
   c = -(mvnpdf(m_i_j,m_ii_jj,S_i_j+S_ii_jj)*inv(S_i_j+S_ii_jj)*(m_i_j - m_ii_jj))';
end