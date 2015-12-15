function dpi = dfdpi(gmms)
%DFDPI returns the first derivative w.r.t pi (weights).
%
%
%   See Also: f

%   $ Hyunwoo J. Kim $  $ 2014/11/01 23:30:25 (CDT) $

    N = length(gmms);
    K = gmms{1}.NComponents;

    dpi = zeros(N-2,K);  % Low vectors
    mydpi = zeros(N-2,K);


%     my11dpi = zeros(N-2,K);
    my111dpi = zeros(N-2,K);
%     my22dpi = zeros(N-2,K);
    my222dpi = zeros(N-2,K);
%     my33dpi = zeros(N-2,K);
    my333dpi = zeros(N-2,K);

    % PComponents is pi
    % Confusing indices. The first gmm and the last gmm are fixed.
    % dpi has two less vectors than #gmms.
    for i = 2:N-1

        for j =1:K
            for jj = 1:K
                dpi(i-1,j) = dpi(i-1,j) + 4*gmms{i}.PComponents(jj)*dCdpi(gmms, i, i, j, jj)...   
                    - 2*gmms{i+1}.PComponents(jj)*dCdpi(gmms, i, i+1, j, jj)... 
                    - 2*gmms{i-1}.PComponents(jj)*dCdpi(gmms, i, i-1, j, jj);   
            end
        end

      
%         for j = 1:K
%           my11dpi(i-1,j) = sum(4*gmms{i}.PComponents'.*(mvnpdf(gmms{i}.mu(j,:),gmms{i}.mu,bsxfun(@plus,gmms{i}.Sigma,gmms{i}.Sigma(:,:,j)))));
%         end  
        
        D = size(gmms{i}.mu,2);
        
        tmpA = num2cell(gmms{i}.mu,[D,K])';
        tmpB = mat2cell(repmat(gmms{i}.mu,[K,1]),K*ones(1,K))';
        tmpC = mat2cell(repmat(gmms{i}.Sigma,[K,1,1])+repmat(reshape(gmms{i}.Sigma,[D,D*K])',[1,1,K]),D*ones(1,K))';
        my111dpi(i-1,:) = sum(bsxfun(@times,cell2mat(cellfun(@mvnpdf,tmpA,tmpB,tmpC,'UniformOutput',false)),4*gmms{i}.PComponents'),1);
        
%         for j = 1:K
%           my22dpi(i-1,j) = sum(-2*gmms{i+1}.PComponents'.*(mvnpdf(gmms{i}.mu(j,:),gmms{i+1}.mu,bsxfun(@plus,gmms{i+1}.Sigma,gmms{i}.Sigma(:,:,j)))));
%         end  
        
        tmpB = mat2cell(repmat(gmms{i+1}.mu,[K,1]),K*ones(1,K))';
        tmpC = mat2cell(repmat(gmms{i+1}.Sigma,[K,1,1])+repmat(reshape(gmms{i}.Sigma,[D,D*K])',[1,1,K]),D*ones(1,K))';
        my222dpi(i-1,:) = sum(bsxfun(@times,cell2mat(cellfun(@mvnpdf,tmpA,tmpB,tmpC,'UniformOutput',false)),-2*gmms{i+1}.PComponents'),1);
        
%         for j = 1:K
%           my33dpi(i-1,j) = sum(-2*gmms{i-1}.PComponents'.*(mvnpdf(gmms{i}.mu(j,:),gmms{i-1}.mu,bsxfun(@plus,gmms{i-1}.Sigma,gmms{i}.Sigma(:,:,j)))));
%         end
        
        tmpB = mat2cell(repmat(gmms{i-1}.mu,[K,1]),K*ones(1,K))';
        tmpC = mat2cell(repmat(gmms{i-1}.Sigma,[K,1,1])+repmat(reshape(gmms{i}.Sigma,[D,D*K])',[1,1,K]),D*ones(1,K))';
        my333dpi(i-1,:) = sum(bsxfun(@times,cell2mat(cellfun(@mvnpdf,tmpA,tmpB,tmpC,'UniformOutput',false)),-2*gmms{i-1}.PComponents'),1);
        
%         mydpi(i-1,:) = my11dpi(i-1,:)+my22dpi(i-1,:)+my33dpi(i-1,:);
        mydpi(i-1,:) = my111dpi(i-1,:)+my222dpi(i-1,:)+my333dpi(i-1,:);
        
    end 
    assert(all(all(sum((mydpi(:)-dpi(:)).^2)<1e-10)))
end

function c = dCdpi(gmms, i, ii, j, jj)
   m_i_j = gmms{i}.mu(j,:);
   m_ii_jj = gmms{ii}.mu(jj,:);
   S_i_j = gmms{i}.Sigma(:,:,j);
   S_ii_jj = gmms{ii}.Sigma(:,:,jj);
   c = mvnpdf(m_i_j,m_ii_jj,S_i_j+S_ii_jj);
end