function fs = myplotgmm1D(gmm, xinterval, varargin)
%MYPLOTGMM1D plot 1D gaussian mixture models using built-in gmdistribution
%object.
%
%   See Also: GMDISTRIBUTION, MYPLOTGMM2D

%   $ Hyunwoo J. Kim $  $ 2014/11/02 20:56:08 (CST) $

    

    if length(varargin) >= 1
        mycolor = varargin{1};
    else
        mycolor = rand(1,3);
    end
    if length(varargin) >= 2
        delta = varargin{2};
    else
        delta = 0.1;
    end

    x = (xinterval(1):delta:xinterval(2))';
    
    
%    figure
    hold on;
    
    for i=1:gmm.NComponents
        
        gm = gmdistribution(gmm.mu(i,:),gmm.Sigma(:,:,i), 1);
        fs = gm.pdf(x)*gmm.PComponents(i);
        if i==1
            h2=plot(x,fs,'-.', 'color', [0.5 0.5 0.5]);
            h3=plot(gmm.mu(i,:),0,'or');
        else
            h2=plot(x,fs,'-.', 'color', [0.5 0.5 0.5]);
            h3=plot(gmm.mu(i,:),0,'or');
        end
    end
    
    gm = gmdistribution(gmm.mu,gmm.Sigma, gmm.PComponents);
    fs = gm.pdf(x);
    h1 = plot(x,fs, '--', 'Color',mycolor);
    legend([h1,h2,h3], 'GMM', 'Components', 'Mus');
    hold off;
end