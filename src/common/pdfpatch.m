function C = pdfpatch(gmm, xrange, yrange, N)
%PDFPATCH
%
%
%   See Also: PDF

%   $ Hyunwoo J. Kim $  $ 2015/03/22 23:34:52 (CDT) $

    minval = [xrange(1), yrange(1)];
    maxval = [xrange(2), yrange(2)];
    offset = (maxval-minval)*0.2;
    ix = linspace(minval(1)-offset(1), maxval(1)+offset(1), N);
    C = zeros(N,N);
    parfor i = 1:N
        for j =1:N
            iy = linspace(minval(2)-offset(2), maxval(2)+offset(2), N);
            C(j,i) = pdf(gmm,[ix(i) iy(j)]);
        end
    end
end