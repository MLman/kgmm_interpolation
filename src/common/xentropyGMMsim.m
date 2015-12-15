function E_f_g = xentropyGMMsim(f, g, nsamples)
%XENTROPYGMMSIM computes cross entropy H(f,g) = - \sum_x f(x) log g(x)
%
%
%   See Also: KLDIVGMMSIM, CROSSENTROPYGMMS

%   $ Hyunwoo J. Kim $  $ 2015/10/13 11:57:00 (CDT) $   

    Xfromf = random(f, nsamples);
    gxfromf = pdf(g,Xfromf);
    E_f_g = -mean(log(gxfromf));
end