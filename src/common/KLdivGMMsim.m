function KLdiv = KLdivGMMsim(f, g, nsamples)
%KLDIVGMMSIM
%
%
%   See Also: XENTROPYGMMSIM, CROSSENTROPYGMMS

%   $ Hyunwoo J. Kim $  $ 2015/10/13 12:09:38 (CDT) $
    Xfromf = random(f, nsamples);
    Xfromg = random(g, nsamples);

    % PDF
    fxfromf = pdf(f,Xfromf);
    gxfromf = pdf(g,Xfromf);

    % Cross entropy is E_f[-log g(x)]
    % KL divergence 
    % 
    
    E_f_g = -mean(log(gxfromf));
    E_f_f = mean(log(fxfromf));

    KLdiv = E_f_g+E_f_f;
end