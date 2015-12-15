function gmmobj = struct2objGMM(mygmm)
%STRUCT2OBJGMM converts my gmm structure to gmm object compatible 
%   with Matlab built-in code.
%
%
%   See Also: OBJ2STRUCTUREGMM

%   $ Hyunwoo J. Kim $  $ 2015/02/24 17:07:09 (CST) $
    gmmobj = gmdistribution(mygmm.mu, mygmm.Sigma, mygmm.PComponents);
end
