function myobj = obj2structGMM(obj)
%OBJ2STRUCTUREGMM converts GMM object to my structure.
%
%   GMM objects don't allow to modify field values. 
%
%   See Also: STRUCT2OBJGMM

%   $ Hyunwoo J. Kim $  $ 2014/10/26 17:51:50 (CDT) $

    if isstruct(obj) % To prevent redundant type conversion.
        myobj = obj;
        return 
    end

    ver = version;
    if ~strcmp(ver(end-5:end-1),'2014b')
        ofields = fieldnames(obj);
        myobj = [];

        for i = 1:length(ofields)
            eval(sprintf('myobj.%s = obj.%s;',ofields{i}, ofields{i}));
        end

    else  % After 2014b, the field names have changed.
        myobj.NDimensions = obj.NumVariables;
        myobj.DistName = obj.DistributionName;
        myobj.NComponents = obj.NumComponents;
        myobj.PComponents = obj.ComponentProportion;
        myobj.SharedCov = obj.SharedCovariance;
        myobj.Iters = obj.NumIterations;
        myobj.RegV = obj.RegularizationValue;
        myobj.NlogL = obj.NegativeLogLikelihood;
        myobj.CovType = obj.CovarianceType;
        myobj.mu = obj.mu;
        myobj.Sigma = obj.Sigma;
        myobj.AIC = obj.AIC;
        myobj.BIC = obj.BIC;
        myobj.Converged = obj.Converged;
    end
end


 