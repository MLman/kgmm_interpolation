function mystructcells = obj2structGMMs(objcells)
%OBJ2STRUCTUREGMM converts GMM object to my structure.
%
%   GMM objects don't allow to modify field values. 
%
%   See Also: STRUCT2OBJGMM

%   $ Hyunwoo J. Kim $  $ 2014/10/26 17:51:50 (CDT) $
    mystructcells = cell(size(objcells));
    
    for i = 1:numel(objcells)
        if isempty(objcells{i})
            continue
        end
        mystructcells{i} = obj2structGMM(objcells{i});
    end
end


 