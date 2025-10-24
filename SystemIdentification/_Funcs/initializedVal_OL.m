function varargout = initializedVal_OL(sys,valData)

valSim = compareRMS(valData,'Initalized Model',sys,inf);
valPrd = compareRMS(valData,'Initalized Model',sys,1);

    if nargout>0
        varargout{1}.sim = valSim;
        varargout{1}.prd = valPrd;
    else
        figure;
        compareRMS(valData,{'Initalized Model'},sys,inf);
        figure;
        compareRMS(valData,{'Initalized Model'},sys,1);
        figure;
        resid(valData,sys);grid minor
    end

end