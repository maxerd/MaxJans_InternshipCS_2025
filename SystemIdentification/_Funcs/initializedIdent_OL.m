function varargout = initializedIdent_OL(data,initSys,opts)

tVec = 0:data.Ts:(length(data.y)*data.Ts-data.Ts);

optSS = ssestOptions('Focus','Prediction');

if strcmp(opts.focus,'sim')
    optSS = ssestOptions('Focus','Simulation');
elseif strcmp(opts.focus,'pred')
    optSS = ssestOptions('Focus','Prediction');
else
    warning('Focus option not known, using default (prediction focus)')
end

sys  = ssest(data,initSys,optSS);

    if nargout>0
        varargout{1} = sys;
    else
        set(gcf,'position',[500 100 900 500])
        sgtitle(['Heater to thermal mass temperature, first order approximation, ',num2str(tVec(end)/3600),' hours'])
            subplot(211)
                simpleBodemag(sys,'Hz',opts.lw ,opts.bodeRange);hold on;grid minor
            subplot(212)
                simpleBodephase(sys,'Hz',opts.lw,'wrap' ,opts.bodeRange);hold on;grid minor
    end

end