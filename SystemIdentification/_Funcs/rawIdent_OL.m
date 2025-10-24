function varargout = rawIdent_OL(data,opts)

tVec = 0:data.Ts:(length(data.y)*data.Ts-data.Ts);
warning('Focus option not possible, using default (prediction focus)')

sys  = ssest(data);

    if nargout>0
        varargout{1} = sys;
    else
        set(gcf,'position',[500 100 900 500])
        sgtitle(['Heater to thermal mass temperature, fixed order approximation, ',num2str(tVec(end)/3600),' hours'])
            subplot(211)
                simpleBodemag(sys,'Hz',opts.lw ,opts.bodeRange);hold on;grid minor
            subplot(212)
                simpleBodephase(sys,'Hz',opts.lw,'wrap' ,opts.bodeRange);hold on;grid minor
    end
end