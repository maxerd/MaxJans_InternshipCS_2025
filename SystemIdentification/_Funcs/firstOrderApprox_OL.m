function varargout = firstOrderApprox_OL(data,opts)
% lw = 1.5;
% bodeRange = logspace(-6,-1,100);

tVec = 0:data.Ts:(length(data.y)*data.Ts-data.Ts);
approxSys = step_sysID(data.u',zeros(size(data.u))',data.y,tVec,360/data.Ts);

    if nargout>0
        varargout{1} = approxSys;
    else
        set(gcf,'position',[500 100 900 500])
        sgtitle(['Heater to thermal mass temperature, first order approximation, ',num2str(tVec(end)/3600),' hours'])
            subplot(211)
                simpleBodemag(approxSys,'Hz',opts.lw ,opts.bodeRange);hold on;grid minor
            subplot(212)
                simpleBodephase(approxSys,'Hz',opts.lw,'wrap' ,opts.bodeRange);hold on;grid minor
    end
end