% function varargout = simpleBodemag_MC(sys,freqUnit,bodeRange,varargin)
function varargout = simpleBodemag_MC(sys,freqUnit,varargin)

n_MC = size(sys,2);

rangeGiven = false;

if nargin>2
    for i=1:(nargin-2)
        if strcmp(class(varargin{i}),'double')
            bodeRange = varargin{i};
            rangeGiven = true;
        end
        if strcmp(class(varargin{i}),'char')
            lineDat = varargin{i};
            if lineDat(1)=='b'
                lineColor = 'b';
                boundColor = [0.8 0.8 1];
            elseif lineDat(1)=='r'
                lineColor = 'r';
                boundColor = [1 0.8 0.8];
            elseif lineDat(1)=='g'
                lineColor = 'g';
                boundColor = [0.87 1 0.87];
            elseif lineDat(1)=='k'
                lineColor = 'k';
                boundColor = [0.8 0.8 0.8];
            elseif lineDat(1)=='m'
                lineColor = 'm';
                boundColor = [1 0.85 1];
            elseif lineDat(1)=='c'
                lineColor = 'c';
                boundColor = [0.85 1 0.95];
            elseif lineDat(1)=='y'
                lineColor = 'y';
                boundColor = [0.97 1 0.86];
            end
        end
    end
else
    lineColor = 'b';
    boundColor = [0.8 0.8 1];
end

% avgMag = zeros(size(bodeRange));
% difMag = zeros(size(bodeRange));

for i=1:n_MC
    if rangeGiven
        [mag,~,wout] = bode(sys{i},bodeRange*2*pi);
    else
        [mag,~,wout] = bode(sys{i});
    end
    resp(:,i) = squeeze(mag);
end
avgMag = mean(resp,2);
difMag = max(abs(resp-repmat(avgMag,1,n_MC))')';
lowerBound = avgMag-difMag;
lowerBound = ((lowerBound>0)).*lowerBound+1e-12;
upperBound = avgMag+difMag;


if strcmp(freqUnit,'Hz')
    wout = wout/2/pi;
end
    fill([(wout);flip((wout))], [db(lowerBound);flip(db(upperBound))], boundColor, 'EdgeColor', 'none');hold on
    plot((wout),db(avgMag),lineColor,LineWidth=1.5);hold on
    % plot(wout,db(lowerBound))
    % plot(wout,abs(avgMag-difMag))
    grid minor
    ylim([db(min(avgMag)) db(max(upperBound))])
    set(gca,'XScale','log')


end