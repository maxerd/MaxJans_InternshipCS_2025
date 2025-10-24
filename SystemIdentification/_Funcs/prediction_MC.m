function [error,varargout] = prediction_MC(validData_id,validData_tVec,sys,varargin)

n_MC = size(sys,2);

if nargin>3
    for i=1:(nargin-3)
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

for i=1:n_MC
    pred = compare(validData_id,sys{i},1);
    pred_val(:,i) = pred.y;
end

avgPred = mean(pred_val,2);
difPred = max(abs(pred_val-repmat(avgPred,1,n_MC))')';

if ~nargout
    fill([(validData_tVec) flip((validData_tVec))], [(avgPred-difPred);flip((avgPred+difPred))], boundColor, 'EdgeColor', 'none');hold on
    plot((validData_tVec),(avgPred),lineColor,LineWidth=1.5);hold on
    plot((validData_tVec),(validData_id.y),'Color',[0 1 0],LineWidth=1.5);hold on
    grid minor
end

error.avg = rmse(avgPred,validData_id.y);
error.all = rmse(pred_val,validData_id.y);

if nargout>1
    varargout{1}.avg = avgPred;
    varargout{1}.all = pred_val;
end



