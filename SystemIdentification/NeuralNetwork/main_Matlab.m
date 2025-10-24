
%% Set Python environment
clc
% pyenv('Version','C:\python27.exe');
pyenv('Version',"C:\Users\maxja\AppData\Local\Programs\Python\Python310\python.exe");


%% --- Train ---
clc
rmse = py.main_NN_V3.train_model('file.csv', 'lstm_model_internship.keras', int32(10), int32(10));
disp(rmse);

%% --- Predict ---
clear y_val x_val X_test

% file = dataImportPython("C:\Users\maxja\Documents\(4)School\Master\Q9_Internship\MaxJans_InternshipCS_2025\SystemIdentification\NeuralNetwork\file.csv", [2, Inf]);
dataDir_OL = 'C:\Users\maxja\Documents\(4)School\Master\Q9_Internship\matlabFiles\measurements\processedData\p__CB_none__ST_ms__SM_24w__SA_12w__DT_250905__MD_9h__WT_no__DS_100.mat'; % Open loop data used for identification

dat = load(dataDir_OL);

udat = dat.Watt;
ydat = dat.tempTM;

% y_val = ydat(length(udat)*0.15:end);
% x_val = udat(length(udat)*0.15:end);
% y_val = ydat(1:length(udat)*0.3);
% x_val = udat(1:length(udat)*0.3);
y_val = ydat;
x_val = udat;

seq_length = 10;
for i=1:(length(x_val)-seq_length-1)
    X_test(i,:,1) = x_val(i+1:i+seq_length+1);
    X_test(i,:,2) = y_val(i:i+seq_length);
end

%%
tic;
% X_test = [samples, seq_length, features] in MATLAB
% pred = py.main_NN_V3.predict_model('lstm_model_internship.keras', X_test);
pred = py.main_NN_V3.predict_model('lstm_model_internship_backup_seq10.keras', X_test);
pred = double(pred);  % convert Python list to MATLAB array
timeVal = toc;
%%


err = y_val(9:end-3)-pred(1:end);
err_offsetAdj = err-err(1);
err = err;


cost_func = 'NRMSE';
fit_offsetAdj = goodnessOfFit(pred'-pred(1),y_val(9:end-3)'-y_val(9),cost_func) ;
fit = goodnessOfFit(pred',y_val(9:end-3)',cost_func) ;


figure(101);clf
    ax1 = subplot(211);grid minor;hold on
        plot(y_val(9:end)-mean(y_val),'r')
        plot(pred-mean(pred),'b')
        plot([0],'b')
        plot([0],'b')
        % xlim([1.39 1.53].*1e4)
            xlabel('Time [s]')
            ylabel('Temperature [degC]')
            title('Offset adjusted temperature esimation using a Neural Network')
            legend('Validation data','NN estimation',['   RMSE: ',num2str(rms(err_offsetAdj),3)],['   NRMSE: ',num2str(rms(fit_offsetAdj),5)])
    ax2 = subplot(212);grid minor;hold on
        plot(y_val(9:end),'r')
        plot(pred,'b')
        plot([y_val(1)],'b')
        plot([y_val(1)],'b')
        % xlim([1.39 1.53].*1e4)
            xlabel('Time [s]')
            ylabel('Temperature [degC]')
            title('Temperature esimation using a Neural Network')
            legend('Validation data','NN estimation',['   RMSE: ',num2str(rms(err),3)],['   NRMSE: ',num2str(rms(fit),5)])
    linkaxes([ax1,ax2],'x')


figure(102);clf;hold on;grid minor
    plot(err,'b')
    plot([err(1)],'b')
    plot([err(1)],'b')
    % xlim([1.39 1.53].*1e4)
        xlabel('Time [s]')
        ylabel('Estimation error [degC]')
        title('Temperature esimation error using a Neural Network')
        legend('NN estimation error',['   RMSE: ',num2str(rms(err),3)],['   NRMSE: ',num2str(rms(fit),5)])

C = xcorr(err(7500:end));
lags = linspace(-length(C)/2,length(C)/2,length(C));

figure(103);clf;
    plot(lags,C/max(C));grid minor
            xlabel('Lag [-]')
            ylabel('Autocorrelation [-]')
            title('Autocorrelation on temperature esimation error')








