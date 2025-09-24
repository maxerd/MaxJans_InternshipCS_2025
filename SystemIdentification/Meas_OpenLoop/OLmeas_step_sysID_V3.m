function H_ident = step_sysID(datIn,datAmb,datOut,tVec,nSettle,varargin)


dc = mean(datOut(end-nSettle:end));
tau_idx = find((datOut-datOut(1))>(0.623*(dc-datOut(1))),1,'first');
tau = tVec(tau_idx);

A = (dc-datOut(1))/(mean(datIn(end-nSettle:end)));

H_ident = tf(A,[tau 1]);
H_ident = ss(H_ident);
% H_ident = ss(H_ident.A,H_ident.B/(1/H_ident.C),H_ident.C*1/H_ident.C,H_ident.D);
H_ident = compreal(H_ident,'o');

H_amb = H_ident;
H_amb.B = [H_ident.B H_ident.B/A];
H_amb = compreal(H_amb,'o');

% x0_H = datOut(1)-datOut(1);
x0_H = datOut(1).*ones(1,1);

% size(tVec)
% y_H   = lsim(H_ident,datIn,tVec,x0_H);
y_H   = lsim(H_amb,[datIn datAmb]',tVec,x0_H);

% figure(901);clf;hold on
%     subplot(211);hold on;grid minor
%         plot(tVec,datOut,'--',lineWidth=1.5)
%         plot(tVec,y_H,lineWidth=1.5)
%             xlabel('Time [s]')
%             ylabel('Temperature [degC]')
%             title('Temperature of input data and simulation using identified model')
%             legend('Input data','Simulation')
%     subplot(212);hold on;grid minor
%         plot(tVec,datOut-y_H)
%             xlabel('Time [s]')
%             ylabel('Temperature difference [degC]')
%             title('Difference between input data and simulation using identified model')

    if (nargin-5)==4
        x0_val = varargin{3}(1);
        y_val   = lsim(H_amb,[varargin{1} varargin{2}],varargin{4},x0_val);
        figure(902);clf;hold on;grid minor
            subplot(211);hold on;grid minor
                plot(varargin{4},varargin{3})
                plot(varargin{4},y_val)
            subplot(212);hold on;grid minor
                plot(varargin{4},varargin{3}-y_val)
        
    end

end








