
%% Part 6

%% Direct method, PBRS

% Use signal length of 511, with 5 periods to get 2555 samples
%   and use 455 samples of sweep for validation
N = 511; % [-]
periods = 5;
N_sweep = 3000-N*periods;

% Generate time vectors, mainly for plotting purposes
t_prbs = 0:N*periods-1;
t_prbs_single = 0:N-1;
t_swep = 0:N_sweep-1;

% Generate the reference signals, 
Range = [Mlower,Mupper];
Band = [0 1];
r_p6_prbs = idinput([N,1,periods],'prbs',Band,Range);
r_p6_prbs_single = idinput([N,1,1],'prbs',Band,Range);
r_p6_swep = chirp(t_swep,0,t_swep(end),0.5*x).*2.2; % Make the sweep over the entire sampled frequency range (up to nyquist frequency)


% Put reference signals into the system to get u and y
[u_p6_prbs,y_p6_prbs] = assignment_sys_04(r_p6_prbs,'closed loop');
[u_p6_swep,y_p6_swep] = assignment_sys_04(r_p6_swep,'closed loop');


% Filter the input data, as u is needed somewhere (?) as input, and not r
u_p6_prbs_filt = lsim(buttWorth,r_p6_prbs,t_prbs);
u_p6_prbs_single_filt = lsim(buttWorth,r_p6_prbs_single,t_prbs_single);
u_p6_swep_filt = lsim(buttWorth,r_p6_swep,t_swep);

%%

% dataCL_direct = iddata(y_p6_prbs,u_p6_prbs,1/fs);

% Do a parametric identification using the (same as in part 4) model
% structure
nb = 1;
nc = 0;
nd = 0;
nf = 3;
nk = 0;
est = bj(sysID.data.CL.train.id, [nb nc nd nf nk]);

% Make the transfer function with the identified parameters
G_CL_dir = tf(est.B,est.F,est.Ts,'variable','q^-1');
H_CL_dir = tf(est.C,est.D,est.Ts,'variable','q^-1');

%% Verification of the data
figure(611);clf;hold on;
subplot(121)
    bode(est,G(1,1));grid minor
        title('Bode plots from different identification methods');
        legend('Closed loop, direct, parametric','Open loop, parametric','Open loop, non-parametric','Location','northwest');
subplot(122)
    bode(H_CL_dir);grid minor
        title('Bode plots from different identification methods');
        legend('Closed loop, direct, parametric','Open loop, parametric','Open loop, non-parametric','Location','northwest');

figure(1611);
    compare(sysID.data.CL.train.id, est)
figure(1612);
    resid(sysID.data.CL.train.id, est)

%% for indirect method prep

% Make the data for identification from r to y (input to the loop to output of the loop)
%   This results in the process sensitivity of the system
dataCL_direct_yr = iddata(y_p6_prbs,r_p6_prbs,1/fs);


% Perform a parametric identification with an OE model structure, since for
%   the indirect methods the noise model do not matter anymore therefore they
%   can be neglected in the model structure as well

% est_yr = oe(dataCL_direct_yr, [nb nf 0]);
est_yr = oe(dataCL_direct_yr, [6 10 0]);
G_yr = tf(est_yr.B,est_yr.F,est_yr.Ts,'variable','q^-1');

% Verification using the residial test, only look at the cross correlation
figure(1613);
    resid(dataCL_direct_yr, est_yr)

%% for indirect method prep

% Make the data for identification from r to u (input to the loop to input of the plant)
%   This results in the input sensitivity of the system
dataCL_direct_ur = iddata(u_p6_prbs,r_p6_prbs,1/fs);


% Perform a parametric identification with an OE model structure, since for
%   the indirect methods the noise model do not matter anymore therefore they
%   can be neglected in the model structure as well

% est_ur = oe(dataCL_direct_ur, [nb nf 0]);
est_ur = oe(dataCL_direct_ur, [12 16 0]);

G_ur = tf(est_ur.B,est_ur.F,est_ur.Ts,'variable','q^-1');

% Verification using the residial test, only look at the cross correlation
figure(1614);
    resid(dataCL_direct_ur, est_ur)

%% Indirect approach 1 - co-prime

% Since G_yr is the process sensitivity and G_ur is the input sensitivity;
%   the division between those two results in the plant model, this is the
%   essence of the co-prime method
G_CL_cp = G_yr/G_ur;

% Verify with previously obtainted model
figure(1615);
    bodemag(G_CL_cp,G_OL_nonParam);
    legend('Co-prime, closed loop','Non-param, open loop')

%% Indirect approach 2 - two-stage

% Simulate the plant input (instead of using the measurement) using the 
%   input sensitivity as defined before
u_r = lsim(G_ur,r_p6_prbs,t_prbs);

% Use the simulated plant input data for the identification in combination
%   for the measured output, known as the two-stage method
dataCL_direct_2s = iddata(y_p6_prbs,u_r,1/fs);

% Only a OE model structure should be enough since the noise model is not
%   needed for the indirect identification methods
est_2s = oe(dataCL_direct_2s, [nb nf nk]);

G_CL_2s = tf(est_2s.B,est_2s.F,est_2s.Ts);

% Verify with previously obtainted model
figure(616);clf;hold on;
    bode(G_CL_2s,G_OL_param,opt)

%% Plot comparing all the different indirect parametric identification methods, including open loop

figure(1617);clf;hold on;
    bode(G_OL_param,'b',opt);
    bode(G_CL_dir,'r--',opt);
    bode(G_CL_2s,'g-.',opt);
    bode(G_CL_cp,'m:',opt);grid minor
        xlim([0.01 10]);
        title('Bode plots of different parametric identifications')
        legend('Open loop','Closed loop, direct','Closed loop, indirect, two-stage','Closed loop, indirect, co-prime','Location','northwest')

%% Comparing the open loop and two stage indentification methods
figure(618);clf;hold on;
    bode(G_OL_param,'b',opt);
    bode(G_CL_2s,'r-.',opt);grid minor
        xlim([0.01 10]);
        title('Bode plots of open and closed loop parametric identifications')
        legend('Open loop','Closed loop (indirect, two-stage)','Location','northwest')


%% Comparing the open loop and two stage indentification methods
figure(619);clf;hold on;
    bode(G_OL_nonParam,'b',opt);
    bode(G_OL_param,'r',opt);
    bode(G_CL_2s,'g-.',opt);grid minor
        xlim([0.01 10]);
        ylim([])
        title('Bode plots of open and closed loop parametric identifications')
        legend('Open loop, non-parametric','Open loop, parametric','Closed loop (indirect, two-stage)','Location','northwest')
    % print('figures\part61_comp_OL_CL_nonPar_Par','-depsc')
























