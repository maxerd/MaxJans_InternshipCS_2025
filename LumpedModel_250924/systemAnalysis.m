

%% Make bode plot
if or(makeFigs==-1,makeFigs==2)

    if or(outputOpt==0,outputOpt==1)
        figure(500);clf;hold on
            set(gcf,'Position',[100,200,600,450])
            for i=1:5
                bode(sysTMC(1,i));grid minor
            end
            xlim([1e-5 1e1])
            title('Bodeplots to top middle lump')
            legend('HT1 to TM','HT2 to TM','HT3 to TM','HT4 to TM','HT15 to TM','Location','best')
    else 
        figure(500);clf;hold on
            set(gcf,'Position',[100,200,600,450])
            for i=1:5
                bode(sysTMC(32,i));grid minor
            end
            xlim([1e-5 1e1])
            title('Bodeplots to top middle lump')
            legend('HT1 to TM','HT2 to TM','HT3 to TM','HT4 to TM','HT15 to TM','Location','best')
    end
end

%% Observability and controllability
obs = Rank(obsv(sysTMC))==size(sysTMC.A,2);
ctl = Rank(ctrb(sysTMC))==size(sysTMC.A,2);

if obs
    disp('System is observable')
end
if ctl
    disp('System is controllable')
end
disp(' ')











