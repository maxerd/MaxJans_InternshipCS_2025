
function TM_Visual(tempSimData,tVec)

% Load in some thermal and mechanical constants
c = constsGen('thermal','mech');

% The TM is divided into 9 lumps in the x direction
%                        7 lumps in the z direction
%                        3 lumps in the y direction
% This gives a total of 63 lumps per layer (layer = per y position)
Nx = 9;
Nz = 7;
Ny = 3;
Nlumps = Nx*Nz;

makeFigs = 0;

manualLumpData_water_V3

figure(700);clf;hold on
    % set(gcf,'Position',[100,200,600*1.5,350*1.5])
    xlabel('X-direction')
    ylabel('Y-direction')
        xlim([-10 270])
        ylim([-10 160])
        % Outer boundry
        rectangle('Position',[0,0,260,150],'EdgeColor','r','LineWidth',3) 
        % Lump boundries
        for i=1:63
            % clr = makeTempColor(tempSimData(j,i),max(max(tempSimData)),min(min(tempSimData)));
            h(i) = rectangle('Position',[eval(['LTM' num2str(i) '.Px']),eval(['LTM' num2str(i) '.Pz']),eval(['LTM' num2str(i) '.Lx']),eval(['LTM' num2str(i) '.Lz'])],'FaceColor',[1 1 1],'EdgeColor','c','LineWidth',1);
            plot(eval(['LTM' num2str(i) '.Px'])+eval(['LTM' num2str(i) '.Lx'])/2,eval(['LTM' num2str(i) '.Pz'])+eval(['LTM' num2str(i) '.Lz'])/2,'cx')
        end 
        % Water boundry - In & Outlet
        wc(1)  = rectangle('Position',[6,117.5,15,15], 'EdgeColor','w','LineWidth',1,'Curvature',1);
        wc(29) = rectangle('Position',[239,17.5,15,15],'EdgeColor','w','LineWidth',1,'Curvature',1);

        % ws = zeros(29,1);
        % wc = zeros(29,1);
        % Water boundry - straight
        curve = 0;
        ws(1) = rectangle('Position',[LTM6.Lx/2-4,120,LTM6.Lx/2+4,10], 'EdgeColor','w','LineWidth',1,'Curvature',curve);
        ws(2) = rectangle('Position',[24.5, 120,LTM13.Lx,10],          'EdgeColor','w','LineWidth',1,'Curvature',curve);
        ws(3) = rectangle('Position',[32.5, 120,LTM20.Lx,10],          'EdgeColor','w','LineWidth',1,'Curvature',curve);
        ws(4) = rectangle('Position',[66.5, 120,LTM27.Lx,10],          'EdgeColor','w','LineWidth',1,'Curvature',curve);
        ws(5) = rectangle('Position',[108.5,120,LTM34.Lx,10],          'EdgeColor','w','LineWidth',1,'Curvature',curve);
        ws(6) = rectangle('Position',[150.5,120,LTM41.Lx,10],          'EdgeColor','w','LineWidth',1,'Curvature',curve);
        ws(7) = rectangle('Position',[192.5,120,LTM48.Lx/2+5,10],      'EdgeColor','w','LineWidth',1,'Curvature',curve);

        ws(12) = rectangle('Position',[32.5+LTM20.Lx/2-5, 95,LTM20.Lx/2+5,10],'EdgeColor','w','LineWidth',1,'Curvature',curve);
        ws(11) = rectangle('Position',[66.5, 95,LTM27.Lx,10],          'EdgeColor','w','LineWidth',1,'Curvature',curve);
        ws(10) = rectangle('Position',[108.5,95,LTM34.Lx,10],          'EdgeColor','w','LineWidth',1,'Curvature',curve);
        ws(9) = rectangle('Position',[150.5,95,LTM41.Lx,10],           'EdgeColor','w','LineWidth',1,'Curvature',curve);
        ws(8) = rectangle('Position',[192.5,95,LTM48.Lx/2+5,10],       'EdgeColor','w','LineWidth',1,'Curvature',curve);

        ws(13) = rectangle('Position',[32.5+LTM20.Lx/2-5, 70,LTM20.Lx/2+5,10],'EdgeColor','w','LineWidth',1,'Curvature',curve);
        ws(14) = rectangle('Position',[66.5, 70,LTM27.Lx,10],          'EdgeColor','w','LineWidth',1,'Curvature',curve);
        ws(15) = rectangle('Position',[108.5,70,LTM34.Lx,10],          'EdgeColor','w','LineWidth',1,'Curvature',curve);
        ws(16) = rectangle('Position',[150.5,70,LTM41.Lx,10],          'EdgeColor','w','LineWidth',1,'Curvature',curve);
        ws(17) = rectangle('Position',[192.5,70,LTM48.Lx/2+5,10],      'EdgeColor','w','LineWidth',1,'Curvature',curve);

        ws(22) = rectangle('Position',[32.5+LTM20.Lx/2-5, 45,LTM20.Lx/2+5,10],'EdgeColor','w','LineWidth',1,'Curvature',curve);
        ws(21) = rectangle('Position',[66.5, 45,LTM27.Lx,10],          'EdgeColor','w','LineWidth',1,'Curvature',curve);
        ws(20) = rectangle('Position',[108.5,45,LTM34.Lx,10],          'EdgeColor','w','LineWidth',1,'Curvature',curve);
        ws(19) = rectangle('Position',[150.5,45,LTM41.Lx,10],          'EdgeColor','w','LineWidth',1,'Curvature',curve);
        ws(18) = rectangle('Position',[192.5,45,LTM48.Lx/2+5,10],      'EdgeColor','w','LineWidth',1,'Curvature',curve);

        ws(23) = rectangle('Position',[32.5+LTM20.Lx/2-5, 20,LTM20.Lx/2+5,10],'EdgeColor','w','LineWidth',1,'Curvature',curve);
        ws(24) = rectangle('Position',[66.5, 20,LTM27.Lx,10],          'EdgeColor','w','LineWidth',1,'Curvature',curve);
        ws(25) = rectangle('Position',[108.5,20,LTM34.Lx,10],          'EdgeColor','w','LineWidth',1,'Curvature',curve);
        ws(26) = rectangle('Position',[150.5,20,LTM41.Lx,10],          'EdgeColor','w','LineWidth',1,'Curvature',curve);
        ws(27) = rectangle('Position',[192.5,20,LTM48.Lx,10],          'EdgeColor','w','LineWidth',1,'Curvature',curve);
        ws(28) = rectangle('Position',[226.5,20,LTM55.Lx,10],          'EdgeColor','w','LineWidth',1,'Curvature',curve);
        ws(29) = rectangle('Position',[234.5,20,LTM62.Lx/2+4,10],      'EdgeColor','w','LineWidth',1,'Curvature',curve);

        % Water boundry - corners
        wc(7) = rectangle('Position',[205,95,10,35/2],    'EdgeColor','w','LineWidth',1,'Curvature',curve);
        wc(8) = rectangle('Position',[205,112.5,10,35/2], 'EdgeColor','w','LineWidth',1,'Curvature',curve);

        wc(12) = rectangle('Position',[45,70,10,35/2],    'EdgeColor','w','LineWidth',1,'Curvature',curve);
        wc(13) = rectangle('Position',[45,87.5,10,35/2],  'EdgeColor','w','LineWidth',1,'Curvature',curve);

        wc(17) = rectangle('Position',[205,45,10,35/2],   'EdgeColor','w','LineWidth',1,'Curvature',curve);
        wc(18) = rectangle('Position',[205,62.5,10,35/2], 'EdgeColor','w','LineWidth',1,'Curvature',curve);

        wc(22) = rectangle('Position',[45,20,10,35/2],    'EdgeColor','w','LineWidth',1,'Curvature',curve);
        wc(23) = rectangle('Position',[45,37.5,10,35/2],  'EdgeColor','w','LineWidth',1,'Curvature',curve);

%% Colormap
% cmap = [summer(256);flip(autumn(256))];
n = 512;   % number of colors

% Define control colors (RGB triplets in [0,1])
blue   = [0 0.5 0.8];
yellow = [1 0 0];    
red    = [1 0.8 0];  

% Split into two segments: blue→red and red→yellow
n1 = round(n/2);
n2 = n - n1;

% Interpolate blue → red
cmap1 = [linspace(blue(1), red(1), n1)', ...
         linspace(blue(2), red(2), n1)', ...
         linspace(blue(3), red(3), n1)'];

% Interpolate red → yellow
cmap2 = [linspace(red(1), yellow(1), n2)', ...
         linspace(red(2), yellow(2), n2)', ...
         linspace(red(3), yellow(3), n2)'];

% Concatenate into final colormap
cmap = [cmap1; cmap2];


%%Colormap end
%%

colorbar()
colormap(cmap)
c = colorbar;
c.Label.String = 'Temperature (°C)';
c.Ticks = linspace(0,1,6);
c.TickLabels = string(round(linspace(min(min(tempSimData(:,1:63)))-0.1,max(max(tempSimData(:,1:63)))+0.1,6),2));


% vidTitle = 'testVid';
% myVideo = VideoWriter(vidTitle,'MPEG-4'); %open video file
% myVideo.FrameRate = 10;  %can adjust this, 5 - 10 works well for me
% open(myVideo)

% for i=1:25:length(tempSimData(:,1))
for i=1:10:length(tempSimData(:,1))
    for j=1:63
        tempVal = makeTempColor(tempSimData(i,j),max(max(tempSimData(:,1:63)))+0.1,min(min(tempSimData(:,1:63)))-0.1);
        h(j).FaceColor = cmap(tempVal,:);
    end
    for k=1:29
        tempVal = makeTempColor(tempSimData(i,end-29+k),max(max(tempSimData(:,1:63)))+0.1,min(min(tempSimData(:,1:63)))-0.1);
        ws(k).FaceColor = cmap(tempVal,:);
    end
    for l=[1 7 8 12 13 17 18 22 23 29]
        tempVal = makeTempColor(tempSimData(i,end-29+l),max(max(tempSimData(:,1:63)))+0.1,min(min(tempSimData(:,1:63)))-0.1);
        wc(l).FaceColor = cmap(tempVal,:);
    end
    title(['Thermal mass lumped representation at t=',num2str(tVec(i)),'s'])
    drawnow;
    % frame = getframe(gcf); %get frame
    % writeVideo(myVideo, frame);
end
% close(myVideo)

end








