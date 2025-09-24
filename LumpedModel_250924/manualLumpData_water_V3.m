

% This file contains all the data that is aquired on the lumps
% All these values have to manually be found and put into the file in case
% anything changes
%   Future work might include automation of this process


%% Lump masses
% The mass of all the seperate lumps, with symmetry included

TM1.mass  = 9.95;  % Only A
TM2.mass  = 18.19; % Only A
TM3.mass  = 20.25; % Only A
TM4.mass  = 19.91; % Only A
TM5.mass  = 23.92; % Not Full A&B
TM6.mass  = 32.61;
TM7.mass  = 18.53;
TM8.mass  = 11.00; % Fill in
TM9.mass  = 11.00; % Fill in
TM10.mass = 11.00; % Fill in
TM11.mass = 11.00; % Fill in
TM12.mass = 11.00; % Fill in
TM13.mass = 11.00; % Fill in
TM14.mass = 11.00; % Fill in
TM15.mass = 26.16;
TM16.mass = 48.94;
TM17.mass = 47.47;
TM18.mass = 46.91;
TM19.mass = 47.52;
TM20.mass = 47.36;
TM21.mass = 26.05;
TM22.mass = 31.82;
TM23.mass = 58.63;
TM24.mass = 58.50;
TM25.mass = 57.74;
TM26.mass = 58.38;
TM27.mass = 58.38;
TM28.mass = 32.48;
TM29.mass = 31.88;
TM30.mass = 58.63;
TM31.mass = 58.50;
TM32.mass = 57.87;
% Symmetry
TM33.mass = TM31.mass;
TM34.mass = TM30.mass;
TM35.mass = TM29.mass;
TM36.mass = TM28.mass;
TM37.mass = TM27.mass;
TM38.mass = TM26.mass;
TM39.mass = TM25.mass;
TM40.mass = TM24.mass;
TM41.mass = TM23.mass;
TM42.mass = TM22.mass;
TM43.mass = TM21.mass;
TM44.mass = TM20.mass;
TM45.mass = TM19.mass;
TM46.mass = TM18.mass;
TM47.mass = TM17.mass;
TM48.mass = TM16.mass;
TM49.mass = TM15.mass;
TM50.mass = TM14.mass;
TM51.mass = TM13.mass;
TM52.mass = TM12.mass;
TM53.mass = TM11.mass;
TM54.mass = TM10.mass;
TM55.mass = TM9.mass;
TM56.mass = TM8.mass;
TM57.mass = TM7.mass;
TM58.mass = TM6.mass;
TM59.mass = TM5.mass;
TM60.mass = TM4.mass;
TM61.mass = TM3.mass;
TM62.mass = TM2.mass;
TM63.mass = TM1.mass;

for i=1:Nlumps
    TMA(i).mass = eval(['TM' num2str(i) '.mass'])/2;
    TMB(i).mass = eval(['TM' num2str(i) '.mass'])/2;
end

% Very large mass for ambient air, assuming the ambient air does not heat
% up
% This value is randomly chosen, but that does not matter as long as it is
% large enough
LAMB.mass = 1e15;

%% Lump dimensions
% [mm]

LTM1.Lx   = 25.0;    LTM1.Lz   = 12.5;    LTM1.Ly   = 12.0;
LTM2.Lx   = 25.0;    LTM2.Lz   = 25.0;    LTM2.Ly   = 12.0;
LTM3.Lx   = 25.0;    LTM3.Lz   = 25.0;    LTM3.Ly   = 12.0;
LTM4.Lx   = 25.0;    LTM4.Lz   = 25.0;    LTM4.Ly   = 12.0;
LTM5.Lx   = 25.0;    LTM5.Lz   = 25.0;    LTM5.Ly   = -1.0;
LTM6.Lx   = 25.0;    LTM6.Lz   = 25.0;    LTM6.Ly   = 24.0;
LTM7.Lx   = 25.0;    LTM7.Lz   = 12.5;    LTM7.Ly   = 24.0;
  
LTM8.Lx   = 8.0;     LTM8.Lz   = 12.5;    LTM8.Ly   = 24.0;
LTM9.Lx   = 8.0;     LTM9.Lz   = 25.0;    LTM9.Ly   = 24.0;
LTM10.Lx  = 8.0;     LTM10.Lz  = 25.0;    LTM10.Ly  = 24.0;
LTM11.Lx  = 8.0;     LTM11.Lz  = 25.0;    LTM11.Ly  = 24.0;
LTM12.Lx  = 8.0;     LTM12.Lz  = 25.0;    LTM12.Ly  = 24.0;
LTM13.Lx  = 8.0;     LTM13.Lz  = 25.0;    LTM13.Ly  = 24.0;
LTM14.Lx  = 8.0;     LTM14.Lz  = 12.5;    LTM14.Ly  = 24.0;
  
LTM15.Lx  = 34.0;    LTM15.Lz  = 12.5;    LTM15.Ly  = 24.0;
LTM16.Lx  = 34.0;    LTM16.Lz  = 25.0;    LTM16.Ly  = 24.0;
LTM17.Lx  = 34.0;    LTM17.Lz  = 25.0;    LTM17.Ly  = 24.0;
LTM18.Lx  = 34.0;    LTM18.Lz  = 25.0;    LTM18.Ly  = 24.0;
LTM19.Lx  = 34.0;    LTM19.Lz  = 25.0;    LTM19.Ly  = 24.0;
LTM20.Lx  = 34.0;    LTM20.Lz  = 25.0;    LTM20.Ly  = 24.0;
LTM21.Lx  = 34.0;    LTM21.Lz  = 12.5;    LTM21.Ly  = 24.0;
  
LTM22.Lx  = 42.0;    LTM22.Lz  = 12.5;    LTM22.Ly  = 24.0;
LTM23.Lx  = 42.0;    LTM23.Lz  = 25.0;    LTM23.Ly  = 24.0;
LTM24.Lx  = 42.0;    LTM24.Lz  = 25.0;    LTM24.Ly  = 24.0;
LTM25.Lx  = 42.0;    LTM25.Lz  = 25.0;    LTM25.Ly  = 24.0;
LTM26.Lx  = 42.0;    LTM26.Lz  = 25.0;    LTM26.Ly  = 24.0;
LTM27.Lx  = 42.0;    LTM27.Lz  = 25.0;    LTM27.Ly  = 24.0;
LTM28.Lx  = 42.0;    LTM28.Lz  = 12.5;    LTM28.Ly  = 24.0;
  
LTM29.Lx  = 42.0;    LTM29.Lz  = 12.5;    LTM29.Ly  = 24.0;
LTM30.Lx  = 42.0;    LTM30.Lz  = 25.0;    LTM30.Ly  = 24.0;
LTM31.Lx  = 42.0;    LTM31.Lz  = 25.0;    LTM31.Ly  = 24.0;
LTM32.Lx  = 42.0;    LTM32.Lz  = 25.0;    LTM32.Ly  = 24.0;
LTM33.Lx  = 42.0;    LTM33.Lz  = 25.0;    LTM33.Ly  = 24.0;
LTM34.Lx  = 42.0;    LTM34.Lz  = 25.0;    LTM34.Ly  = 24.0;
LTM35.Lx  = 42.0;    LTM35.Lz  = 12.5;    LTM35.Ly  = 24.0;
  
LTM36.Lx  = 42.0;    LTM36.Lz  = 12.5;    LTM36.Ly  = 24.0;
LTM37.Lx  = 42.0;    LTM37.Lz  = 25.0;    LTM37.Ly  = 24.0;
LTM38.Lx  = 42.0;    LTM38.Lz  = 25.0;    LTM38.Ly  = 24.0;
LTM39.Lx  = 42.0;    LTM39.Lz  = 25.0;    LTM39.Ly  = 24.0;
LTM40.Lx  = 42.0;    LTM40.Lz  = 25.0;    LTM40.Ly  = 24.0;
LTM41.Lx  = 42.0;    LTM41.Lz  = 25.0;    LTM41.Ly  = 24.0;
LTM42.Lx  = 42.0;    LTM42.Lz  = 12.5;    LTM42.Ly  = 24.0;
  
LTM43.Lx  = 34.0;    LTM43.Lz  = 12.5;    LTM43.Ly  = 24.0;
LTM44.Lx  = 34.0;    LTM44.Lz  = 25.0;    LTM44.Ly  = 24.0;
LTM45.Lx  = 34.0;    LTM45.Lz  = 25.0;    LTM45.Ly  = 24.0;
LTM46.Lx  = 34.0;    LTM46.Lz  = 25.0;    LTM46.Ly  = 24.0;
LTM47.Lx  = 34.0;    LTM47.Lz  = 25.0;    LTM47.Ly  = 24.0;
LTM48.Lx  = 34.0;    LTM48.Lz  = 25.0;    LTM48.Ly  = 24.0;
LTM49.Lx  = 34.0;    LTM49.Lz  = 12.5;    LTM49.Ly  = 24.0;
  
LTM50.Lx  = 8.0;     LTM50.Lz  = 12.5;    LTM50.Ly  = 24.0;
LTM51.Lx  = 8.0;     LTM51.Lz  = 25.0;    LTM51.Ly  = 24.0;
LTM52.Lx  = 8.0;     LTM52.Lz  = 25.0;    LTM52.Ly  = 24.0;
LTM53.Lx  = 8.0;     LTM53.Lz  = 25.0;    LTM53.Ly  = 24.0;
LTM54.Lx  = 8.0;     LTM54.Lz  = 25.0;    LTM54.Ly  = 24.0;
LTM55.Lx  = 8.0;     LTM55.Lz  = 25.0;    LTM55.Ly  = 24.0;
LTM56.Lx  = 8.0;     LTM56.Lz  = 12.5;    LTM56.Ly  = 24.0;
  
LTM57.Lx  = 25.0;    LTM57.Lz  = 12.5;    LTM57.Ly  = 24.0;
LTM58.Lx  = 25.0;    LTM58.Lz  = 25.0;    LTM58.Ly  = 24.0;
LTM59.Lx  = 25.0;    LTM59.Lz  = 25.0;    LTM59.Ly  = -1.0;
LTM60.Lx  = 25.0;    LTM60.Lz  = 25.0;    LTM60.Ly  = 12.0;
LTM61.Lx  = 25.0;    LTM61.Lz  = 25.0;    LTM61.Ly  = 12.0;
LTM62.Lx  = 25.0;    LTM62.Lz  = 25.0;    LTM62.Ly  = 12.0;
LTM63.Lx  = 25.0;    LTM63.Lz  = 12.5;    LTM63.Ly  = 12.0;

%% Define Lump position (left bottom corner)
% [mm]

LTM1.Px   = 0;                   LTM1.Pz   = 0;                   LTM1.Py   = 0;
LTM2.Px   = 0;                   LTM2.Pz   = LTM1.Lz;             LTM2.Py   = 0;
LTM3.Px   = 0;                   LTM3.Pz   = LTM2.Pz+LTM2.Lz;     LTM3.Py   = 0;
LTM4.Px   = 0;                   LTM4.Pz   = LTM3.Pz+LTM3.Lz;     LTM4.Py   = 0;
LTM5.Px   = 0;                   LTM5.Pz   = LTM4.Pz+LTM4.Lz;     LTM5.Py   = 0;
LTM6.Px   = 0;                   LTM6.Pz   = LTM5.Pz+LTM5.Lz;     LTM6.Py   = 0;
LTM7.Px   = 0;                   LTM7.Pz   = LTM6.Pz+LTM6.Lz;     LTM7.Py   = 0;

LTM8.Px   = LTM1.Lx;             LTM8.Pz   = 0;                   LTM8.Py   = 0;
LTM9.Px   = LTM1.Lx;             LTM9.Pz   = LTM1.Lz;             LTM9.Py   = 0;
LTM10.Px  = LTM1.Lx;             LTM10.Pz  = LTM9.Pz+LTM2.Lz;     LTM10.Py  = 0;
LTM11.Px  = LTM1.Lx;             LTM11.Pz  = LTM10.Pz+LTM3.Lz;    LTM11.Py  = 0;
LTM12.Px  = LTM1.Lx;             LTM12.Pz  = LTM11.Pz+LTM4.Lz;    LTM12.Py  = 0;
LTM13.Px  = LTM1.Lx;             LTM13.Pz  = LTM12.Pz+LTM5.Lz;    LTM13.Py  = 0;
LTM14.Px  = LTM1.Lx;             LTM14.Pz  = LTM13.Pz+LTM6.Lz;    LTM14.Py  = 0;

LTM15.Px  = LTM8.Px +LTM8.Lx;    LTM15.Pz  = 0;                   LTM15.Py  = 0;
LTM16.Px  = LTM9.Px +LTM8.Lx;    LTM16.Pz  = LTM1.Lz;             LTM16.Py  = 0;
LTM17.Px  = LTM10.Px+LTM8.Lx;    LTM17.Pz  = LTM9.Pz+LTM2.Lz;     LTM17.Py  = 0;
LTM18.Px  = LTM11.Px+LTM8.Lx;    LTM18.Pz  = LTM10.Pz+LTM3.Lz;    LTM18.Py  = 0;
LTM19.Px  = LTM12.Px+LTM8.Lx;    LTM19.Pz  = LTM11.Pz+LTM4.Lz;    LTM19.Py  = 0;
LTM20.Px  = LTM13.Px+LTM8.Lx;    LTM20.Pz  = LTM12.Pz+LTM5.Lz;    LTM20.Py  = 0;
LTM21.Px  = LTM14.Px+LTM8.Lx;    LTM21.Pz  = LTM13.Pz+LTM6.Lz;    LTM21.Py  = 0;

LTM22.Px  = LTM15.Px+LTM15.Lx;   LTM22.Pz  = 0;                   LTM22.Py  = 0;
LTM23.Px  = LTM16.Px+LTM15.Lx;   LTM23.Pz  = LTM1.Lz;             LTM23.Py  = 0;
LTM24.Px  = LTM17.Px+LTM15.Lx;   LTM24.Pz  = LTM9.Pz+LTM2.Lz;     LTM24.Py  = 0;
LTM25.Px  = LTM18.Px+LTM15.Lx;   LTM25.Pz  = LTM10.Pz+LTM3.Lz;    LTM25.Py  = 0;
LTM26.Px  = LTM19.Px+LTM15.Lx;   LTM26.Pz  = LTM11.Pz+LTM4.Lz;    LTM26.Py  = 0;
LTM27.Px  = LTM20.Px+LTM15.Lx;   LTM27.Pz  = LTM12.Pz+LTM5.Lz;    LTM27.Py  = 0;
LTM28.Px  = LTM21.Px+LTM15.Lx;   LTM28.Pz  = LTM13.Pz+LTM6.Lz;    LTM28.Py  = 0;

LTM29.Px  = LTM22.Px+LTM22.Lx;   LTM29.Pz  = 0;                   LTM29.Py  = 0;
LTM30.Px  = LTM23.Px+LTM22.Lx;   LTM30.Pz  = LTM1.Lz;             LTM30.Py  = 0;
LTM31.Px  = LTM24.Px+LTM22.Lx;   LTM31.Pz  = LTM9.Pz+LTM2.Lz;     LTM31.Py  = 0;
LTM32.Px  = LTM25.Px+LTM22.Lx;   LTM32.Pz  = LTM10.Pz+LTM3.Lz;    LTM32.Py  = 0;
LTM33.Px  = LTM26.Px+LTM22.Lx;   LTM33.Pz  = LTM11.Pz+LTM4.Lz;    LTM33.Py  = 0;
LTM34.Px  = LTM27.Px+LTM22.Lx;   LTM34.Pz  = LTM12.Pz+LTM5.Lz;    LTM34.Py  = 0;
LTM35.Px  = LTM28.Px+LTM22.Lx;   LTM35.Pz  = LTM13.Pz+LTM6.Lz;    LTM35.Py  = 0;

LTM36.Px  = LTM29.Px+LTM29.Lx;   LTM36.Pz  = 0;                   LTM36.Py  = 0;
LTM37.Px  = LTM30.Px+LTM30.Lx;   LTM37.Pz  = LTM1.Lz;             LTM37.Py  = 0;
LTM38.Px  = LTM31.Px+LTM31.Lx;   LTM38.Pz  = LTM9.Pz+LTM2.Lz;     LTM38.Py  = 0;
LTM39.Px  = LTM32.Px+LTM32.Lx;   LTM39.Pz  = LTM10.Pz+LTM3.Lz;    LTM39.Py  = 0;
LTM40.Px  = LTM33.Px+LTM33.Lx;   LTM40.Pz  = LTM11.Pz+LTM4.Lz;    LTM40.Py  = 0;
LTM41.Px  = LTM34.Px+LTM34.Lx;   LTM41.Pz  = LTM12.Pz+LTM5.Lz;    LTM41.Py  = 0;
LTM42.Px  = LTM35.Px+LTM35.Lx;   LTM42.Pz  = LTM13.Pz+LTM6.Lz;    LTM42.Py  = 0;

LTM43.Px  = LTM36.Px+LTM36.Lx;   LTM43.Pz  = 0;                   LTM43.Py  = 0;
LTM44.Px  = LTM37.Px+LTM37.Lx;   LTM44.Pz  = LTM1.Lz;             LTM44.Py  = 0;
LTM45.Px  = LTM38.Px+LTM38.Lx;   LTM45.Pz  = LTM9.Pz+LTM2.Lz;     LTM45.Py  = 0;
LTM46.Px  = LTM39.Px+LTM39.Lx;   LTM46.Pz  = LTM10.Pz+LTM3.Lz;    LTM46.Py  = 0;
LTM47.Px  = LTM40.Px+LTM40.Lx;   LTM47.Pz  = LTM11.Pz+LTM4.Lz;    LTM47.Py  = 0;
LTM48.Px  = LTM41.Px+LTM41.Lx;   LTM48.Pz  = LTM12.Pz+LTM5.Lz;    LTM48.Py  = 0;
LTM49.Px  = LTM42.Px+LTM42.Lx;   LTM49.Pz  = LTM13.Pz+LTM6.Lz;    LTM49.Py  = 0;

LTM50.Px  = LTM43.Px+LTM43.Lx;   LTM50.Pz  = 0;                   LTM50.Py  = 0;
LTM51.Px  = LTM44.Px+LTM44.Lx;   LTM51.Pz  = LTM1.Lz;             LTM51.Py  = 0;
LTM52.Px  = LTM45.Px+LTM45.Lx;   LTM52.Pz  = LTM9.Pz+LTM2.Lz;     LTM52.Py  = 0;
LTM53.Px  = LTM46.Px+LTM46.Lx;   LTM53.Pz  = LTM10.Pz+LTM3.Lz;    LTM53.Py  = 0;
LTM54.Px  = LTM47.Px+LTM47.Lx;   LTM54.Pz  = LTM11.Pz+LTM4.Lz;    LTM54.Py  = 0;
LTM55.Px  = LTM48.Px+LTM48.Lx;   LTM55.Pz  = LTM12.Pz+LTM5.Lz;    LTM55.Py  = 0;
LTM56.Px  = LTM49.Px+LTM49.Lx;   LTM56.Pz  = LTM13.Pz+LTM6.Lz;    LTM56.Py  = 0;

LTM57.Px  = LTM50.Px+LTM50.Lx;   LTM57.Pz  = 0;                   LTM57.Py  = 0;
LTM58.Px  = LTM51.Px+LTM51.Lx;   LTM58.Pz  = LTM1.Lz;             LTM58.Py  = 0;
LTM59.Px  = LTM52.Px+LTM52.Lx;   LTM59.Pz  = LTM9.Pz+LTM2.Lz;     LTM59.Py  = 0;
LTM60.Px  = LTM53.Px+LTM53.Lx;   LTM60.Pz  = LTM10.Pz+LTM3.Lz;    LTM60.Py  = 0;
LTM61.Px  = LTM54.Px+LTM54.Lx;   LTM61.Pz  = LTM11.Pz+LTM4.Lz;    LTM61.Py  = 0;
LTM62.Px  = LTM55.Px+LTM55.Lx;   LTM62.Pz  = LTM12.Pz+LTM5.Lz;    LTM62.Py  = 0;
LTM63.Px  = LTM56.Px+LTM56.Lx;   LTM63.Pz  = LTM13.Pz+LTM6.Lz;    LTM63.Py  = 0;

%% Make a point cloud from the lump positions
PxVec = [LTM1.Px];
PzVec = [LTM1.Pz];
for i=2:35
    PxVec = [PxVec;eval(['LTM' num2str(i) '.Px'])];
    PzVec = [PzVec;eval(['LTM' num2str(i) '.Pz'])];
end

%% Visualize the lumps, including the water-path placement
if or(makeFigs==-1,nonzeros(makeFigs==1))

figure(700);clf;hold on
    set(gcf,'Position',[100,200,600,350])
    xlabel('X-direction')
    ylabel('Y-direction')
    title('Thermal mass lumped representation')
        xlim([-10 270])
        ylim([-10 160])
        % Outer boundry
        rectangle('Position',[0,0,260,150],'EdgeColor','r','LineWidth',3) 
        % Lump boundries
        for i=1:63
            rectangle('Position',[eval(['LTM' num2str(i) '.Px']),eval(['LTM' num2str(i) '.Pz']),eval(['LTM' num2str(i) '.Lx']),eval(['LTM' num2str(i) '.Lz'])],'EdgeColor','c','LineWidth',1)
            plot(eval(['LTM' num2str(i) '.Px'])+eval(['LTM' num2str(i) '.Lx'])/2,eval(['LTM' num2str(i) '.Pz'])+eval(['LTM' num2str(i) '.Lz'])/2,'bx')
        end 
        % Water boundry - In & Outlet
        rectangle('Position',[6,117.5,15,15],'EdgeColor','m','LineWidth',1,'Curvature',1)
        rectangle('Position',[239,17.5,15,15],'EdgeColor','m','LineWidth',1,'Curvature',1)
        % Water boundry - straight
        curve = 1;
        rectangle('Position',[8.5,120,206.5,10],'EdgeColor','m','LineWidth',1,'Curvature',curve)
        rectangle('Position',[45,95,170,10],'EdgeColor','m','LineWidth',1,'Curvature',curve)
        rectangle('Position',[45,70,170,10],'EdgeColor','m','LineWidth',1,'Curvature',curve)
        rectangle('Position',[45,45,170,10],'EdgeColor','m','LineWidth',1,'Curvature',curve)
        rectangle('Position',[45,20,206.5,10],'EdgeColor','m','LineWidth',1,'Curvature',curve)
        % Water boundry - corners
        rectangle('Position',[45,70,10,35],'EdgeColor','m','LineWidth',1,'Curvature',curve)
        rectangle('Position',[45,20,10,35],'EdgeColor','m','LineWidth',1,'Curvature',curve)
        rectangle('Position',[205,45,10,35],'EdgeColor','m','LineWidth',1,'Curvature',curve)
        rectangle('Position',[205,95,10,35],'EdgeColor','m','LineWidth',1,'Curvature',curve)
end

%% Define Conductive conductance - Top to bottom
LTM1.GTB  = 5; % Fill in
LTM2.GTB  = 5; % Fill in
LTM3.GTB  = 5; % Fill in
LTM4.GTB  = 5; % Fill in
LTM5.GTB  = 5; % Fill in
LTM6.GTB  = 5; % Fill in
LTM7.GTB  = 5; % Fill in
LTM8.GTB  = 5; % Fill in
LTM9.GTB  = 5; % Fill in
LTM10.GTB = 5; % Fill in
LTM11.GTB = 5; % Fill in
LTM12.GTB = 5; % Fill in
LTM13.GTB = 5; % Fill in
LTM14.GTB = 5; % Fill in
LTM15.GTB = 2.7654;
LTM16.GTB = 4.3814;
LTM17.GTB = 4.4854;
LTM18.GTB = 4.4194;
LTM19.GTB = 4.4767;
LTM20.GTB = 4.1905;
LTM21.GTB = 2.7498;
LTM22.GTB = 3.3493;
LTM23.GTB = 5.3308;
LTM24.GTB = 5.6555;
LTM25.GTB = 5.5668;
LTM26.GTB = 5.6360;
LTM27.GTB = 5.3432;
LTM28.GTB = 3.4239;
LTM29.GTB = 3.3649;
LTM30.GTB = 5.3348;
LTM31.GTB = 5.6546;
LTM32.GTB = 5.5828;
% Symmetry
LTM33.GTB = LTM31.GTB;
LTM34.GTB = LTM30.GTB;
LTM35.GTB = LTM29.GTB;
LTM36.GTB = LTM28.GTB;
LTM37.GTB = LTM27.GTB;
LTM38.GTB = LTM26.GTB;
LTM39.GTB = LTM25.GTB;
LTM40.GTB = LTM24.GTB;
LTM41.GTB = LTM23.GTB;
LTM42.GTB = LTM22.GTB;
LTM43.GTB = LTM21.GTB;
LTM44.GTB = LTM20.GTB;
LTM45.GTB = LTM19.GTB;
LTM46.GTB = LTM18.GTB;
LTM47.GTB = LTM17.GTB;
LTM48.GTB = LTM16.GTB;
LTM49.GTB = LTM15.GTB;
LTM50.GTB = LTM14.GTB;
LTM51.GTB = LTM13.GTB;
LTM52.GTB = LTM12.GTB;
LTM53.GTB = LTM11.GTB;
LTM54.GTB = LTM10.GTB;
LTM55.GTB = LTM9.GTB;
LTM56.GTB = LTM8.GTB;
LTM57.GTB = LTM7.GTB;
LTM58.GTB = LTM6.GTB;
LTM59.GTB = LTM5.GTB;
LTM60.GTB = LTM4.GTB;
LTM61.GTB = LTM3.GTB;
LTM62.GTB = LTM2.GTB;
LTM63.GTB = LTM1.GTB;

% Multiply the conductance by 2, to divide the total height of the lumps by
% 2 to make two seperate lumps. The new height is exactly twice as small
for i=1:Nlumps
    TMA(i).GTB = eval(['LTM' num2str(i) '.GTB'])*2*1.18;
    TMB(i).GTB = eval(['LTM' num2str(i) '.GTB'])*2*1.18;
    % TMC(i).GTB = eval(['LTM' num2str(i) '.GTB'])/98;
    % TMC(i).GTB = eval(['LTM' num2str(i) '.GTB'])/103;
    TMC(i).GTB = eval(['LTM' num2str(i) '.GTB'])/75;
    % TMC(i).GTB = eval(['LTM' num2str(i) '.GTB'])/10;
    % TMC(i).GTB = eval(['LTM' num2str(i) '.GTB'])/150;
end


%% Define Conductive conductance - X-direction
LTM1.GX  = 2; % Fill in
LTM2.GX  = 2; % Fill in
LTM3.GX  = 2; % Fill in
LTM4.GX  = 2; % Fill in
LTM5.GX  = 2; % Fill in
LTM6.GX  = 2; % Fill in
LTM7.GX  = 2; % Fill in
LTM8.GX  = 2; % Fill in
LTM9.GX  = 2; % Fill in
LTM10.GX = 2; % Fill in
LTM11.GX = 2; % Fill in
LTM12.GX = 2; % Fill in
LTM13.GX = 2; % Fill in
LTM14.GX = 2; % Fill in
LTM15.GX = 1.3701; 
LTM16.GX = 2.5248; 
LTM17.GX = 2.4274; 
LTM18.GX = 2.3721; 
LTM19.GX = 2.4372; 
LTM20.GX = 2.5879; 
LTM21.GX = 1.3550; 
LTM22.GX = 1.0741; 
LTM23.GX = 2.0928; 
LTM24.GX = 2.0683; 
LTM25.GX = 2.0171; 
LTM26.GX = 2.0612; 
LTM27.GX = 2.0728; 
LTM28.GX = 1.1077;
LTM29.GX = 1.0711;
LTM30.GX = 2.0921;
LTM31.GX = 2.0684;
LTM32.GX = 2.0261;
% Symmetry
LTM33.GX = LTM31.GX;
LTM34.GX = LTM30.GX;
LTM35.GX = LTM29.GX;
LTM36.GX = LTM28.GX;
LTM37.GX = LTM27.GX;
LTM38.GX = LTM26.GX;
LTM39.GX = LTM25.GX;
LTM40.GX = LTM24.GX;
LTM41.GX = LTM23.GX;
LTM42.GX = LTM22.GX;
LTM43.GX = LTM21.GX;
LTM44.GX = LTM20.GX;
LTM45.GX = LTM19.GX;
LTM46.GX = LTM18.GX;
LTM47.GX = LTM17.GX;
LTM48.GX = LTM16.GX;
LTM49.GX = LTM15.GX;
LTM50.GX = LTM14.GX;
LTM51.GX = LTM13.GX;
LTM52.GX = LTM12.GX;
LTM53.GX = LTM11.GX;
LTM54.GX = LTM10.GX;
LTM55.GX = LTM9.GX;
LTM56.GX = LTM8.GX;
LTM57.GX = LTM7.GX;
LTM58.GX = LTM6.GX;
LTM59.GX = LTM5.GX;
LTM60.GX = LTM4.GX;
LTM61.GX = LTM3.GX;
LTM62.GX = LTM2.GX;
LTM63.GX = LTM1.GX;

% Divide the conductance by 2, to get the conductance to the edge of the
% lumps
for i=1:Nlumps
    TMA(i).GX = eval(['LTM' num2str(i) '.GX'])/2*1.18;
    TMB(i).GX = eval(['LTM' num2str(i) '.GX'])/2*1.18;
end

%% Conductive conductance - Z-direction
% G = (k*A)/L

LTM1.GZ  = 5; % Fill in
LTM2.GZ  = 5; % Fill in
LTM3.GZ  = 5; % Fill in
LTM4.GZ  = 5; % Fill in
LTM5.GZ  = 5; % Fill in
LTM6.GZ  = 5; % Fill in
LTM7.GZ  = 5; % Fill in
LTM8.GZ  = 5; % Fill in
LTM9.GZ  = 5; % Fill in
LTM10.GZ = 5; % Fill in
LTM11.GZ = 5; % Fill in
LTM12.GZ = 5; % Fill in
LTM13.GZ = 5; % Fill in
LTM14.GZ = 5; % Fill in
LTM15.GZ = 10.016;
LTM16.GZ = 4.6398;
LTM17.GZ = 4.4125;
LTM18.GZ = 4.2971;
LTM19.GZ = 4.3890;
LTM20.GZ = 4.2166;
LTM21.GZ = 10.066;
LTM22.GZ = 12.285;
LTM23.GZ = 5.2314;
LTM24.GZ = 5.1821;
LTM25.GZ = 5.0920;
LTM26.GZ = 5.2266;
LTM27.GZ = 5.2299;
LTM28.GZ = 12.748;
LTM29.GZ = 12.328;
LTM30.GZ = 5.2321;
LTM31.GZ = 5.1826;
LTM32.GZ = 5.1062;
% Symmetry
LTM33.GZ = LTM31.GZ;
LTM34.GZ = LTM30.GZ;
LTM35.GZ = LTM29.GZ;
LTM36.GZ = LTM28.GZ;
LTM37.GZ = LTM27.GZ;
LTM38.GZ = LTM26.GZ;
LTM39.GZ = LTM25.GZ;
LTM40.GZ = LTM24.GZ;
LTM41.GZ = LTM23.GZ;
LTM42.GZ = LTM22.GZ;
LTM43.GZ = LTM21.GZ;
LTM44.GZ = LTM20.GZ;
LTM45.GZ = LTM19.GZ;
LTM46.GZ = LTM18.GZ;
LTM47.GZ = LTM17.GZ;
LTM48.GZ = LTM16.GZ;
LTM49.GZ = LTM15.GZ;
LTM50.GZ = LTM14.GZ;
LTM51.GZ = LTM13.GZ;
LTM52.GZ = LTM12.GZ;
LTM53.GZ = LTM11.GZ;
LTM54.GZ = LTM10.GZ;
LTM55.GZ = LTM9.GZ;
LTM56.GZ = LTM8.GZ;
LTM57.GZ = LTM7.GZ;
LTM58.GZ = LTM6.GZ;
LTM59.GZ = LTM5.GZ;
LTM60.GZ = LTM4.GZ;
LTM61.GZ = LTM3.GZ;
LTM62.GZ = LTM2.GZ;
LTM63.GZ = LTM1.GZ;

% Divide the conductance by 2, to get the conductance to the edge of the
% lumps
for i=1:Nlumps
    TMA(i).GZ = eval(['LTM' num2str(i) '.GZ'])/2*1.18;
    TMB(i).GZ = eval(['LTM' num2str(i) '.GZ'])/2*1.18;
end

%% Convection

% Area = Top area + Side Area - Overlap
LTMA.A(1)  = 625.00+300;
LTMB.A(1)  = 0;
LTMA.A(2)  = 561.38+300;
LTMB.A(2)  = 0;
LTMA.A(3)  = 625.00+300;
LTMB.A(3)  = 0;
LTMA.A(4)  = 614.47+300;
LTMB.A(4)  = 0;
LTMA.A(5)  = 559.34+300-189.07; % This is an overlap, the 189.07mm^2 is that overlap
LTMB.A(5)  = 158.64+432.74;
LTMA.A(6)  = 0+300;
LTMB.A(6)  = 619.69+300;
LTMA.A(7)  = 0+424.25;
LTMB.A(7)  = 256.87+424.25;
LTMA.A(8)  = 0; % Fill in
LTMB.A(8)  = 0; % Fill in
LTMA.A(9)  = 0; % Fill in
LTMB.A(9)  = 0; % Fill in
LTMA.A(10) = 0; % Fill in
LTMB.A(10) = 0; % Fill in
LTMA.A(11) = 0; % Fill in
LTMB.A(11) = 0; % Fill in
LTMA.A(12) = 0; % Fill in
LTMB.A(12) = 0; % Fill in
LTMA.A(13) = 0; % Fill in
LTMB.A(13) = 0; % Fill in
LTMA.A(14) = 0; % Fill in
LTMB.A(14) = 0; % Fill in
LTMA.A(15) = 0+408;
LTMB.A(15) = 377.84+408;
LTMA.A(16) = 0+0;
LTMB.A(16) = 0+844.69+0;
LTMA.A(17) = 0+0;
LTMB.A(17) = 783.01+0;
LTMA.A(18) = 0+0;
LTMB.A(18) = 760.04+0;
LTMA.A(19) = 0+0;
LTMB.A(19) = 785.17+0;
LTMA.A(20) = 0+0;
LTMB.A(20) = 850+0;
LTMA.A(21) = 0+408;
LTMB.A(21) = 374.73+408;
LTMA.A(22) = 0+504;
LTMB.A(22) = 448.8+504;
LTMA.A(23) = 0+0;
LTMB.A(23) = 1050+0;
LTMA.A(24) = 0+0;
LTMB.A(24) = 1024.87+0;
LTMA.A(25) = 0+0;
LTMB.A(25) = 992.28+0;
LTMA.A(26) = 0+0;
LTMB.A(26) = 1017.41+0;
LTMA.A(27) = 0+0;
LTMB.A(27) = 1050+0;
LTMA.A(28) = 0+504;
LTMB.A(28) = 474.73+504;
LTMA.A(29) = 0+504;
LTMB.A(29) = 450.4+504;
LTMA.A(30) = 0+0;
LTMB.A(30) = 1050+0;
LTMA.A(31) = 0+0;
LTMB.A(31) = 1024.87+0;
LTMA.A(32) = 0+0;
LTMB.A(32) = 999.73+0;
% Symmetry
for i = (ceil(Nlumps/2)+1):Nlumps
    % disp([num2str(i),'==',num2str(Nlumps-i+1)]) % For debugging
    LTMA.A(i) = LTMA.A(Nlumps-i+1);
    LTMB.A(i) = LTMB.A(Nlumps-i+1);
end

% Convert the area from mm^2 to m^2
LTM.A = [LTMB.A LTMA.A]*1e-6;
 
% Multiply the area with the convection coefficient to get the conductance
% to the ambient air
LTM.GA = c.therm.Conv.airNatural.nom.*LTM.A;

%% Water onvection

% Area = Top area + Side Area - Overlap
LTMA.WA(1)  = 0;
LTMB.WA(1)  = 0;
LTMA.WA(2)  = 0;
LTMB.WA(2)  = 0;
LTMA.WA(3)  = 0;
LTMB.WA(3)  = 0;
LTMA.WA(4)  = 0;
LTMB.WA(4)  = 0;
LTMA.WA(5)  = 0;
LTMB.WA(5)  = 0;
LTMA.WA(6)  = 43.660;
LTMB.WA(6)  = 314.33;
LTMA.WA(7)  = 0;
LTMB.WA(7)  = 0;
LTMA.WA(8)  = 0;
LTMB.WA(8)  = 0;
LTMA.WA(9)  = 0;
LTMB.WA(9)  = 0;
LTMA.WA(10) = 0;
LTMB.WA(10) = 0;
LTMA.WA(11) = 0;
LTMB.WA(11) = 0;
LTMA.WA(12) = 0;
LTMB.WA(12) = 0;
LTMA.WA(13) = 251.33;
LTMB.WA(13) = 251.33;
LTMA.WA(14) = 0;
LTMB.WA(14) = 0;
LTMA.WA(15) = 0;
LTMB.WA(15) = 0;
LTMA.WA(16) = 395.97;
LTMB.WA(16) = 395.97;
LTMA.WA(17) = 395.97;
LTMB.WA(17) = 395.97;
LTMA.WA(18) = 395.97;
LTMB.WA(18) = 395.97;
LTMA.WA(19) = 395.97;
LTMB.WA(19) = 395.97;
LTMA.WA(20) = 534.07;
LTMB.WA(20) = 534.07;
LTMA.WA(21) = 0;
LTMB.WA(21) = 0;
LTMA.WA(22) = 0;
LTMB.WA(22) = 0;
LTMA.WA(23) = 659.73;
LTMB.WA(23) = 659.73;
LTMA.WA(24) = 659.73;
LTMB.WA(24) = 659.73;
LTMA.WA(25) = 659.73;
LTMB.WA(25) = 659.73;
LTMA.WA(26) = 659.73;
LTMB.WA(26) = 659.73;
LTMA.WA(27) = 659.73;
LTMB.WA(27) = 659.73;
LTMA.WA(28) = 0;
LTMB.WA(28) = 0;
LTMA.WA(29) = 0;
LTMB.WA(29) = 0;
LTMA.WA(30) = 659.73;
LTMB.WA(30) = 659.73;
LTMA.WA(31) = 659.73;
LTMB.WA(31) = 659.73;
LTMA.WA(32) = 659.73;
LTMB.WA(32) = 659.73;
% Symmetry
for i = (ceil(Nlumps/2)+1):Nlumps
    % disp([num2str(i),'==',num2str(Nlumps-i+1)]) % For debugging
    LTMA.WA(i) = LTMA.WA(Nlumps-i+1);
    LTMB.WA(i) = LTMB.WA(Nlumps-i+1);
end

% Convert the area from mm^2 to m^2
LTM.WA = [LTMB.WA LTMA.WA]*1e-6;
 
% Multiply the area with the convection coefficient to get the convection
% to the flowing water
LTM.GWA = c.therm.Conv.waterForced.nom.*LTM.WA;
LTM.GWA_cond = c.therm.heatCond.water.*LTM.WA/0.005;




