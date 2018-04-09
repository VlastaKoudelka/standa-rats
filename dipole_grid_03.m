%% Møížka v mozku = dipóly a spoèítat pro nì leadfield pro 12 elektrod
clc;clear variables;close all;
ft_defaults;

%% Naètení meshe a headmodelu
load('mesh4down');
load('headmodel4down');

%% Elektrody - NUDZ uspoøádání
electrodeTable = readtable('SRC/ratElectrode12.txt'); %soubor s reálnými elektrodami
elPositions = [electrodeTable.X,electrodeTable.Y,electrodeTable.Z]; % Pozice  elektrod

elec = [];
elec.label = electrodeTable.Name;
elec.pnt = [electrodeTable.X,electrodeTable.Y,electrodeTable.Z];
elec.elecpos = elec.pnt;
elec.unit = 'mm';

%% Grid
step    = 1;
grid_x  = (min(mesh.pos(:,1)): step : max(mesh.pos(:,1)))';
grid_y  = (min(mesh.pos(:,2)): step : max(mesh.pos(:,2)))';
grid_z  = (min(mesh.pos(:,3)): step : max(mesh.pos(:,3)))';
[X,Y,Z] = meshgrid(grid_x,grid_y,grid_z);
npoints = size(X,1)*size(X,2)*size(X,3);
grid    = [reshape(X,[npoints,1]),reshape(Y,[npoints,1]),reshape(Z,[npoints,1])];
[inside] = ft_inside_vol(grid, headmodel);
gridinside = grid(inside==1,:);
scatter3(gridinside(:,1),gridinside(:,2),gridinside(:,3));
%% Source model = headmodel + pozice elektrod + pozice dopólù
cfg                 = [];
cfg.headmodel       = headmodel;
cfg.sourceunits     = headmodel.unit;
cfg.elec            = elec;
cfg.grid.pos        = gridinside;
sourcemodel         = ft_prepare_sourcemodel(cfg);
  
%% Prepare leadfield = Výpoèet potenciálù ve 3 osách pro každou elektrodu
cfg       = [];
cfg.vol   = headmodel;  
cfg.elec  = elec;
cfg.grid  = sourcemodel;
leadfield = ft_prepare_leadfield(cfg);

%% Výbìr dipólù pro zobrazení
dipole = 780; % Poøadové èíslo dipólu

potencial = leadfield.leadfield{dipole};
potencial = sqrt(potencial(:,1).^2 + potencial(:,2).^2 + potencial(:,3).^2);

leadfield.cfg.elec.chanpos(:,3) = leadfield.cfg.elec.chanpos(:,3) + 0.7;
figure('Name','Mesh + eldy + topo + osy potencial')
hold on
ft_plot_topo3d(leadfield.cfg.elec.chanpos,potencial,'facealpha',0.6,'refine',5);
ft_plot_sens(elec,'style','r.');
ft_plot_mesh(mesh,'surfaceonly','yes','facealpha',0.2) ;
line(linspace(-13,13,200),ones(1,200)*sourcemodel.pos(dipole,2),ones(1,200)*sourcemodel.pos(dipole,3),'LineWidth',2,'Color','r')
line(ones(1,200)*sourcemodel.pos(dipole,1),linspace(-20,15,200),ones(1,200)*sourcemodel.pos(dipole,3),'LineWidth',2,'Color','g')
line(ones(1,200)*sourcemodel.pos(dipole,1),ones(1,200)*sourcemodel.pos(dipole,2),linspace(-10,10,200),'LineWidth',2,'Color','b')

