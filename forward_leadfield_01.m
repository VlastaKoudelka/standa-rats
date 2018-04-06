clc;
clear;
close all;
%% Nastavení cesty k fieldtrip toolboxu
% restoredefaultpath
% addpath C:\fieldtrip-20180309 %Notebook
ft_defaults
tic
%% Naètení MRI mozku potkana
filePath = 'SRC/NIfTILow/'; % Cesta do adresáøe souboru
file     = 'volume_Brain.nii'; % Název souboru
mri      = ft_read_mri([filePath file]); % Naètení dat MRI

%% Manuální MRI segmentace mozku 
brain                = zeros(size(mri.anatomy)); % Alokace pamìti pro vysegmentovaný mozek
brain(mri.anatomy>0) = 1; % Když hodnota voxelu > 0, tak se tam nachází mozek (logická 1)
brain                = logical(brain); % Pøetypování na logickou promìnnou

segMri       = mri; % Vytvoøení struktury segMri z mri
segMri       = rmfield(segMri,'anatomy'); % Odstranìní pole anatomy
segMri.unit  = 'mm'; % Urèení jednotek
segMri.brain = brain; % Pøeøazení matice vysegmentovaného mozku

%% Podvzorkování MRI
% cfg            = [];
% cfg.downsample = 2; % Podvzorkování 2x ... 4 voxely -> 1 voxel
% [mridown]      = ft_volumedownsample(cfg,segMri); % Pro použitelnou velikost 

%% Prepare mesh - vytvoøení meshe z MRI
cfg              = [];
cfg.method       = 'tetrahedral';
cfg.downsample   = 4;
mesh             = ft_prepare_mesh(cfg,segMri);
mesh.labels      = ones(size(mesh.tet,1),1);
mesh.tissue      = ones(size(mesh.tet,1),1);
mesh.tissuelabel = {'brain'};
mesh.tet         = [mesh.tet(:,1:2) mesh.tet(:,4) mesh.tet(:,3)];
figure('Name','Mesh + souøadné osy')
hold on
ft_plot_mesh(mesh,'surfaceonly','yes','facecolor','brain') 
line(linspace(-13,13,200),zeros(1,200),zeros(1,200),'LineWidth',2,'Color','r')
line(zeros(1,200),linspace(-20,15,200),zeros(1,200),'LineWidth',2,'Color','g')
line(zeros(1,200),zeros(1,200),linspace(-10,10,200),'LineWidth',2,'Color','b')
hold off

%% Vodivostní model = head model
cfg              = [];
cfg.method       = 'simbio'; % Metoda výpoètu
cfg.conductivity = [0.33] ; % Vodivost v S/m
headmodel        = ft_prepare_headmodel(cfg,mesh);
headmodel.unit   = 'mm'; % Pøepis jednotek na mm

%% Elektrody
electrodeTable = readtable('SRC/ratElectrode12.txt'); %soubor s reálnými elektrodami
elPositions = [electrodeTable.X,electrodeTable.Y,electrodeTable.Z]; % Pozice  elektrod

elec = [];
elec.label = electrodeTable.Name;
elec.pnt = [electrodeTable.X,electrodeTable.Y,electrodeTable.Z];
elec.elecpos = elec.pnt;
elec.unit = 'mm';

% [headmodel, elec] = ft_prepare_vol_sens(headmodel, elec); % Úprava pozic elektrod pro pøiléhání

figure('Name','Mesh + elektrody')
hold on
ft_plot_mesh(mesh,'surfaceonly','yes','facecolor','brain');
ft_plot_sens(elec,'style','r.');
hold off

%% Source model = headmodel + pozice elektrod + pozice dopólù
cfg                 = [];
cfg.headmodel       = headmodel;
cfg.sourceunits     = headmodel.unit;
cfg.elec            = elec;
cfg.grid.inside     = [1]; % Dipól je uvnitø
cfg.grid.pos        = [1 -4 5]; % Pozice dipólu
sourcemodel         = ft_prepare_sourcemodel(cfg);
  
%% Prepare leadfield = Výpoèet potenciálù ve 3 osách pro každou elektrodu
cfg       = [];
cfg.vol   = headmodel;  
cfg.elec  = elec;
cfg.grid  = sourcemodel;
leadfield = ft_prepare_leadfield(cfg);

%% Výpoèet potenciálù na elektrodách - Geometrický souèet potenciálù v osách x,y,z
potencial = leadfield.leadfield{1};
potencial = sqrt(potencial(:,1).^2 + potencial(:,2).^2 + potencial(:,3).^2);

%% Rozložení potenciálu po povrchu mozku
leadfield.cfg.elec.chanpos(:,3) = leadfield.cfg.elec.chanpos(:,3) + 0.7;
figure('Name','Mesh + eldy + topo + osy potencial')
hold on
ft_plot_topo3d(leadfield.cfg.elec.chanpos,potencial,'facealpha',0.6,'refine',5);
ft_plot_sens(elec,'style','r.');
ft_plot_mesh(mesh,'surfaceonly','yes','facealpha',0.2,'edgealpha',0.2) ;
line(linspace(-13,13,200),ones(1,200)*sourcemodel.pos(2),ones(1,200)*sourcemodel.pos(3),'LineWidth',2,'Color','r')
line(ones(1,200)*sourcemodel.pos(1),linspace(-20,15,200),ones(1,200)*sourcemodel.pos(3),'LineWidth',2,'Color','g')
line(ones(1,200)*sourcemodel.pos(1),ones(1,200)*sourcemodel.pos(2),linspace(-10,10,200),'LineWidth',2,'Color','b')
toc

