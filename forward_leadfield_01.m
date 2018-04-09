clc;
clear;
close all;
%% Nastaven� cesty k fieldtrip toolboxu
% restoredefaultpath
% addpath C:\fieldtrip-20180309 %Notebook
ft_defaults
tic
%% Na�ten� MRI mozku potkana
filePath = 'SRC/NIfTILow/'; % Cesta do adres��e souboru
file     = 'volume_Brain.nii'; % N�zev souboru
mri      = ft_read_mri([filePath file]); % Na�ten� dat MRI

%% Manu�ln� MRI segmentace mozku 
brain                = zeros(size(mri.anatomy)); % Alokace pam�ti pro vysegmentovan� mozek
brain(mri.anatomy>0) = 1; % Kdy� hodnota voxelu > 0, tak se tam nach�z� mozek (logick� 1)
brain                = logical(brain); % P�etypov�n� na logickou prom�nnou

segMri       = mri; % Vytvo�en� struktury segMri z mri
segMri       = rmfield(segMri,'anatomy'); % Odstran�n� pole anatomy
segMri.unit  = 'mm'; % Ur�en� jednotek
segMri.brain = brain; % P�e�azen� matice vysegmentovan�ho mozku

%% Podvzorkov�n� MRI
% cfg            = [];
% cfg.downsample = 2; % Podvzorkov�n� 2x ... 4 voxely -> 1 voxel
% [mridown]      = ft_volumedownsample(cfg,segMri); % Pro pou�itelnou velikost 

%% Prepare mesh - vytvo�en� meshe z MRI
cfg              = [];
cfg.method       = 'tetrahedral';
cfg.downsample   = 4;
mesh             = ft_prepare_mesh(cfg,segMri);
mesh.labels      = ones(size(mesh.tet,1),1);
mesh.tissue      = ones(size(mesh.tet,1),1);
mesh.tissuelabel = {'brain'};
mesh.tet         = [mesh.tet(:,1:2) mesh.tet(:,4) mesh.tet(:,3)];
figure('Name','Mesh + sou�adn� osy')
hold on
ft_plot_mesh(mesh,'surfaceonly','yes','facecolor','brain') 
line(linspace(-13,13,200),zeros(1,200),zeros(1,200),'LineWidth',2,'Color','r')
line(zeros(1,200),linspace(-20,15,200),zeros(1,200),'LineWidth',2,'Color','g')
line(zeros(1,200),zeros(1,200),linspace(-10,10,200),'LineWidth',2,'Color','b')
hold off

%% Vodivostn� model = head model
cfg              = [];
cfg.method       = 'simbio'; % Metoda v�po�tu
cfg.conductivity = [0.33] ; % Vodivost v S/m
headmodel        = ft_prepare_headmodel(cfg,mesh);
headmodel.unit   = 'mm'; % P�epis jednotek na mm

%% Elektrody
electrodeTable = readtable('SRC/ratElectrode12.txt'); %soubor s re�ln�mi elektrodami
elPositions = [electrodeTable.X,electrodeTable.Y,electrodeTable.Z]; % Pozice  elektrod

elec = [];
elec.label = electrodeTable.Name;
elec.pnt = [electrodeTable.X,electrodeTable.Y,electrodeTable.Z];
elec.elecpos = elec.pnt;
elec.unit = 'mm';

% [headmodel, elec] = ft_prepare_vol_sens(headmodel, elec); % �prava pozic elektrod pro p�il�h�n�

figure('Name','Mesh + elektrody')
hold on
ft_plot_mesh(mesh,'surfaceonly','yes','facecolor','brain');
ft_plot_sens(elec,'style','r.');
hold off

%% Source model = headmodel + pozice elektrod + pozice dop�l�
cfg                 = [];
cfg.headmodel       = headmodel;
cfg.sourceunits     = headmodel.unit;
cfg.elec            = elec;
cfg.grid.inside     = [1]; % Dip�l je uvnit�
cfg.grid.pos        = [1 -4 5]; % Pozice dip�lu
sourcemodel         = ft_prepare_sourcemodel(cfg);
  
%% Prepare leadfield = V�po�et potenci�l� ve 3 os�ch pro ka�dou elektrodu
cfg       = [];
cfg.vol   = headmodel;  
cfg.elec  = elec;
cfg.grid  = sourcemodel;
leadfield = ft_prepare_leadfield(cfg);

%% V�po�et potenci�l� na elektrod�ch - Geometrick� sou�et potenci�l� v os�ch x,y,z
potencial = leadfield.leadfield{1};
potencial = sqrt(potencial(:,1).^2 + potencial(:,2).^2 + potencial(:,3).^2);

%% Rozlo�en� potenci�lu po povrchu mozku
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

