clc;
clear all;
close all;
%% Nastavení cesty k fieldtrip toolboxu
restoredefaultpath
addpath C:\fieldtrip-20180309 %Notebook
ft_defaults

%% Naètení modelu mozku
mesh = load('src/HexahedralBrain_300k.mat');
mesh = mesh.mesh;
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
headmodel.unit   = 'mm';
% headmodel = load('src/headmodel.mat');
% headmodel = headmodel.vol;

%% Elektrody
elec = [];
elec.pnt = [-3,5,6;3,5,6;
            -6,0,7;0,0,7;6,0,7;
            -6,-6,7;0,-6,7;6,-6,7;
            -6,-13,7;0,-13,7;6,-13,7;]; %Vymyšlené pozice 11 elektrod
for i = 1:11
elec.label{i} = sprintf('%02d', i); % Oznaèení elektrod
end
elec.elecpos = elec.pnt;
elec.unit = 'mm';

[headmodel, elec] = ft_prepare_vol_sens(headmodel, elec); % Úprava pozic elektrod pro pøiléhání

figure('Name','Mesh + elektrody')
hold on
ft_plot_mesh(mesh,'surfaceonly','yes','facecolor','brain') 
ft_plot_sens(elec,'style','r.');
hold off

%% Source model = headmodel + pozice elektrod + pozice dopólù
cfg                 = [];
cfg.headmodel       = headmodel;
cfg.sourceunits     = headmodel.unit;
cfg.elec            = elec;
cfg.grid.pos        = [1 -4 6]; % Pozice dipólu
sourcemodel         = ft_prepare_sourcemodel(cfg);

%% Prepare leadfield = Výpoèet potenciálù ve 3 osách pro každou elektrodu
cfg       = [];
cfg.vol   = headmodel;  
cfg.elec  = elec;
cfg.grid  = sourcemodel;
leadfield = ft_prepare_leadfield(cfg);

%% Výpoèet potenciálù na elektrodách
potencial = leadfield.leadfield{1};
potencial = sqrt(potencial(:,1).^2 + potencial(:,2).^2 + potencial(:,3).^2);

%% Rozložení potenciálu po povrchu mozku
leadfield.cfg.elec.chanpos(:,3) = leadfield.cfg.elec.chanpos(:,3) + 0.8;
figure('Name','Mesh + eldy + topo + osy potencial')
hold on
ft_plot_topo3d(leadfield.cfg.elec.chanpos,potencial);
ft_plot_sens(elec,'style','r.');
ft_plot_mesh(mesh,'surfaceonly','yes','facecolor','brain') ;
line(linspace(-13,13,200),ones(1,200)*sourcemodel.pos(2),ones(1,200)*sourcemodel.pos(3),'LineWidth',2,'Color','r')
line(ones(1,200)*sourcemodel.pos(1),linspace(-20,15,200),ones(1,200)*sourcemodel.pos(3),'LineWidth',2,'Color','g')
line(ones(1,200)*sourcemodel.pos(1),ones(1,200)*sourcemodel.pos(2),linspace(-10,10,200),'LineWidth',2,'Color','b')



