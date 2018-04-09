clc;clear variables;close all;
ft_defaults;
load('mesh4down');
% load('surfleadfield');
load('headmodel4down');
elec = [];
nsurfpoint = 5870; % Poèet povrchových bodù
skip = 20; % Kolikátý bod zahrnout do výpoètu
k = 1;
elec.elecpos = [];
for i=1:skip:nsurfpoint
    elec.label{k} = sprintf('%04d', i);
    elec.elecpos = [elec.elecpos;mesh.pos(i,:)];
    k=k+1;
end
elec.label = elec.label';
elec.pnt = elec.elecpos;
elec.unit = 'mm';

cfg                 = [];
cfg.headmodel       = headmodel;
cfg.sourceunits     = headmodel.unit;
cfg.elec            = elec;
cfg.grid.inside     = [1]; % Dipól je uvnitø
cfg.grid.pos        = [1 -4 2]; % Pozice dipólu
sourcemodel1         = ft_prepare_sourcemodel(cfg);

cfg       = [];
cfg.vol   = headmodel;  
cfg.elec  = elec;
cfg.grid  = sourcemodel1;
leadfield1 = ft_prepare_leadfield(cfg);

potencial = leadfield1.leadfield{1};
potencial = sqrt(potencial(:,1).^2 + potencial(:,2).^2 + potencial(:,3).^2);
% dipMom = [0,0,-1];
% potencial = dipMom * leadfield1.leadfield{1}';

% [headmodel, elec] = ft_prepare_vol_sens(headmodel, elec); % Úprava pozic elektrod pro pøiléhání

% figure('Name','Mesh + elektrody')
% hold on
% ft_plot_mesh(mesh,'surfaceonly','yes','facealpha',0.2,'edgealpha',0.2);
% ft_plot_sens(elec1,'style','r.','MarkerSize',1);
% hold off
figure('Name','Mesh + eldy + topo + osy potencial')
hold on
ft_plot_topo3d(leadfield1.cfg.elec.chanpos,potencial,'facealpha',0.6);
% ft_plot_sens(elec1,'style','r.');
% ft_plot_mesh(mesh,'surfaceonly','yes','facealpha',0.2,'edgealpha',0.2) ;
line(linspace(-13,13,200),ones(1,200)*sourcemodel1.pos(2),ones(1,200)*sourcemodel1.pos(3),'LineWidth',2,'Color','r')
line(ones(1,200)*sourcemodel1.pos(1),linspace(-20,15,200),ones(1,200)*sourcemodel1.pos(3),'LineWidth',2,'Color','g')
line(ones(1,200)*sourcemodel1.pos(1),ones(1,200)*sourcemodel1.pos(2),linspace(-10,10,200),'LineWidth',2,'Color','b')

hold off