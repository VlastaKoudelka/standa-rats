clc;clear;close all;
ft_defaults;
%% Na�ten� pot�ebn�ch soubor�
load('mesh4down'); % mesh
load('headmodel4down'); % headmodel
load('elec'); % NUDZ elektrody
load('gridinside'); % Grid uvnit� mozku
load('SRC/CBWJ13_P80_indexed_volume/lut.mat'); % lut pro atlas
load('atlas4down'); % atlas podvzorkovan� 4x
load('sourcemodel'); % sourcemodel
load('sourcemodel_atlas'); % sourcemodel + atlas = indexovan� zdroje
load('leadfield'); % leadfield pro v�echny zdroje (NUDZ elektrody)

%% Na�ten� atlasu
% file = 'SRC/CBWJ13_P80_indexed_volume/CBWJ13_P80_indexed_volume.nii';
% atlas = ft_read_atlas(file);
% atlas.brick0label  = cellstr(labels.Name');
%% Podvzorkov�n� atlasu
% cfg            = [];
% cfg.downsample = 4; % Podvzorkov�n� 4x ... 16 voxely -> 1 voxel
% [atlas4down]      = ft_volumedownsample(cfg,atlas); % Pro pou�itelnou velikost 
atlas4down.unit='mm';
%% Zobrazen� atlasu
% cfg = [];
% cfg.atlas = atlas4down;
% cfg.funparameter = 'brick0';
% cfg.funcolormap = 'lines';
% ft_sourceplot(cfg, atlas4down);
% 
% imagesc(atlas4down.brick0(:,:,20)); % 1 �ez kolm� na osu z
% 
% figure
% for i = 1:size(atlas4down.brick0,3)
%     hold on
%     imagesc(atlas4down.brick0(:,:,i));
%     pause(0.1);
% end


%% P�i�azen� atlasu k sourcemodelu = indexov�n� zdroj�
% cfg = []; 
% cfg.interpmethod = 'nearest'; 
% cfg.parameter = 'brick0'; 
% sourcemodel_atlas = ft_sourceinterpolate(cfg, atlas4down, sourcemodel); 

%% Headmodel + zdroje vybran� oblasti
oblast = 11;
indx = find(sourcemodel_atlas.brick0==oblast);
figure('Name',labels.Name{oblast});
hold on;
scatter3(gridinside(indx,1),gridinside(indx,2),gridinside(indx,3),...
        'MarkerEdgeColor','k',...
        'MarkerFaceColor',[labels.Red(oblast)/255 labels.Green(oblast)/255 labels.Blue(oblast)/255]);
ft_plot_mesh(mesh,'surfaceonly','yes','facealpha',0.2,'edgealpha',0.2);
hold off;
%% Grid barevn� odli�en� oblasti
figure('Name','Barevn� odli�en� zdroje podle oblast�');
hold on;
for i=1:size(labels,1)
    index = find(sourcemodel_atlas.brick0==i);
    scatter3(gridinside(index,1),gridinside(index,2),gridinside(index,3),...
    'MarkerEdgeColor','None','MarkerFaceColor',...
    [labels.Red(i)/255 labels.Green(i)/255 labels.Blue(i)/255]);
end
zlim([-25 15]);
xlim([-9 9]);
ylim([-25 15]);
hold off;

%% Generov�n� cosinusovky ve vybran� oblasti
oblast = 11; % index oblasti (11...hippocampus)

fvz    = 1000; % vzorkovac� frekvence
f      = 35; % frekvence sign�lu
T      = 60; % d�lka sign�lu
t      = 0:1/fvz:T-1/fvz; % �asov� vektor
A      = 1; % amplituda
signal = A * cos(2*pi*f*t)+(rand(1,fvz*T)-.5); % sinusovka

indx             = find(sourcemodel_atlas.brick0==oblast); % V�b�r v�ech zdroj� podle oblasti
zdroj_pozice     = randi(length(indx)); % V�b�r konkr�tn�ho zdroje z oblasti n�hodn�
zdroj            = indx(zdroj_pozice);
zdroj_leadfield  = leadfield.leadfield(zdroj); % Pro NUDZ elektrody
zdroj_leadfield  = zdroj_leadfield{1};
signal_leadfield = cell(1,length(signal));

maska = [1 0 0]; % V jak�ch os�ch jde sign�l [x y z];

for i=1:length(signal)
    signal_leadfield{i} = signal(i)*zdroj_leadfield;
end

figure('Name',labels.Name{oblast});
hold on
subplot(2,2,1)
scatter3(gridinside(zdroj,1),gridinside(zdroj,2),gridinside(zdroj,3),...
'MarkerEdgeColor','k',...
'MarkerFaceColor',[labels.Red(oblast)/255 labels.Green(oblast)/255 labels.Blue(oblast)/255]);
line(linspace(-6,6,200),zeros(1,200),zeros(1,200),'LineWidth',2,'Color','r');
line(zeros(1,200),linspace(-10,6,200),zeros(1,200),'LineWidth',2,'Color','g');
text(6,0,0,'+X');
text(-8,0,0,'-X');
text(0.5,6,0,'+Y');
text(0,-11,0,'-Y');


for i=1:length(signal)
    potencial = maska(1)*signal_leadfield{i}(:,1) +...
                maska(2)*signal_leadfield{i}(:,2) + maska(3)*signal_leadfield{i}(:,3);      
    subplot(2,2,1)
    title('potenci�l mezi elektrodami + osy + zdroj');
    ft_plot_topo3d(leadfield.cfg.elec.chanpos,potencial,'facealpha',0.6,'refine',2); % NUDZ
    view(0,90);
    subplot(2,2,2)
    hold on
    title('Hloubka zdroje');
    ft_plot_topo3d(leadfield.cfg.elec.chanpos,potencial,'facealpha',0.6,'refine',2); % NUDZ
    text(0,6,8,'+Y');
    text(0,-11.5,8,'-Y');
    
    if i<90
    view(i,90-i);
    end
    
    line(zeros(1,200),linspace(-10,6,200),8*ones(1,200),'LineWidth',2,'Color','g');
    scatter3(gridinside(zdroj,1),gridinside(zdroj,2),gridinside(zdroj,3),...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor',[labels.Red(oblast)/255 labels.Green(oblast)/255 labels.Blue(oblast)/255]);

    subplot(2,2,3:4);
    if i>100
    plot(i-50:i+50,signal(i-50:i+50));
    title('Pr�b�h zdroje');
    line([i i],[-max(signal) max(signal)],'Color','r');
    xlim([i-50 i+50]);
    end
    pause(0.1);
end






