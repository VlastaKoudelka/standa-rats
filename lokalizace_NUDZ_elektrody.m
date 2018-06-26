%% Lokalizace zdrojù s více simulovanými povrchovými elektrodami
clc;clear;close all;
ft_defaults;
%% Naètení potøebných souborù
load('mesh4down'); % mesh
load('headmodel4down'); % headmodel
load('gridinside'); % Grid uvnitø mozku
load('SRC/CBWJ13_P80_indexed_volume/lut.mat'); % lut pro atlas
load('atlas4down'); % atlas podvzorkovaný 4x
load('sourcemodel_atlas'); % sourcemodel + atlas = indexované zdroje
load('mri_segmented'); % mri mozku potkana segmentované
load('mri'); % nesegmentované mri mozku potkana
load('mri4down'); % podvzorkované nesegmentované mri mozku potkana
load('elec'); % NUDZ elektrody
load('leadfield'); % leadfield pro všechny zdroje (NUDZ elektrody)
load('sourcemodel'); % sourcemodel

%% Generování cosinusovky se SNR
oblast = 11; % index oblasti (11...hippocampus)
cela   = 1; % 1 pro celou oblast/ 0 pro 1 zdroj
moment = [1 0 0]; % V jakých osách jde signál [x y z];

fvz    = 1000; % vzorkovací frekvence
f      = 35; % frekvence signálu
T      = 60; % délka signálu
t      = 0:1/fvz:T-1/fvz; % èasový vektor
A      = 1; % amplituda
SNR    = 200; % Pomìr signál šum v dB
sinus  = A * cos(2*pi*f*t); % cosinusovka
signal = awgn(sinus,SNR); % Add white Gaussian noise to a signal.


indx   = find(sourcemodel_atlas.brick0==oblast); % Výbìr všech zdrojù podle oblasti

switch cela % Výbìr jestli jeden zdroj (0) nebo celá oblast zdrojù (1)
    case 0
        zdroj_pozice     = randi(length(indx)); % Výbìr konkrétního zdroje z oblasti náhodnì
        zdroj            = indx(zdroj_pozice); % Index zdroje na základì náhodného výbìru
        zdroj_leadfield  = leadfield.leadfield(zdroj); % Pro NUDZ elektrody
        zdroj_leadfield  = zdroj_leadfield{1};
        signal_leadfield = cell(1,length(signal)); % Definování promìnné typu cell
        potencial        = zeros(12,length(signal)); % Alokace pamìti pro èasovou øadu ne * nt

        for i=1:length(signal)
            signal_leadfield{i} = signal(i)*zdroj_leadfield; % Pronásobení vzorkù signálu s leadfield
            potencial(:,i)      = moment*(signal_leadfield{i})'; % Pronásobení maskou pro urèení momentù
        end
        
        % Umístìní simulovaného zdroje a elektrod
                figure('Name','Mesh + eldy + osy potencial')
                hold on
                ft_plot_sens(elec,'style','r.');
                ft_plot_mesh(mesh,'surfaceonly','yes','facealpha',0.2,'edgealpha',0.2) ;
                line(linspace(-13,13,200),ones(1,200)*sourcemodel.pos(zdroj,2),ones(1,200)*sourcemodel.pos(zdroj,3),'LineWidth',2,'Color','r');
                line(ones(1,200)*sourcemodel.pos(zdroj,1),linspace(-20,15,200),ones(1,200)*sourcemodel.pos(zdroj,3),'LineWidth',2,'Color','g');
                line(ones(1,200)*sourcemodel.pos(zdroj,1),ones(1,200)*sourcemodel.pos(zdroj,2),linspace(-10,10,200),'LineWidth',2,'Color','b');
                for i=1:length(elec.label)
                text(elec.pnt(i,1),elec.pnt(i,2),elec.pnt(i,3),elec.label(i),'HorizontalAlignment','left','FontSize',8);
                end
                view(0,90);

                hold off
    case 1
                zdroj_leadfield  = leadfield.leadfield(indx); % Pro NUDZ elektrody
                potencial        = zeros(12,length(signal)); % Alokace pamìti pro èasovou øadu ne * nt
                
                for i = 1:length(zdroj_leadfield) % Projde leadfield všech zdrojù vybrané oblasti
                   disp(i);
                    for j = 1:length(signal) % Všechny vzorky signálu
                       signal_leadfield = signal(j)*zdroj_leadfield{i}; % Pronásobení vzorkù signálu s aktuálním leadfield
                       potencial(:,j)        = potencial(:,j) + (moment*(signal_leadfield)')'; % Pronásobení vektorem momentù a pøiètení potenciálu k pøedchozím
                   end
                end
                figure('Name',labels.Name{oblast});
                hold on;
                scatter3(gridinside(indx,1),gridinside(indx,2),gridinside(indx,3),...
                  'MarkerEdgeColor','k',...
                 'MarkerFaceColor',[labels.Red(oblast)/255 labels.Green(oblast)/255 labels.Blue(oblast)/255]);
                ft_plot_mesh(mesh,'surfaceonly','yes','facealpha',0.2,'edgealpha',0.2);
                                view(0,90);

                hold off;
                          
end
%% Fieldtrip formát
time = (1:length(potencial))/fvz; % Celkový èasový vektor

data         = [];
data.label   = elec.label; % oznaèení elektrod
data.fsample = 1000; % vzorkovací frekvence
data.trial{1}= potencial; % data = elektrody * èasové vzorky
data.time{1} = time; % èasový vektor = 1* èasové vzorky

cfg         = [];
cfg.length  = 2; % rozdìlení dat na x sekundové úseky
data_trials = ft_redefinetrial(cfg, data);

cfg          = [];
data_preproc = ft_preprocessing(cfg,data_trials);

%% Frekvenèní analýza
cfg           = [];
cfg.method    = 'mtmfft'; % použití fft
cfg.output    = 'powandcsd'; % výstupem power a cross spektra
cfg.tapsmofrq = 4;
cfg.foi       = 35; % která frekvence nás zajímá
cfg.pad       = 'nextpow2'; % prodloužení poètu vzorkù na mocninu 2 pro fft
freq          = ft_freqanalysis(cfg, data_preproc);


%% Analýza zdrojù
cfg                   = []; 
cfg.method            = 'dics';
cfg.frequency         = 35;  
cfg.grid              = sourcemodel; 
cfg.headmodel         = headmodel;
cfg.dics.projectnoise  = 'yes';
cfg.dics.lambda        = '0.000000001%';
cfg.elec              = elec;

sourcePost_nocon      = ft_sourceanalysis(cfg, freq);

%% Interpolace ze sourcemodel na MRI (volume)
cfg = [];
cfg.interpmethod = 'nearest'; % metoda interpolace
cfg.parameter = 'pow'; % power spektrum jako interpolované hodnoty
sourceInt = ft_sourceinterpolate(cfg,sourcePost_nocon,mridown);


%% Maska - pro zobrazení oblasti použitím mask jako funparameter na anatomy mridown
cfg = [];
cfg.inputcoord = 'unknown';
cfg.atlas = atlas4down;
cfg.roi = atlas4down.brick0label(:);
mask = ft_volumelookup(cfg, sourceInt); %Vytvoøení masky pro zobrazení zdrojù pouze v oblasti mozku
sourceInt.mask = mask;

%% Vymezení zdrojù pouze na mozek - chyba interpolace
indexy = find(sourceInt.mask(:,:,:) == false);
sourceInt.pow(indexy)=0;

%% Zobrazení do MRI
cfg = [];
%cfg.location = [x y z]
cfg.method = 'slice';
cfg.funparameter = 'pow'; % parcel-rozparcelované zdroje/pow-originál interpolované
ft_sourceplot(cfg,sourceInt,mridown)
hold on
text(100,200,'lambda 0.000000001%','Color','White')
hold off


%% Zobrazení ve zdrojích
cfg = [];
cfg.method = 'vertex';
cfg.funparameter = 'pow'; % parcel-rozparcelované zdroje/pow-originál interpolované
hold on
ft_sourceplot(cfg,sourcePost_nocon)
ft_plot_mesh(mesh,'surfaceonly','yes','facealpha',0.1,'edgealpha',0.1) ;
if zdroj
line(linspace(-13,13,200),ones(1,200)*sourcemodel.pos(zdroj,2),ones(1,200)*sourcemodel.pos(zdroj,3),'LineWidth',2,'Color','r');
line(ones(1,200)*sourcemodel.pos(zdroj,1),linspace(-20,15,200),ones(1,200)*sourcemodel.pos(zdroj,3),'LineWidth',2,'Color','g');
line(ones(1,200)*sourcemodel.pos(zdroj,1),ones(1,200)*sourcemodel.pos(zdroj,2),linspace(-10,10,200),'LineWidth',2,'Color','b');
end
hold off

% %% Zdroje s velkými chybovými hodnotami
% zdroj_err = sourcePost_nocon.avg.pow;
% max_err = max(zdroj_err);
% err_ind = find(zdroj_err >=0.00001* max_err);
% sourcePost_nocon.avg.pow(err_ind) = 0;
