%% Lokalizace zdroj� s v�ce simulovan�mi povrchov�mi elektrodami
clc;clear;close all;
ft_defaults;
%% Na�ten� pot�ebn�ch soubor�
load('mesh4down'); % mesh
load('headmodel4down'); % headmodel
load('gridinside'); % Grid uvnit� mozku
load('SRC/CBWJ13_P80_indexed_volume/lut.mat'); % lut pro atlas
load('atlas4down'); % atlas podvzorkovan� 4x
load('sourcemodel_atlas'); % sourcemodel + atlas = indexovan� zdroje
load('mri_segmented'); % mri mozku potkana segmentovan�
load('mri'); % nesegmentovan� mri mozku potkana
load('mri4down'); % podvzorkovan� nesegmentovan� mri mozku potkana
load('elec'); % NUDZ elektrody
load('leadfield'); % leadfield pro v�echny zdroje (NUDZ elektrody)
load('sourcemodel'); % sourcemodel

%% Generov�n� cosinusovky se SNR
oblast = 11; % index oblasti (11...hippocampus)
cela   = 1; % 1 pro celou oblast/ 0 pro 1 zdroj
moment = [1 0 0]; % V jak�ch os�ch jde sign�l [x y z];

fvz    = 1000; % vzorkovac� frekvence
f      = 35; % frekvence sign�lu
T      = 60; % d�lka sign�lu
t      = 0:1/fvz:T-1/fvz; % �asov� vektor
A      = 1; % amplituda
SNR    = 200; % Pom�r sign�l �um v dB
sinus  = A * cos(2*pi*f*t); % cosinusovka
signal = awgn(sinus,SNR); % Add white Gaussian noise to a signal.


indx   = find(sourcemodel_atlas.brick0==oblast); % V�b�r v�ech zdroj� podle oblasti

switch cela % V�b�r jestli jeden zdroj (0) nebo cel� oblast zdroj� (1)
    case 0
        zdroj_pozice     = randi(length(indx)); % V�b�r konkr�tn�ho zdroje z oblasti n�hodn�
        zdroj            = indx(zdroj_pozice); % Index zdroje na z�klad� n�hodn�ho v�b�ru
        zdroj_leadfield  = leadfield.leadfield(zdroj); % Pro NUDZ elektrody
        zdroj_leadfield  = zdroj_leadfield{1};
        signal_leadfield = cell(1,length(signal)); % Definov�n� prom�nn� typu cell
        potencial        = zeros(12,length(signal)); % Alokace pam�ti pro �asovou �adu ne * nt

        for i=1:length(signal)
            signal_leadfield{i} = signal(i)*zdroj_leadfield; % Pron�soben� vzork� sign�lu s leadfield
            potencial(:,i)      = moment*(signal_leadfield{i})'; % Pron�soben� maskou pro ur�en� moment�
        end
        
        % Um�st�n� simulovan�ho zdroje a elektrod
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
                potencial        = zeros(12,length(signal)); % Alokace pam�ti pro �asovou �adu ne * nt
                
                for i = 1:length(zdroj_leadfield) % Projde leadfield v�ech zdroj� vybran� oblasti
                   disp(i);
                    for j = 1:length(signal) % V�echny vzorky sign�lu
                       signal_leadfield = signal(j)*zdroj_leadfield{i}; % Pron�soben� vzork� sign�lu s aktu�ln�m leadfield
                       potencial(:,j)        = potencial(:,j) + (moment*(signal_leadfield)')'; % Pron�soben� vektorem moment� a p�i�ten� potenci�lu k p�edchoz�m
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
%% Fieldtrip form�t
time = (1:length(potencial))/fvz; % Celkov� �asov� vektor

data         = [];
data.label   = elec.label; % ozna�en� elektrod
data.fsample = 1000; % vzorkovac� frekvence
data.trial{1}= potencial; % data = elektrody * �asov� vzorky
data.time{1} = time; % �asov� vektor = 1* �asov� vzorky

cfg         = [];
cfg.length  = 2; % rozd�len� dat na x sekundov� �seky
data_trials = ft_redefinetrial(cfg, data);

cfg          = [];
data_preproc = ft_preprocessing(cfg,data_trials);

%% Frekven�n� anal�za
cfg           = [];
cfg.method    = 'mtmfft'; % pou�it� fft
cfg.output    = 'powandcsd'; % v�stupem power a cross spektra
cfg.tapsmofrq = 4;
cfg.foi       = 35; % kter� frekvence n�s zaj�m�
cfg.pad       = 'nextpow2'; % prodlou�en� po�tu vzork� na mocninu 2 pro fft
freq          = ft_freqanalysis(cfg, data_preproc);


%% Anal�za zdroj�
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
cfg.parameter = 'pow'; % power spektrum jako interpolovan� hodnoty
sourceInt = ft_sourceinterpolate(cfg,sourcePost_nocon,mridown);


%% Maska - pro zobrazen� oblasti pou�it�m mask jako funparameter na anatomy mridown
cfg = [];
cfg.inputcoord = 'unknown';
cfg.atlas = atlas4down;
cfg.roi = atlas4down.brick0label(:);
mask = ft_volumelookup(cfg, sourceInt); %Vytvo�en� masky pro zobrazen� zdroj� pouze v oblasti mozku
sourceInt.mask = mask;

%% Vymezen� zdroj� pouze na mozek - chyba interpolace
indexy = find(sourceInt.mask(:,:,:) == false);
sourceInt.pow(indexy)=0;

%% Zobrazen� do MRI
cfg = [];
%cfg.location = [x y z]
cfg.method = 'slice';
cfg.funparameter = 'pow'; % parcel-rozparcelovan� zdroje/pow-origin�l interpolovan�
ft_sourceplot(cfg,sourceInt,mridown)
hold on
text(100,200,'lambda 0.000000001%','Color','White')
hold off


%% Zobrazen� ve zdroj�ch
cfg = [];
cfg.method = 'vertex';
cfg.funparameter = 'pow'; % parcel-rozparcelovan� zdroje/pow-origin�l interpolovan�
hold on
ft_sourceplot(cfg,sourcePost_nocon)
ft_plot_mesh(mesh,'surfaceonly','yes','facealpha',0.1,'edgealpha',0.1) ;
if zdroj
line(linspace(-13,13,200),ones(1,200)*sourcemodel.pos(zdroj,2),ones(1,200)*sourcemodel.pos(zdroj,3),'LineWidth',2,'Color','r');
line(ones(1,200)*sourcemodel.pos(zdroj,1),linspace(-20,15,200),ones(1,200)*sourcemodel.pos(zdroj,3),'LineWidth',2,'Color','g');
line(ones(1,200)*sourcemodel.pos(zdroj,1),ones(1,200)*sourcemodel.pos(zdroj,2),linspace(-10,10,200),'LineWidth',2,'Color','b');
end
hold off

% %% Zdroje s velk�mi chybov�mi hodnotami
% zdroj_err = sourcePost_nocon.avg.pow;
% max_err = max(zdroj_err);
% err_ind = find(zdroj_err >=0.00001* max_err);
% sourcePost_nocon.avg.pow(err_ind) = 0;
