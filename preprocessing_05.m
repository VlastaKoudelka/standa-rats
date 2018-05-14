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
load('mri_segmented'); % mri mozku potkana segmentovan�
load('mri'); % nesegmentovan� mri mozku potkana
load('mri4down'); % podvzorkovan� nesegmentovan� mri mozku potkana

%% Generov�n� cosinusovky se SNR
oblast = 11; % index oblasti (11...hippocampus)
cela   = 1; % 1 pro celou oblast/ 0 pro 1 zdroj
moment = [0 0 1]; % V jak�ch os�ch jde sign�l [x y z];

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

end
%% Fieldtrip format
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
cfg.method            = 'mne';
cfg.frequency         = 35;  
cfg.grid              = sourcemodel; 
cfg.headmodel         = headmodel;
cfg.mne.projectnoise  = 'yes';
cfg.mne.lambda        = 0;
cfg.elec              = elec;

sourcePost_nocon      = ft_sourceanalysis(cfg, freq);

%% Interpolace ze sourcemodel na MRI (volume)

cfg = [];
cfg.interpmethod = 'nearest'; % metoda interpolace
cfg.parameter = 'pow'; % power spektrum jako interpolovan� hodnoty
sourceInt = ft_sourceinterpolate(cfg,sourcePost_nocon,atlas4down);
sourceInt.unit = 'mm';
%% Parcelace zdroj� podle oblast� atlasu

cfg = [];
cfg.method = 'mean'; % st�edn� hodnota
cfg.parcellation = 'brick0';
parcel = ft_sourceparcellate(cfg, sourceInt, atlas4down); % Parcelace zdroj�
sourceInt.parcel = zeros(length(sourceInt.pow),1); % Vytvo�en� nov� prom�nn� pro ulo�en� parcelovan�ch zdroj�
for i=1:length(parcel.pow)
      sourceInt.parcel(find(atlas4down.brick0==i))=parcel.pow(i); % P�i�azen� st�edn� hodnoty oblasti do v�ech voxel� oblasti
end

%% Skutecne zdroje

cfg = [];
cfg.method = 'mean'; % st�edn� hodnota
cfg.parcellation = 'brick0';
parcel = ft_sourceparcellate(cfg, sourceInt, atlas4down); % Parcelace zdroj�
sourceInt.trueSource = zeros(length(sourceInt.pow),1); % Vytvo�en� nov� prom�nn� pro ulo�en� umelych zdroju

sourceInt.trueSource(find(atlas4down.brick0==oblast))=1; % P�i�azen� st�edn� hodnoty oblasti do v�ech voxel� oblasti

%% Maska

cfg = [];
cfg.inputcoord = 'unknown';
cfg.atlas = atlas4down;
cfg.roi = atlas4down.brick0label(:);
mask = ft_volumelookup(cfg, sourceInt); %Vytvo�en� masky pro zobrazen� zdroj� pouze v oblasti mozku

sourceInt.mask = mask;


%% Zobrazen�
cfg = [];
cfg.method = 'ortho';
cfg.funparameter = 'parcel'; % parcel-rozparcelovan� zdroje/pow-origin�l interpolovan�
cfg.maskparameter = 'mask'; % maska
cfg.atlas = atlas4down;
%cfg.mesh = mesh;
ft_sourceplot(cfg,sourceInt)

%% Zobrazen�
cfg = [];
cfg.method = 'ortho';
cfg.funparameter = 'trueSource'; % parcel-rozparcelovan� zdroje/pow-origin�l interpolovan�
cfg.maskparameter = 'mask'; % maska
cfg.atlas = atlas4down;
%cfg.mesh = mesh;
ft_sourceplot(cfg,sourceInt)




