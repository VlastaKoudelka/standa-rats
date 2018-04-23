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

atlas4down.unit='mm';

%% Generov�n� cosinusovky se SNR
oblast = 11; % index oblasti (11...hippocampus)

fvz    = 1000; % vzorkovac� frekvence
f      = 35; % frekvence sign�lu
T      = 60; % d�lka sign�lu
t      = 0:1/fvz:T-1/fvz; % �asov� vektor
A      = 1; % amplituda
SNR    = 20; % Pom�r sign�l �um v dB
sinus  = A * cos(2*pi*f*t); % cosinusovka
signal = awgn(sinus,SNR); % Add white Gaussian noise to a signal.


indx             = find(sourcemodel_atlas.brick0==oblast); % V�b�r v�ech zdroj� podle oblasti
zdroj_pozice     = randi(length(indx)); % V�b�r konkr�tn�ho zdroje z oblasti n�hodn�
zdroj            = indx(zdroj_pozice); % Index zdroje na z�klad� n�hodn�ho v�b�ru
zdroj_leadfield  = leadfield.leadfield(zdroj); % Pro NUDZ elektrody
zdroj_leadfield  = zdroj_leadfield{1};
signal_leadfield = cell(1,length(signal)); % Definov�n� prom�nn� typu cell

maska = [1 1 0]; % V jak�ch os�ch jde sign�l [x y z];

potencial = zeros(12,length(signal)); % Alokace pam�ti pro �asovou �adu ne * nt

for i=1:length(signal)
    signal_leadfield{i} = signal(i)*zdroj_leadfield; % Pron�soben� vzork� sign�lu s leadfield
    potencial(:,i)      = maska*(signal_leadfield{i})'; % Pron�soben� maskou pro ur�en� moment�
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
cfg.method            = 'dics';
cfg.frequency         = 35;  
cfg.grid              = sourcemodel; 
cfg.headmodel         = headmodel;
cfg.dics.projectnoise = 'yes';
cfg.dics.lambda       = 0;
cfg.elec              = elec;

sourcePost_nocon       = ft_sourceanalysis(cfg, freq);

%% Interpolace ze sourcemodel na MRI (volume)

cfg = [];
cfg.interpmethod = 'nearest';
cfg.parameter = 'pow';
cfg.downsample = 4;

sourceInt = ft_sourceinterpolate(cfg,sourcePost_nocon,mri);

%% Zobrazeni

cfg = [];
cfg.funparameter = 'pow';
ft_sourceplot(cfg,sourceInt)





