clc;clear;close all;
ft_defaults;
%% Naètení potøebných souborù
load('mesh4down'); % mesh
load('headmodel4down'); % headmodel
load('elec'); % NUDZ elektrody
load('gridinside'); % Grid uvnitø mozku
load('SRC/CBWJ13_P80_indexed_volume/lut.mat'); % lut pro atlas
load('atlas4down'); % atlas podvzorkovaný 4x
load('sourcemodel'); % sourcemodel
load('sourcemodel_atlas'); % sourcemodel + atlas = indexované zdroje
load('leadfield'); % leadfield pro všechny zdroje (NUDZ elektrody)
load('mri_segmented'); % mri mozku potkana segmentované
load('mri'); % nesegmentované mri mozku potkana

atlas4down.unit='mm';

%% Generování cosinusovky se SNR
oblast = 11; % index oblasti (11...hippocampus)

fvz    = 1000; % vzorkovací frekvence
f      = 35; % frekvence signálu
T      = 60; % délka signálu
t      = 0:1/fvz:T-1/fvz; % èasový vektor
A      = 1; % amplituda
SNR    = 20; % Pomìr signál šum v dB
sinus  = A * cos(2*pi*f*t); % cosinusovka
signal = awgn(sinus,SNR); % Add white Gaussian noise to a signal.


indx             = find(sourcemodel_atlas.brick0==oblast); % Výbìr všech zdrojù podle oblasti
zdroj_pozice     = randi(length(indx)); % Výbìr konkrétního zdroje z oblasti náhodnì
zdroj            = indx(zdroj_pozice); % Index zdroje na základì náhodného výbìru
zdroj_leadfield  = leadfield.leadfield(zdroj); % Pro NUDZ elektrody
zdroj_leadfield  = zdroj_leadfield{1};
signal_leadfield = cell(1,length(signal)); % Definování promìnné typu cell

maska = [1 1 0]; % V jakých osách jde signál [x y z];

potencial = zeros(12,length(signal)); % Alokace pamìti pro èasovou øadu ne * nt

for i=1:length(signal)
    signal_leadfield{i} = signal(i)*zdroj_leadfield; % Pronásobení vzorkù signálu s leadfield
    potencial(:,i)      = maska*(signal_leadfield{i})'; % Pronásobení maskou pro urèení momentù
end

%% Fieldtrip format
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





