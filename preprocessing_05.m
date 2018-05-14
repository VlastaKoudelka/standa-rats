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
load('mri4down'); % podvzorkované nesegmentované mri mozku potkana

%% Generování cosinusovky se SNR
oblast = 11; % index oblasti (11...hippocampus)
cela   = 1; % 1 pro celou oblast/ 0 pro 1 zdroj
moment = [0 0 1]; % V jakých osách jde signál [x y z];

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
cfg.parameter = 'pow'; % power spektrum jako interpolované hodnoty
sourceInt = ft_sourceinterpolate(cfg,sourcePost_nocon,atlas4down);
sourceInt.unit = 'mm';
%% Parcelace zdrojù podle oblastí atlasu

cfg = [];
cfg.method = 'mean'; % støední hodnota
cfg.parcellation = 'brick0';
parcel = ft_sourceparcellate(cfg, sourceInt, atlas4down); % Parcelace zdrojù
sourceInt.parcel = zeros(length(sourceInt.pow),1); % Vytvoøení nové promìnné pro uložení parcelovaných zdrojù
for i=1:length(parcel.pow)
      sourceInt.parcel(find(atlas4down.brick0==i))=parcel.pow(i); % Pøiøazení støední hodnoty oblasti do všech voxelù oblasti
end

%% Skutecne zdroje

cfg = [];
cfg.method = 'mean'; % støední hodnota
cfg.parcellation = 'brick0';
parcel = ft_sourceparcellate(cfg, sourceInt, atlas4down); % Parcelace zdrojù
sourceInt.trueSource = zeros(length(sourceInt.pow),1); % Vytvoøení nové promìnné pro uložení umelych zdroju

sourceInt.trueSource(find(atlas4down.brick0==oblast))=1; % Pøiøazení støední hodnoty oblasti do všech voxelù oblasti

%% Maska

cfg = [];
cfg.inputcoord = 'unknown';
cfg.atlas = atlas4down;
cfg.roi = atlas4down.brick0label(:);
mask = ft_volumelookup(cfg, sourceInt); %Vytvoøení masky pro zobrazení zdrojù pouze v oblasti mozku

sourceInt.mask = mask;


%% Zobrazení
cfg = [];
cfg.method = 'ortho';
cfg.funparameter = 'parcel'; % parcel-rozparcelované zdroje/pow-originál interpolované
cfg.maskparameter = 'mask'; % maska
cfg.atlas = atlas4down;
%cfg.mesh = mesh;
ft_sourceplot(cfg,sourceInt)

%% Zobrazení
cfg = [];
cfg.method = 'ortho';
cfg.funparameter = 'trueSource'; % parcel-rozparcelované zdroje/pow-originál interpolované
cfg.maskparameter = 'mask'; % maska
cfg.atlas = atlas4down;
%cfg.mesh = mesh;
ft_sourceplot(cfg,sourceInt)




