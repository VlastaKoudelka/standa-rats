%% read MRI and create mesh
clear all; close all; ft_defaults

filePath = 'D:\_PROJECTS\Potkani GACR\Inputs\NIfTILow\';
file = 'volume_Brain.nii';

%% read MRI
mri = ft_read_mri([filePath file]);

cfg = [];
ft_sourceplot(cfg,mri)

%% automatic MRI segmentation 
cfg =[];
cfg.downsample = 2;
cfg.template = [filePath file];
cfg.output = 'brain';

segMri = ft_volumesegment(cfg,mri);

%% manual MRI segmentation 
brain = zeros(size(mri.anatomy));
brain(mri.anatomy>0) = 1;
brain = logical(brain);

segManMri = mri;
segManMri = rmfield(segManMri,'anatomy');
segManMri.unit = 'mm';
segManMri.brain = brain;

%% mesh
cfg = [];
cfg.method = 'hexahedral';
cfg.downsample = 2;

mesh = ft_prepare_mesh(cfg,segManMri);

ft_plot_mesh(mesh, 'surfaceonly', 'yes');

save Mesh_2400k mesh