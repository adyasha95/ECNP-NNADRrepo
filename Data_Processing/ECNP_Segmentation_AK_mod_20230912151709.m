%Todo: Please adapt the following path to your directory structure
your_path = 'D:\Forschungsmodul - Praezisionspsychiatrie\SummerschoolData\Depression Project\Resting State MRI\MRI-Images\ECNP/';


%%%======================%%%
%%% Open SPM12 and CAT12 %%%
%%%======================%%%

spm_path = [your_path 'spm12/'];
addpath(spm_path)

%Open SPM12
spm pet

%start CAT12 in expert mode
cat12('expert')

%%%===================%%%
%%% Define data paths %%%
%%%===================%%%

data_path = ['D:\Forschungsmodul - Praezisionspsychiatrie\SummerschoolData\Depression Project\Resting State MRI\MRI-Images\Structural MRI/'];

data_paths = readtable(['D:\Forschungsmodul - Praezisionspsychiatrie\SummerschoolData\Depression Project\Resting State MRI\MRI-Images\ECNP\Preprocessing_sMRI.xlsx']);%%here
data_paths = table2cell(data_paths);
data_paths = data_paths(:,2);

%%%==============%%%
%%% Segmentation %%%
%%%==============%%%

clear matlabbatch

% Volumes: highres raw data (e.g. T1 images), either nifti or compressed nifti
matlabbatch{1}.spm.tools.cat.estwrite.data = data_paths;

% hidden in GUI
matlabbatch{1}.spm.tools.cat.estwrite.data_wmh = {''};

%split job into separate processes (multithreading)
matlabbatch{1}.spm.tools.cat.estwrite.nproc = 0;

%% Options for initial SPM processing
%======================================

% hidden in GUI
matlabbatch{1}.spm.tools.cat.estwrite.useprior = '';

% Tissue Probability Map (commented out in standalone)
matlabbatch{1}.spm.tools.cat.estwrite.opts.tpm = {[your_path 'spm12/tpm/TPM.nii']};
matlabbatch{1}.spm.tools.cat.estwrite.opts.affreg = 'mni';
matlabbatch{1}.spm.tools.cat.estwrite.opts.biasstr = 0.5;
% SPM processing accuracy: 0.5: average (in most images good enough)
matlabbatch{1}.spm.tools.cat.estwrite.opts.accstr = 0.5;

%% Extended options for CAT12 preprocessing
%===========================================
%% Segmentation options
matlabbatch{1}.spm.tools.cat.estwrite.extopts.segmentation.restypes.optimal = [1 0.3];
matlabbatch{1}.spm.tools.cat.estwrite.extopts.segmentation.setCOM = 1;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.segmentation.APP = 1070;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.segmentation.affmod = 0;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.segmentation.NCstr = -Inf;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.segmentation.spm_kamap = 0;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.segmentation.LASstr = 0.5;
% Strength of Local Adaptive Segmentation based on the assumption of an equally thick cortex (in development)
%matlabbatch{1}.spm.tools.cat.estwrite.extopts.segmentation.LASmyostr = 0;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.segmentation.gcutstr = 2;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.segmentation.cleanupstr = 0.5;
% hidden in GUI
matlabbatch{1}.spm.tools.cat.estwrite.extopts.segmentation.BVCstr = 0.5;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.segmentation.WMHC = 2;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.segmentation.SLC = 0;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.segmentation.mrf = 1;

%% Spatial registration options
% Shooting registration selected by default
%matlabbatch{1}.spm.tools.cat.estwrite.extopts.registration.regmethod.shooting.shootingtpm = {[your_path 'spm12/toolbox/cat12/templates_MNI152NLin2009cAsym/Template_0_GS.nii']};
% Optimized shoothing - standard (0.5)
%matlabbatch{1}.spm.tools.cat.estwrite.extopts.registration.regmethod.shooting.regstr = 0.5;
% DARTEL registration
matlabbatch{1}.spm.tools.cat.estwrite.extopts.registration.regmethod.dartel.darteltpm = {[your_path 'spm12/toolbox/cat12/templates_MNI152NLin2009cAsym/Template_1_Dartel.nii']};

matlabbatch{1}.spm.tools.cat.estwrite.extopts.registration.vox = 1.5;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.registration.bb = 12;

%% Surface options
matlabbatch{1}.spm.tools.cat.estwrite.extopts.surface.pbtres = 0.5;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.surface.pbtmethod = 'pbt2x';
matlabbatch{1}.spm.tools.cat.estwrite.extopts.surface.SRP = 22;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.surface.reduce_mesh = 1;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.surface.vdist = 2;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.surface.scale_cortex = 0.7;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.surface.add_parahipp = 0.1;
% Initial morphological closing of parahippocampus (to minimize risk of large cuts of 
% parahippocampal gyrus after topology correction)
% may lead to poorer quality of topology correction for other data -> use only if large cuts in 
% parahippocampal areas occur
matlabbatch{1}.spm.tools.cat.estwrite.extopts.surface.close_parahipp = 0; % default 1
% hidden in GUI
matlabbatch{1}.spm.tools.cat.estwrite.extopts.admin.experimental = 0;
% hidden in GUI
matlabbatch{1}.spm.tools.cat.estwrite.extopts.admin.new_release = 0;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.admin.lazy = 0;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.admin.ignoreErrors = 1;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.admin.verb = 2;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.admin.print = 2;

%% Writing options
%===================
% Use BIDS directory structure? -> No selected
matlabbatch{1}.spm.tools.cat.estwrite.output.BIDS.BIDSno = 1;

%% Surface and thickness estimation
matlabbatch{1}.spm.tools.cat.estwrite.output.surface = 0;

% hidden by GUI
matlabbatch{1}.spm.tools.cat.estwrite.output.surf_measures = 1;

%% Process volume ROIs: ROI atlas maps
matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.neuromorphometrics = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.lpba40 = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.cobra = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.hammers = 1;
matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.thalamus = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.ibsr = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.aal3 = 1;
matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.mori = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.anatomy3 = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.julichbrain = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.Schaefer2018_100Parcels_17Networks_order = 1;
matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.Schaefer2018_200Parcels_17Networks_order = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.Schaefer2018_400Parcels_17Networks_order = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.Schaefer2018_600Parcels_17Networks_order = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.ownatlas = {''};

%% Writing options (see cat_defaults for the description of parameters)
%   native    0/1     (none/yes)
%   warped    0/1     (none/yes)
%   mod       0/1/2/3 (none/affine+nonlinear/nonlinear only/both)
%   dartel    0/1/2/3 (none/rigid/affine/both)

%% Options to save grey matter images
% The native space option allows you to save a tissue class image (p*) that is in alignment 
% with the original image
matlabbatch{1}.spm.tools.cat.estwrite.output.GM.native = 1; % default 0
% Write image in normlized space without any modulation.
matlabbatch{1}.spm.tools.cat.estwrite.output.GM.warped = 1; % default 0
matlabbatch{1}.spm.tools.cat.estwrite.output.GM.mod = 1;
matlabbatch{1}.spm.tools.cat.estwrite.output.GM.dartel = 0;

%% Options to save white matter images
% The native space option allows you to save a tissue class image (p*) that is in alignment 
% with the original image
matlabbatch{1}.spm.tools.cat.estwrite.output.WM.native = 0; % default 0
% Write image in normlized space without any modulation.
matlabbatch{1}.spm.tools.cat.estwrite.output.WM.warped = 0; % default 0
matlabbatch{1}.spm.tools.cat.estwrite.output.WM.mod = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.WM.dartel = 0;

%% Options to save cerebro-spinal fluid images
% The native space option allows you to save a tissue class image (p*) that is in alignment 
% with the original image
matlabbatch{1}.spm.tools.cat.estwrite.output.CSF.native = 0; % default 0
% Write image in normlized space without any modulation.
matlabbatch{1}.spm.tools.cat.estwrite.output.CSF.warped = 0; % default 0
% Modulation is to compensate for the effect of spatial normalization 
matlabbatch{1}.spm.tools.cat.estwrite.output.CSF.mod = 0; % default 0
matlabbatch{1}.spm.tools.cat.estwrite.output.CSF.dartel = 0;

%% 'ct' options hidden in GUI
%matlabbatch{1}.spm.tools.cat.estwrite.output.ct.native = 0;
%matlabbatch{1}.spm.tools.cat.estwrite.output.ct.warped = 0;
%matlabbatch{1}.spm.tools.cat.estwrite.output.ct.dartel = 0;

%% Options to save percentage position maps (experimental)
%matlabbatch{1}.spm.tools.cat.estwrite.output.pp.native = 0;
%matlabbatch{1}.spm.tools.cat.estwrite.output.pp.warped = 0;
%matlabbatch{1}.spm.tools.cat.estwrite.output.pp.dartel = 0;


%% White matter hyperintensities (WMHs) - in development
matlabbatch{1}.spm.tools.cat.estwrite.output.WMH.native = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.WMH.warped = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.WMH.mod = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.WMH.dartel = 0;

%% Stroke lesions (SLs) - in development
matlabbatch{1}.spm.tools.cat.estwrite.output.SL.native = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.SL.warped = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.SL.mod = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.SL.dartel = 0;

%% Tissue Probability Map Classes (option to save SPM tissue class 4 to 6)
matlabbatch{1}.spm.tools.cat.estwrite.output.TPMC.native = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.TPMC.warped = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.TPMC.mod = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.TPMC.dartel = 0;

%% Atlas label maps: option to save the selected atlas maps in native space
% The native space option allows you to save a tissue class image (p*) that is in alignment 
% with the original image
matlabbatch{1}.spm.tools.cat.estwrite.output.atlas.native = 1; % default 0

%% Partial volume estimation (PVE) label images
matlabbatch{1}.spm.tools.cat.estwrite.output.label.native = 1;
% Write image in normlized space without any modulation.
matlabbatch{1}.spm.tools.cat.estwrite.output.label.warped = 1; % default 0
matlabbatch{1}.spm.tools.cat.estwrite.output.label.dartel = 0;

% hidden in GUI
matlabbatch{1}.spm.tools.cat.estwrite.output.labelnative = 1;

%% Bias, noise and global intensity corrected T1 image
% The native space option allows you to save a tissue class image (p*) that is in alignment 
% with the original image
matlabbatch{1}.spm.tools.cat.estwrite.output.bias.native = 1; % default 0
matlabbatch{1}.spm.tools.cat.estwrite.output.bias.warped = 1;
matlabbatch{1}.spm.tools.cat.estwrite.output.bias.dartel = 0;

%% Bias, noise and local intensity corrected T1 image
% The native space option allows you to save a tissue class image (p*) that is in alignment 
% with the original image
matlabbatch{1}.spm.tools.cat.estwrite.output.las.native = 1; % default 0
% Write image in normlized space without any modulation.
matlabbatch{1}.spm.tools.cat.estwrite.output.las.warped = 1; % default 0
matlabbatch{1}.spm.tools.cat.estwrite.output.las.dartel = 0;

%% Option to save Jacobian determinant, which exppresses local volume changes
matlabbatch{1}.spm.tools.cat.estwrite.output.jacobianwarped = 1; % default 0

%% Deformation fields 
% Deformation fields can be saved and used by the Deformation Utility and/or applied to 
% coregistered data from other modalities (e.g. fMRI). For spatially normalising images to 
% MNI , you will need the forward deformation, whereas for spatially normalizing (eg) GIFTI 
% surface files, you?ll need the inverse. It is also possible to transform data to MNI space on 
% the individual subject, which also requires inverse transform.
% [1 0]: Image->Template, [0 1]: Template->Image, [1 1]: Inverse + forward
matlabbatch{1}.spm.tools.cat.estwrite.output.warps = [1 1]; % default [1 0]

% Registration matrix: 
% Deformation matrixes (affine and rigid) can be saved and used by the SPM Reorient 
% Images Utility and/or applied to coregistered data from other modalities  (e.g. fMRI). 
% (see previous comment)
% Transformation saved as .mat files (containing transformation matrix)
matlabbatch{1}.spm.tools.cat.estwrite.output.rmat = 1; % default 0

spm_jobman('run',matlabbatch);