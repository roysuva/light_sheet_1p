
%%
% ########################### Analysis script #############################
%           CodeBase: CNMFE - Large_data - 1p 
% #########################################################################

%% Set up utils and data path 
 
clear;

 
carch = computer('arch'); 

% Set path directives
initpath = '/Users/username/light_sheet_1p/';  
 
cnmfepath = [initpath 'CNMF_E/']; 
normcorrepath = [initpath 'NoRMCorre-master/'];
curfold = cd(cnmfepath);
cnmfe_setup(cnmfepath);
cd(curfold);
addpath(genpath(normcorrepath));
addpath(genpath([initpath 'Internal_utils/']));
addpath(genpath([initpath 'External_utils/']));


% Set full path for raw data, scripts and analyzed data and figures
datanam = 'Data_200001';
tiffpath = [initpath 'example_data/'];
fijiprojim = [];
if ~exist([initpath 'Matlab_outputs'],'dir')
    mkdir([initpath 'Matlab_outputs']);
end
savepath = [initpath 'analyzed_data/'];
 
% Set path for notebooks containing info about experiments and stimulus
notebookpath = [initpath 'example_data/'];

% Assign omnipresent variables
global fps bitdepth;

% Do you have stimulus?
stim.exist = true;
stim.source = 'DLP'; % 'LED', 'DLP'
stim.type = 'CH'; % CH: Chirp 
stim.filepath = [initpath 'example_data/stimuli/'];  
stim.movie_file = [];
stim.repeat = true;
 

% Do you want to save videos?
make_video = false;
 
 
 
%% Set up parallel workers, gpu arrays, located mex files and screen preference 
 
% Setup parallel workers 
v = ver; 
popt = any(strcmp({v.Name}, 'Parallel Computing Toolbox')); 
if popt 
    pp = gcp('nocreate'); 
    ncores = maxNumCompThreads('automatic'); 
    if ncores > 1 && isempty(pp)
        parpool('local', ncores); 
    end
end
 
% Set up gpu flag (if it exists)
if ~gpuDeviceCount
    gpuflag = true;  
else
    gpuflag = false;  
end
 


%% Read image file, remove outliers, denoise images 



% ------------------------------ Step 1 -----------------------------------
% Read image properties
tiffdirinfo = dir(fullfile(tiffpath, '*.tif'));
chunklim = zeros(length(tiffdirinfo),2);
for fln = 1:length(tiffdirinfo)
    info = imfinfo([tiffpath,tiffdirinfo(fln).name]);
    if fln==1
        chunklim(fln,1:2) = [1 numel(info)];
    else
        chunklim(fln,1:2) = chunklim(fln-1,2)+[1 numel(info)];
    end
end
num_images = chunklim(end,end);
prompt = {'Total # of images: ','Enter range of images to analyze: '};
dlgtitle = 'Number of images';
dims = [1 50]; % size of dialogue box in pixels
definput = {num2str(num_images),sprintf('%d - %d',1,num_images)};
answer = inputdlg(prompt,dlgtitle,dims,definput);
image_rng = cellfun(@str2double,strsplit(answer{2},'-'));
num_images = diff(image_rng)+1;
bitdepth = info(1).BitDepth;
dims = [info(1).Height info(2).Width];

strextrct = info(1).ImageDescription(regexpi(info(1).ImageDescription,'Exposure'):regexpi(info(1).ImageDescription,'Exposure')+18);
fps = 1/str2double(strextrct(regexpi(strextrct,'=')+2:regexpi(strextrct,'=')+8)); % frame rate (/sec)
decaytconst = 300/1000; % sec (GCaMP7f: from Dana et.al 2018)

% ------------------------------ Step 2 -----------------------------------
% Calcium image files can be big, we want to minimize redundancy in the
% image set for analysis as much as possible
clm = find(chunklim(:,1)<=image_rng(1),1,'last'):find(chunklim(:,2)>=image_rng(2),1,'first');
im_uint16_temp = cell(1,length(clm));
parfor fln=1:length(clm)
    if image_rng(1)>=chunklim(clm(fln),1)
        init_fr = image_rng(1)-chunklim(clm(fln),1)+1;
    else
        init_fr = 1;
    end
    if image_rng(2)<=chunklim(clm(fln),2) && fln==1
        num_images_temp = diff(image_rng)+1;
    else
        num_images_temp = min([chunklim(clm(fln),2) image_rng(2)]) - max([chunklim(clm(fln),1) image_rng(1)]) + 1;
    end
    [im_uint16_temp{fln},~,~] = smod_bigread2(fullfile(tiffpath,tiffdirinfo(fln).name), init_fr, num_images_temp);
end
im_uint16 = cat(3,im_uint16_temp{:});
clear im_uint16_temp;


% Set range of median estimation (max frames = 3000)
FR=[1 5000];
if size(im_uint16,3)<diff(FR)+1
    FR=[1 size(im_uint16,3)];
end

% Get std and median projections
medZim = median(im_uint16(:,:,FR(1):FR(2)),3);  % Median projection
stdZim = std(double(im_uint16(:,:,FR(1):FR(2))),0,3); % std projection
qtl = 0.75; % upper quantile value
range_med = get_intensity_range(medZim);
range_std = get_intensity_range(stdZim);
hf1 = fig_pintomonitor();
ha11 = subplot(1,2,1);
imshow(medZim, range_med, 'Parent', ha11); title('Median image');
ha12 = subplot(1,2,2);
imshow(stdZim, range_std, 'Parent', ha12); title('Std image');
close(hf1);


% Ask user for selecting Area of Interest
if exist(AOI_path,'file')
    replyAOI = questdlg('AOI path provided. Do you want to load AOI from another data?',...
        'AOI selection','Yes','No','No');
    if strcmpi(replyAOI,'Yes')
        importline = true;
    else
        importline = false;
    end
else
    importline = false;
end
if importline
    [AOI_x,AOI_y] = reload_AOI(AOI_path);
    hf2 = figure;
    fig_pintomonitor(hf2,'aspect_ratio_x_y',3/4,'fracx',0.7);
    range_med = get_intensity_range(medZim);
    imshow(medZim, range_med); hold on;
    plot(AOI_x,AOI_y,'-y','linewidth',3); legend('Area for analysis');
else
    hf2 = figure;
    fig_pintomonitor(hf2,'aspect_ratio_x_y',3/4,'fracx',0.7);
    range_med = get_intensity_range(medZim);
    imshow(medZim, range_med); imcontrast;  hold on; title('Select the vertices of rectangular Area of Interest');
    [x,y] = deal([]);
    for i=1:2
        [x(i),y(i)] = ginput(1);
        x(i) = round(x(i));
        y(i) = round(y(i));
        if x(i)<1; x(i)=1; end
        if y(i)<1; y(i)=1; end
        if x(i)>size(medZim,2); x(i)=size(medZim,2); end
        if y(i)>size(medZim,1); y(i)=size(medZim,1); end
        h = plot(x(i),y(i),'or','markersize',8,'linewidth',3); hold on;
        set(get(get(h,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    end
    AOI_x = [x(1) x(2) x(2) x(1) x(1)]; % Area of Interest
    AOI_y = [y(1) y(1) y(2) y(2) y(1)];
    plot(AOI_x,AOI_y,'-y','linewidth',3); legend('Area for analysis');
end
close(hf2);


% Generate a smaller set of image
im_uint16_red = im_uint16(min(AOI_y):max(AOI_y),min(AOI_x):max(AOI_x),:);
dims_red = [size(im_uint16_red,1) size(im_uint16_red,2)];
medZim = median(im_uint16_red(:,:,FR(1):FR(2)),3); % median projection

% ------------------------------ Step 3 -----------------------------------
% Extract radius of putative ROIs - cell bodies, dendritic spines etc.
% This will be used for background brightness adjustment and by CaImAn (as
% the standard deviation of Gaussian kernel for spatial segmentation; std =
% half the size of neuron)
% Limitation: Current version allows only 1 size of ROIs.


hf3 = figure;
fig_pintomonitor(hf3,'aspect_ratio_x_y',3/4,'fracx',0.7);
range_medZim = get_intensity_range(medZim);
imshow(medZim,range_medZim); imcontrast; hold on;
title('Zoom into a region for selecting ROIs, then press Enter!');
zoom(hf3,'on');
pause();
title('Select the outline of an ROI using left click, and press Enter when done!');
zoom reset;
clear x y;
[x_,y_] = deal([]);
but = 1;
while ~isempty(but)
    [x,y,but] = ginput(1);
    x_ = [x_ x];
    y_ = [y_ y];
    plot(x_,y_,'.-m','markersize',8,'linewidth',1); hold on;
end
title('Template ROI is now selected');
pshape = polyshape(x_,y_,'simplify',true);
rad = round(sqrt(pshape.area/pi)); % fit a circle to a polygon
close(hf3);


% ------------------------------ Step 4 -------------------------------
denoiseK = questdlg( 'Do you want to denoise images?','Denoise using Kalman.',...
    'Yes','No','No');
if strcmpi(denoiseK,'Yes')
    im_double_temp = double(im_uint16_red);
    im_double_kalmanfilt = Kalman_Stack_Filter(im_double_temp, 0.75);
    im_uint16_red = uint16(im_double_kalmanfilt);
    delete im_double_temp im_double_kalmanfilt

    medZim = median(im_uint16_red(:,:,FR(1):FR(2)),3); % median projection
    hf4 = figure;
    fig_pintomonitor(hf4,'aspect_ratio_x_y',3/4,'fracx',0.7);
    range_medZim = get_intensity_range(medZim);
    imshow(medZim,range_medZim,'initialmagnification','fit');
    close(hf4);
end




% ------------------------------ Step 5 -----------------------------------
% Save processed images as tiff file
processedim_savepath = fullfile([matfigsavepath,datanam,'_processed.tif']);
if exist(processedim_savepath,'file'); delete(processedim_savepath); end
if exist('tiffObj','var'); close(tiffObj); clear tiffObj; end
clear tagstruct;
tiffObj = Tiff(processedim_savepath,'w8');
tagstruct.Compression = Tiff.Compression.None;
tagstruct.BitsPerSample = info(1).BitsPerSample;
tagstruct.SamplesPerPixel = info(1).SamplesPerPixel;
tagstruct.SampleFormat = Tiff.SampleFormat.UInt;
tagstruct.RowsPerStrip = size(im_uint16_red,1);
tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.(info(1).PlanarConfiguration);
tagstruct.ImageLength = size(im_uint16_red,1);
tagstruct.ImageWidth = size(im_uint16_red,2);
tagstruct.Photometric = Tiff.Photometric.MinIsBlack;
fprintf('Writing processed tiff files\n');
parfor_progress(num_images);
for tf=1:num_images
    tiffObj.setTag(tagstruct);
    tiffObj.write(im_uint16_red(:, :, tf));
    if tf ~= num_images
        tiffObj.writeDirectory();
    end
    parfor_progress();
end
parfor_progress(0);
tiffObj.close();


% ------------------------------ Step 6 -----------------------------------
% Get stimulus info
if stim.exist
    if strcmpi(stim.source,'LED')
        stim = get_stimprops(notebookpath, datanam, datadat, 'num_images', num_images, 'fps', fps);
    elseif strcmpi(stim.source,'DLP')
        stim = get_stimprops(stim, datanam);
    end
else
    stim = [];
end

% ----------------------------- Step 7 --------------------------------
% Clear memory
if exist('im_uint16_red','var')
    clear im_uint16_red;
end
if exist('im_uint16','var')
    clear im_uint16;
end


% ------------------------------ Step 8 ------------------------------
% Save results so far
paramsfile = fullfile([savepath,datanam,'_extracted_params.mat']);
if exist(paramsfile,'file')
    delete(paramsfile);
end
pathvars = who;
for nv=1:length(pathvars)
    if nv==1
        save(paramsfile,pathvars{nv},'-v7.3');
    else
        save(paramsfile,pathvars{nv},'-append');
    end
end





%% choose data and assign parameter values 

neuron = Sources2D();
nam = processedim_savepath;
nam = neuron.select_data(nam);  %if nam is [], then select data interactively



% parameters
% -------------------------    COMPUTATION    -------------------------  %
pars_envs = struct('memory_size_to_use', 120, ...   % GB, memory space you allow to use in MATLAB
    'memory_size_per_patch', 1, ...   % GB, space for loading data within one patch
    'patch_dims', [64, 64]);  %GB, patch size [64, 64]

% -------------------------      SPATIAL      -------------------------  %
gSig = 1.5*rad;         % pixel, gaussian width of a gaussian kernel for filtering the data. 0 means no filtering
gSiz = (2*gSig+1);  % pixel, neuron diameter
ssub = 1;           % spatial downsampling factor
with_dendrites = false;   % with dendrites or not
if with_dendrites
    % determine the search locations by dilating the current neuron shapes
    updateA_search_method = 'dilate';  %#ok<UNRCH>
    updateA_bSiz = 5;
    updateA_dist = neuron.options.dist;
else
    % determine the search locations by selecting a round area
    updateA_search_method = 'ellipse'; %#ok<UNRCH>
    updateA_dist = 5;
    updateA_bSiz = neuron.options.dist;
end
spatial_constraints = struct('connected', true, 'circular', true);  % you can include following constraints: 'circular'
spatial_algorithm = 'hals_thresh';

% -------------------------      TEMPORAL     -------------------------  %
Fs = fps;                                   % frame rate
decaytconst = 300/1000;                     % in sec (GCaMP6f: 300ms, GCaMP7f: 180ms, GCaMP8f: 100ms) 
tsub = 1;                                   % temporal downsampling factor
deconv_flag = true;                         % run deconvolution or not 
deconv_options = struct('type', 'ar1', ...  % model of the calcium traces. {'ar1', 'ar2'}
    'method', 'foopsi', ...                 % method for running deconvolution {'foopsi', 'constrained', 'thresholded'}
    'smin', -2.5, ...                         % minimum spike size. When the value is negative, the actual threshold is abs(smin)*noise level
    'optimize_pars', true, ...              % optimize AR coefficients
    'optimize_b', true, ...                 % optimize the baseline);
    'lambda',0.99,...
    'max_tau', ceil(fps*10*decaytconst));   % maximum decay time (unit: frame); Set approx 10 times the decay constant of Calcium Indicator

nk = 3;                                     % detrending the slow fluctuation. usually 1 is fine (no detrending)
                                            % when changed, try some integers smaller than total_frame/(Fs*30)
detrend_method = 'spline';                  % compute the local minimum as an estimation of trend.

% -------------------------     BACKGROUND    -------------------------  %
bg_model = 'ring';                          % model of the background {'ring', 'svd'(default), 'nmf'}
nb = 1;                                     % number of background sources for each patch (only be used in SVD and NMF model)
bg_neuron_factor = 1.4; 
ring_radius = round(bg_neuron_factor*gSiz); % when the ring model used, it is the radius of the ring used in the background model. (default: 16)
                                            % otherwise, it's just the width of the overlapping area (default: 18)
num_neighbors = [];                         % number of neighbors for each neuron
bg_ssub = 2;                                % downsample background for a faster speed 

% -------------------------      MERGING      -------------------------  %
show_merge = false;     % if true, manually verify the merging step
merge_thr = 0.70;       % thresholds for merging neurons; [spatial overlap ratio, temporal correlation of calcium traces, spike correlation]
method_dist = 'center';    % method for computing neuron distances {'mean', 'max'} : mean uses COM, max uses location of max activity
dmin = 2*gSig;        % minimum distances between two neurons. it is used together with merge_thr
dmin_only = dmin;       % merge neurons if their distances are smaller than dmin_only (default: 2)
merge_thr_spatial = [0.9, 0.7, -inf];  % merge components with highly correlated spatial shapes (corr=0.8) and small temporal correlations (corr=0.1)

% -------------------------  INITIALIZATION   -------------------------  %
K = [];                     % maximum number of neurons per patch. when K=[], take as many as possible.
init_method = 'greedy';     % method of initialization ('greedy','sparse_NMF','both') - sparse_NMF for dendrites/axons 
min_corr = 0.80;             % minimum local correlation for a seeding pixel
min_pnr = 9;               % minimum peak-to-noise ratio for a seeding pixel
min_pixel = gSig;           % minimum number of nonzero pixels for each neuron (default: 5)
bd = 0;                     % number of rows/columns to be ignored in the boundary (mainly for motion corrected data)
frame_range = [];           % when [], uses all frames
save_initialization = false;% save the initialization procedure as a video 
use_parallel = true;        % use parallel computation for parallel computing
show_init = true;           % show initialization results
choose_params = true;       % manually choose parameters
center_psf = true;          % set as true when background fluctuation is large (usually 1p), false when background fluctuation is small (2p)
df_window = []; 


% -------------------------  Residual   -------------------------  %
min_corr_res = 0.75;
min_pnr_res = 9;
seed_method_res = 'auto';   % method for initializing neurons from the residual
update_sn = true;

% ----------------------  WITH MANUAL INTERVENTION  --------------------  %
with_manual_intervention = true;

% -------------------------  FINAL RESULTS   -------------------------  %
show_demixed_video = false; 
save_demixed = false;        % save the demixed file or not
kt = 3;                     % frame intervals

% -------------------------    UPDATE ALL    -------------------------  %
neuron.updateParams('gSig', gSig, ...       % -------- spatial --------
    'gSiz', gSiz, ...
    'ring_radius', ring_radius, ...
    'ssub', ssub, ...
    'search_method', updateA_search_method, ...
    'bSiz', updateA_bSiz, ...
    'dist', updateA_bSiz, ...
    'spatial_constraints', spatial_constraints, ...
    'spatial_algorithm', spatial_algorithm, ...
    'tsub', tsub, ...                       % -------- temporal --------                
    'deconv_flag', deconv_flag, ...
    'deconv_options', deconv_options, ...
    'nk', nk, ...
    'detrend_method', detrend_method, ...
    'background_model', bg_model, ...       % -------- background --------
    'nb', nb, ...
    'ring_radius', ring_radius, ...
    'num_neighbors', num_neighbors, ...
    'bg_ssub', bg_ssub, ...
    'merge_thr', merge_thr, ...             % -------- merging ---------
    'dmin', dmin, ...
    'method_dist', method_dist, ...
    'min_corr', min_corr, ...               % ----- initialization -----
    'init_method', init_method, ...
    'min_pnr', min_pnr, ...
    'min_pixel', min_pixel, ...
    'bd', bd, ...
    'center_psf', center_psf,...
    'df_window', df_window,...              % ----- de-trending dF/F0 -----
    'df_prctile', 20);
neuron.Fs = Fs;


%% Run CNMF-E 

echo on; 

% -------------------------------------------------------------------------
% distribute data and be ready to run source extraction
neuron.getReady(pars_envs);


% -------------------------------------------------------------------------
% initialize neurons from the video data within a selected temporal range
if choose_params
    % change parameters for optimized initialization
    [gSig, gSiz, ring_radius, min_corr, min_pnr] = neuron.set_parameters();
end

[center, Cn, PNR] = neuron.initComponents_parallel(K, frame_range, save_initialization, use_parallel);
neuron.PNR = PNR; 
neuron.compactSpatial();
if show_init
    figure();
    ax_init= axes();
    imagesc(Cn, [0, 1]); colormap gray;
    hold on;
    plot(center(:, 2), center(:, 1), '.r', 'markersize', 10);
end


% -------------------------------------------------------------------------
% estimate the background components
neuron.update_background_parallel(use_parallel);
neuron_init = neuron.copy();

% -------------------------------------------------------------------------
%  merge neurons and update spatial/temporal components
neuron.merge_neurons_dist_corr(show_merge);
neuron.merge_high_corr(show_merge, merge_thr_spatial);

% -------------------------------------------------------------------------
% update spatial components

% -------------------------------------------------------------------------
% pick neurons from the residual
[center_res, Cn_res, PNR_res] =neuron.initComponents_residual_parallel([], save_initialization, use_parallel, min_corr_res, min_pnr_res, seed_method_res);
if show_init
    axes(ax_init);
    plot(center_res(:, 2), center_res(:, 1), '.g', 'markersize', 10);
end
neuron_init_res = neuron.copy();

% -------------------------------------------------------------------------
% update spatial&temporal components, delete false positives and merge neurons
if update_sn
    neuron.update_spatial_parallel(use_parallel, true);
    udpate_sn = false;
else
    neuron.update_spatial_parallel(use_parallel);
end
% merge neurons based on correlations 
neuron.merge_high_corr(show_merge, merge_thr_spatial);

for m=1:2
    % update temporal
    neuron.update_temporal_parallel(use_parallel);
    
    % delete bad neurons
    neuron.remove_false_positives();
    
    % merge neurons based on temporal correlation + distances 
    neuron.merge_neurons_dist_corr(show_merge);
end


% -------------------------------------------------------------------------
% add a manual intervention and run the whole procedure for a second time
if with_manual_intervention
    show_merge = true;
    neuron.orderROIs('snr');   % order neurons in different ways {'snr', 'decay_time', 'mean', 'circularity'}
    neuron.viewNeurons([], neuron.C_raw);
    
    % merge closeby neurons
    neuron.merge_close_neighbors(true, dmin_only); 
    
    % delete neurons
    tags = neuron.tag_neurons_parallel();  % find neurons with fewer nonzero pixels than min_pixel and silent calcium transients
    ids = find(tags>0); 
    if ~isempty(ids)
        neuron.viewNeurons(ids, neuron.C_raw); 
    end
end

% -------------------------------------------------------------------------
% run more iterations
neuron.update_background_parallel(use_parallel);
neuron.update_spatial_parallel(use_parallel);
neuron.update_temporal_parallel(use_parallel);

K = size(neuron.A,2);
tags = neuron.tag_neurons_parallel();  % find neurons with fewer nonzero pixels than min_pixel and silent calcium transients
neuron.remove_false_positives();
neuron.merge_neurons_dist_corr(show_merge);
neuron.merge_high_corr(show_merge, merge_thr_spatial);

if K~=size(neuron.A,2)
    neuron.update_spatial_parallel(use_parallel);
    neuron.update_temporal_parallel(use_parallel);
    neuron.remove_false_positives();
end

% -------------------------------------------------------------------------
% final check for overlapping neurons 
neuron.merge_neurons_dist_corr(show_merge);
neuron.merge_close_neighbors(true, dmin_only);
neuron.merge_high_corr(show_merge, merge_thr_spatial);
if K~=size(neuron.A,2)
    neuron.remove_false_positives();
    neuron.update_spatial_parallel(use_parallel);
    neuron.update_temporal_parallel(use_parallel);
end

% -------------------------------------------------------------------------
% save the workspace for future analysis
neuron.orderROIs('snr');

% -------------------------------------------------------------------------
% show neuron contours
Coor = neuron.show_contours(0.9);

% -------------------------------------------------------------------------
% save workspace and neuron obj as a structure 
cnmfe_path = neuron.saveall_to_struct(); 
 

%% Load stimulus information 

if ~exist('stim', 'var') || (exist('stim', 'var') && ~isfield(stim,'time_stamps'))
    stim = get_stimprops(stim, datanam, initpath);
end

%% Align frames to stim 

[stim2frame, stim] = get_stim2frame_association(stim, fps); 


%% Get response properties: for eg. direction selectivity for MB or receptive field for RN or response traces and PCAs for CH 

[resp, StimComb, stim] = get_resp_properties(stim, neuron, stim2frame, 'detrend',true); 


