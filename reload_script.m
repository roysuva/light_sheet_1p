
%%
% ##################### Reload analyzed data script #######################
%           CodeBase: CNMFE - Large_data - 1p 
% #########################################################################

%% Set up utils and data path 
 
 
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

 
 
%% Load saved data (should load 'neuron' object as a structure)

datanam = 'Data_200001';
initpath_cur = initpath; 
paramspath = fullfile([initpath 'analyzed_data/Data_200001/Data_200001_extracted_params.mat']);
dpath = fullfile(initpath 'analyzed_data/Data_200001/Data_200001_processed_source_extraction/frames_1_5000/LOGS_02-Nov_12_21_02']);  
tformpath_for_bw = [];

load(paramspath); 
neuron = Sources2D();
neuron = loadall_to_class(neuron, dpath);

if ~strcmpi(initpath, initpath_cur)
    initpath = initpath_cur; 
end

%% Load stimulus information 

if ~exist('stim', 'var') || (exist('stim', 'var') && ~isfield(stim,'time_stamps'))
    stim = get_stimprops(stim, datanam, initpath);
end

% check for movie file for RN stim 
if strcmpi(stim.type,'RN')
    if ~isfield(stim,'movie_file') || (isfield(stim,'movie_file') && isempty(stim.movie_file)) 
        [ffile, fpath] = uigetfile('*.xml','Select the movie file for this RN stim'); 
        stim.movie_file = strcat(fpath, ffile); 
    end
end

%% Align frames to stim 

[stim2frame, stim] = get_stim2frame_association(stim, fps); 


%% Get response properties: for eg. direction selectivity for MB or receptive field for RN or response traces and PCAs for CH 

[resp, StimComb, stim] = get_resp_properties(stim, neuron, stim2frame, 'detrend',true); 


%% Generate and save figures 

if strcmpi(stim.type,'MB')
    
    % Fig: Example detrended trace + filtered detrended trace of a DS cell 
    ds_ind = 1; 
    hf = figure; set(gcf,'position',[146 390 1735 360],'color','w');
    plot(resp.t_axis, resp.detrended_trace(resp.ds_ids(ds_ind),:),'-k'); hold on; 
    plot(resp.t_axis, resp.detrended_filt_trace(resp.ds_ids(ds_ind),:),'-r');
    legend_curate({'Detrended trace','Detrended filtered trace'},'markershow',false,...
        'lineshow',false,'text_color',{'k','r'},'box',false); 
    xlabel('Time'); ylabel('Calcium signal'); title(sprintf('DS cell id %d; DSI: %2.3f',resp.ds_ids(ds_ind),...
        resp.ds_estimates.mag(resp.ds_ids(ds_ind))));
    set(gca,'fontsize',10); 
    
    % Fig: DS classification
    hf = figure;  set(gcf,'color','w');
    plot(1:length(neuron.ids), resp.ds_estimates.mag,'ok','markerfacecolor','k'); hold on;
    plot(resp.ds_ids, resp.ds_estimates.mag(resp.ds_ids),'ok','markerfacecolor','r'); hold on;
    xlabel('Cell #'); ylabel('DSI'); set(gca,'xlim',[-1 inf]);
    legend_curate({'All RGCs','DS RGCs'},'markershow',false,...
        'lineshow',false,'text_color',{'k','r'},'box',false,'font_size',10);

    % Fig: Preferred directions figure
    hf = plot_pref_dirs(resp, StimComb);

    % Fig: Traces along trials with polar plot at the center
    plot_ds_polar_mb(stim, StimComb, resp, matfigsavepath,  false);


elseif strcmpi(stim.type,'CH')

 
    % Cluster comparisons
    cluster_GMM = cluster_traces(resp, neuron,'cluster_method','GMM','numclust',8,'which_trials','all');
    cluster_HAC = cluster_traces(resp, neuron,'cluster_method','Hierarchical','numclust',8,'which_trials','all');
    cluster_SC = cluster_traces(resp, neuron,'cluster_method','Spectral','numclust',8,'which_trials','all');


    % Gather all normalized psth (exclude outliers) in random order
    outliers = find(isnan(cluster_GMM.clust_id));
    include = setxor(1:numel(resp.cell),outliers);
    tracemat = [];
    for k=1:numel(include)
        tracemat = [tracemat; resp.trace_psth(:,k)'./norm(resp.trace_psth(:,k))];
    end
    tracemat = tracemat(randperm(length(include),length(include)),:);


    % Gather silhouette coefficients for N runs

    sil_gather_GMM = zeros(cluster_GMM.silhouette_stat.niter, ...
        numel(cluster_GMM.silhouette_stat.evals{1}.CriterionValues));
    sil_clust_rng = cluster_GMM.silhouette_stat.evals{1}.InspectedK;
    for ni=1:cluster_GMM.silhouette_stat.niter
        sil_gather_GMM(ni,:)  = cluster_GMM.silhouette_stat.evals{3}.CriterionValues;
    end
    sil_gather_GMM_med = median(sil_gather_GMM,1);
    sil_gather_GMM_mad = mad(sil_gather_GMM,[],1);
    idnan = find(isnan(sil_gather_GMM_med));
    sil_gather_GMM_med(idnan) = mean([sil_gather_GMM_med(idnan-1) sil_gather_GMM_med(idnan+1)]);
    sil_gather_GMM_mad(idnan) = mean([sil_gather_GMM_mad(idnan-1) sil_gather_GMM_mad(idnan+1)]);


    sil_gather_HAC = zeros(cluster_HAC.silhouette_stat.niter, ...
        numel(cluster_HAC.silhouette_stat.evals{1}.CriterionValues));
    sil_clust_rng = cluster_HAC.silhouette_stat.evals{1}.InspectedK;
    for ni=1:cluster_HAC.silhouette_stat.niter
        sil_gather_HAC(ni,:)  = cluster_HAC.silhouette_stat.evals{3}.CriterionValues;
    end
    sil_gather_HAC_med = median(sil_gather_HAC,1);
    sil_gather_HAC_mad = mad(sil_gather_HAC,[],1);
    idnan = find(isnan(sil_gather_HAC_med));
    sil_gather_HAC_med(idnan) = mean([sil_gather_HAC_med(idnan-1) sil_gather_HAC_med(idnan+1)]);
    sil_gather_HAC_mad(idnan) = mean([sil_gather_HAC_mad(idnan-1) sil_gather_HAC_mad(idnan+1)]);


    sil_gather_SC = zeros(cluster_SC.silhouette_stat.niter, ...
        numel(cluster_SC.silhouette_stat.evals{1}.CriterionValues));
    sil_clust_rng = cluster_SC.silhouette_stat.evals{1}.InspectedK;
    for ni=1:cluster_SC.silhouette_stat.niter
        sil_gather_SC(ni,:)  = cluster_SC.silhouette_stat.evals{3}.CriterionValues;
    end
    sil_gather_SC_med = median(sil_gather_SC,1);
    sil_gather_SC_mad = mad(sil_gather_SC,[],1);
    idnan = find(isnan(sil_gather_GMM_med));
    sil_gather_SC_med(idnan) = mean([sil_gather_SC_med(idnan-1) sil_gather_SC_med(idnan+1)]);
    sil_gather_SC_mad(idnan) = mean([sil_gather_SC_mad(idnan-1) sil_gather_SC_mad(idnan+1)]);

end



%% Generate spatial and temporal images of example neurons 

d1 = neuron.options.d1;  
d2 = neuron.options.d2;  
d3 = size(neuron.C,2); 

 

% Spatial footprint 
for rc=1
    ln_cl = 'm';
    txt_cl = 'y';
    ln_wid = 1; 
    ind_show = 1:size(neuron.C,1);
    hf102 = fig_pintomonitor('','fracx',0.9,'aspect_ratio_x_y',1);
    ha102 = axes('Parent',hf102);
    clrmap = [linspace(0,1,100)' linspace(0,1,100)' linspace(0,1,100)'];
    if isempty(medZim)
        medZim = median(double(im_uint16_red(:,:,FR(1):FR(2))),3); 
    end
    show_rois = 'all'; % 'all', or indices, or 'none' [8 10 11 13 24]
    [Coor] = contour_plot_simple(neuron.A,medZim,neuron,true,ln_cl,ln_wid,txt_cl,clrmap,ha102,show_rois);
    title('Spatial components'); imcontrast;
end



% All temporal traces (heatmap)
rng_rec = [prctile(tracemat(:),2) prctile(tracemat(:),99)];
hf = fig_pintomonitor('','fracx',0.5,'aspect_ratio_x_y',4/3);
set(hf,'color','w', 'renderer','painters');
imshow(tracemat,rng_rec); colormap parula; colorbar;

 
