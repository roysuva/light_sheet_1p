%%
%--------------------------------------------------------------------------
% Simple script for extracting dF/F0 for manually selected ROIs
% Data: 2022-07-28-0/data_10 (Ai148;PVCre)
%--------------------------------------------------------------------------
%% 

% Load raw calcium images 
fs = filesep; 
if regexpi(carch,'win')
    initpath = ['Z:' fs 'lab' fs];
else
    initpath = [fs 'Volumes' fs 'dusom_fieldlab' fs 'All_Staff' fs 'lab' fs];  
end
dpath = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Imaging/Light_sheet/Analysis/2022-07-28-0/LE_ventral/Data_1000001/Data_1000001.tif'; 

[im_uint16_temp,~,~] = smod_bigread2(fullfile([initpath, 'Experiments/Imaging/Light_sheet/Analysis/2022-07-28-0/LE_ventral/Data_1000001/Data_1000001.tif']));

% generate some figures 
median_img_new = median(im_uint16_temp, 3); 
rng_med = [prctile(median_img_new(:),2) prctile(median_img_new(:),99)];
std_img_new = std(double(im_uint16_temp), 0, 3);
rng_std = [prctile(std_img_new(:),2) prctile(std_img_new(:),99)];

std_img_filtered = medfilt2(std_img_new,[2 2]); 
rng_std_filtered = [prctile(std_img_filtered(:),2) prctile(std_img_filtered(:),99)];


hf1 = figure('color','w','renderer','painters'); 
imshow(median_img_new,rng_med,'InitialMagnification','fit'); axis equal; hold on;
imcontrast; 

hf2 = figure('color','w','renderer','painters'); 
imshow(std_img_new,[25.1007 70.5116],'InitialMagnification','fit'); axis equal; hold on;
imcontrast; 

hf3 = figure('color','w','renderer','painters'); 
imshow(std_img_filtered,rng_std_filtered,'InitialMagnification','fit'); axis equal; hold on;
imcontrast; 

% fps and nrepeats
fps = 5; % 200ms exposure 
nrepeats = 6; 

% Now estimate dF/F0 trace for freehand selected region 
hf2 = figure('color','w','renderer','painters','position',[210  17  1394 952]); 
imshow(std_img_new,[25.1007 70.5116],'InitialMagnification','fit'); axis equal; hold on;
sz = size(im_uint16_temp); 
boundary_x = [1 sz(2) sz(2) 1 1]; 
boundary_y = [1 1 sz(1) sz(1) 1];  
imcontrast; 
hold on; 
plot(boundary_x, boundary_y,'-r'); 
contained_px = deal(cell(1,5)); 
roi = struct(); 
nrois = 5; 
for i=1:nrois
    dfh = drawfreehand(gca);
    roi(i).position = dfh.Position; 
    mask = createMask(dfh); 
    contained_px{i} = find(mask==1); 
    [rc,cc] = ind2sub(sz(1:2), contained_px{i}); 
    for k=1:size(im_uint16_temp, 3)
        tempA = median(im_uint16_temp(rc, cc,k)); 
        roi(i).f(k) = double(median(tempA(:))); 
    end
end

% calculate moving baseline, dF/F0
for i=1:nrois
    df_prctile = 10; % use low fluorescence values (excludes stimulus induced transients)
    df_twind = 20; % seconds
    df_fwind = df_twind*fps;
    C_raw_baseline_drift = running_percentile(roi(i).f, df_fwind, df_prctile)';
    roi(i).f0 =  C_raw_baseline_drift;     
    roi(i).df_f0 = (roi(i).f - roi(i).f0)./roi(i).f0; 
end


% generate figures 

% 1. full trace of individual dendrites (subplot)
taxis = (1/fps):(1/fps):96;
tl = length(taxis); 
figure('color','w','renderer','painters','position',[596    68   795   901]); 
for i=1:nrois 
    subplot(nrois,1,i);
    %plot(taxis, roi(i).f,'-'); hold on;     
    plot(taxis, roi(i).df_f0(1:tl),'-'); hold on;  
end
set(gca,'box','off'); 

% 2. full trace of individual dendrites (together)
taxis = (1/fps):(1/fps):96;
tl = length(taxis); 
figure('color','w','renderer','painters'); 
for i=1:nrois 
    %plot(taxis, roi(i).f,'-'); hold on;     
    plot(taxis, roi(i).df_f0(1:tl),'-'); hold on;  
end
set(gca,'box','off','xlim',[0 96]); legend(); 

% 3. Trial by trial traces of individual dendrites 
time_per_trial = 96/nrepeats; % sec 
taxis_ = (1/fps):(1/fps):time_per_trial;
tl_ = length(taxis_); 
figure('color','w','renderer','painters','position',[596    68   795   901]); 
for i=1:nrois 
    subplot(nrois,1,i);
    %plot(taxis, roi(i).f,'-'); hold on; 
    accum = []; 
    for k=1:nrepeats
        plot(taxis_, roi(i).df_f0((k-1)*tl_+1:k*tl_), '-','color',0.3.*[1 1 1]); hold on; 
        accum = [accum; roi(i).df_f0((k-1)*tl_+1:k*tl_)]; 
    end
    plot(taxis_, median(accum,1),'-k','linewidth',3); 
    set(gca,'box','off','ylim',[0 0.15]); 
    if i~=5
        set(gca,'box','off','ytick',[],'ycolor','w','xtick',[],'xcolor','w'); 
    end
end



% 4. Overlay selected ROIs on the image 
clx = [250 196 35]./255; 
figure('color','w','renderer','painters'); 
imshow(std_img_new,[25.1007 70.5116],'InitialMagnification','fit'); axis equal; hold on; 
for i=1:nrois 
    pf = fill(roi(i).position(:,1),roi(i).position(:,2),clx);
    pf.EdgeColor = 'none'; 
    pf.FaceColor = clx; 
    pf.FaceAlpha = 0.85; 
end








