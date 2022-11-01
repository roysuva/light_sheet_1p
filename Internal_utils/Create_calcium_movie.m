%% Add file paths 

addpath(genpath('/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Imaging/Light_sheet/Analysis/Matlab/Internal_utils/'));
addpath(genpath('/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Imaging/Light_sheet/Analysis/Matlab/External_utils/'));


%% Set tiff filepath 
% impath = ['/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Imaging/Light_sheet/Analysis/2022-01-07-0/',...
%     '/RE_ventral/Data_000001/'];
clear;
impath = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Imaging/Light_sheet/Analysis/2021-08-17-0/LE_ventral/Data_100001/';

%% Read all tiff files 
imall= read_tiff_files(impath); 

%% Apply Kalman filtering to reduce noise 

imall_ = double(imall); 
imall_filt = Kalman_Stack_Filter(imall_,0.75);
% imall_filt = zeros(size(imall_)); 
% parfor k=1:size(imall_,3)
%     imall_filt(:,:,k) = medfilt2(imall_(:,:,k),[1 1],'symmetric');
% end


%% Select sub-area 

if ~exist('reng_rec','var')
    rng_rec = [prctile(imall_filt(:),2) prctile(imall_filt(:),99)]; 
end

d1 = size(imall_filt,1);
d2 = size(imall_filt,2);
d3 = size(imall_filt,3);

ui = input('Do you want to select sub-area?  y/n  >  ','s'); 

if strcmpi(ui,'y')
    figure;
    imshow(imall_filt(:,:,10),rng_rec,'InitialMagnification','fit');  hold on;  
    [x,y] = deal([]);
    for i=1:2
        [x(i),y(i)] = ginput(1);
        x(i) = round(x(i));
        y(i) = round(y(i));
        if x(i)<1; x(i)=1; end
        if y(i)<1; y(i)=1; end
        if x(i)>size(imall_filt,2); x(i)=size(imall_filt,2); end
        if y(i)>size(imall_filt,1); y(i)=size(imall_filt,1); end
        h = plot(x(i),y(i),'or','markersize',8,'linewidth',3); hold on;
        set(get(get(h,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    end
    AOI_x = [x(1) x(2) x(2) x(1) x(1)]; % Area of Interest
    AOI_y = [y(1) y(1) y(2) y(2) y(1)];
    plot(AOI_x,AOI_y,'-y','linewidth',3); legend('Area for analysis');
else
    AOI_x = [1 d2 d2 1 1];
    AOI_y = [1 1 d1 d1 1];
end

imall_filt = imall_filt(min(AOI_y):max(AOI_y),min(AOI_x):max(AOI_x),:);  
 
 
% figure
% %rng_rec = [prctile(imall_filt(:),5) prctile(imall_filt(:),99.5)];
% imshow(imall_filt(:,:,10),rng_rec,'InitialMagnification','fit'); 
% imcontrast

save('/home/circuit/Downloads/temp_movie_data.mat','imall_filt','rng_rec','-v7.3');

%% Save video 

savefile = '/home/circuit/Downloads/GCaMP_movie2.mp4'; 
if exist(savefile,'file').mp4
    delete(savefile);
end
v = VideoWriter(savefile,'MPEG-4'); 
v.FrameRate = 120; 
v.Quality = 100; 
open(v);  

hf = figure;  set(hf,'color','w','renderer','painters'); 
%rng_rec = [prctile(imall_filt(:),2) prctile(imall_filt(:),99)];
N = size(imall_filt,3); 
parfor_progress(N);
for k=1:N
    imshow(imall_filt(min(AOI_y):max(AOI_y),min(AOI_x):max(AOI_x),k),rng_rec);  
    frame = getframe(hf); 
    writeVideo(v,frame); 
    parfor_progress();
end
parfor_progress(0);
close(v); 
close(hf); 









