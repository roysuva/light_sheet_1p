function [im_uint16] = read_tiff_files(tiffpath)

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


% Read image files
clm = find(chunklim(:,1)<=image_rng(1),1,'last'):find(chunklim(:,2)>=image_rng(2),1,'first');
im_uint16_temp = cell(1,length(clm));
for fln=1:length(clm)
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
