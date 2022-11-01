function subimages = readMultipageTiff(tiffpath)


% Check if there are multiple multipage tiffs
tiffdirinfo = dir(fullfile(tiffpath, '*.tif'));  
num_mptfs = length(tiffdirinfo); 

% Process images 
initfr = 0; cnt =1 ; 
subimages = zeros([576 576 1],'uint16'); 
for fln=1%:num_mptfs
    t = Tiff(fullfile(tiffpath,tiffdirinfo(fln).name),'r');
    subimages(:,:,1+initfr) = t.read(); 
    if t.lastDirectory()
        return; 
    end
    t.nextDirectory();
    while true
        subimages(:,:,end+1) = t.read();
        if t.lastDirectory()
            break; 
        else
            t.nextDirectory(); 
        end
        cnt
        cnt = cnt+1;
    end
    lastfr = size(subimages,3); 
%     if lastfr>framerange(2)
%         break;
%     else
%         initfr = lastfr; 
%     end
end
% 
% % Extract specified range 
% subimages(:,:,framerange(2)+1:end) = []; 
% subimages(:,:,1:framerange(1)-1) = []; 

