function [bw_images, stim] = get_bw_images(stim)
    
    if ~isfield(stim.stimulus,'refresh_rate')
        stim.stimulus.refresh_rate = 60; 
    end
    init_delay = stim.stimulus.delay_frames/stim.stimulus.refresh_rate; % sec 
    
    trigtimes = stim.time_stamps.stim_time_stamps(2:end) - init_delay; 
    nframes = stim.stimulus.frames; 
    
    nimages = floor(nframes/stim.stimulus.refresh); 
    bw_images = get_movie(stim.movie_file, trigtimes(1:length(trigtimes)),nimages); 
    bw_images = permute(bw_images,[2 1 3 4]);
    
    fprintf(sprintf('Extracted %d images with refresh = %d\n',nimages,stim.stimulus.refresh)); 
end
 