function [AOI, FOV, median_im] = get_AOI(AOIpath)

    % This function extracts the size of the full image (FOV) and the 
    % subregion (AOI) used for analysis
    % 
    % Inputs: 
    %       AOIpath: filepath for matlab .fig of AOI selection (look up
    %       CNMF-E analysis script
    %
    % Outputs: 
    %       AOI: Area of interest 
    %       FOV: Full field of view 
    %       median_im: median image 
    
    figh = openfig(AOIpath,'invisible');
    hl = findobj(gca,'type','line');
    AOI.x = hl(1).XData;
    AOI.y = hl(1).YData;

    fco = findobj(get(figh,'Children'),'Type','Image');
    median_im = fco(1).CData; 
    sz = size(median_im); 

    FOV.x = [1 sz(2) sz(2) 1 1];
    FOV.y = [1 1 sz(1) sz(1) 1];

end







