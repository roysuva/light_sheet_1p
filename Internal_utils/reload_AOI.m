function [AOI_x, AOI_y] = reload_AOI(filepath)

    [fpath, fname, fext] = fileparts(filepath); 
    if ~strcmpi(fext, '.fig')
        error('File needs to be Matlab .fig');
    else 
        uiopen(filepath, true); 
        hl = findobj('Type','line'); 
        AOI_x = hl(1).XData;
        AOI_y = hl(1).YData;
        
    end
    
    return;
end