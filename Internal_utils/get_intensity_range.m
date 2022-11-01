function [intensity_range] = get_intensity_range(im, varargin)
    

    % Specify list of optional parameters
    p = inputParser;
    addParameter(p,'range', [0.1 0.9], @isnumeric);
    p.parse(varargin{:});
    prms = p.Results;
    
    qtl1 = prms.range(1);
    qtl2 = prms.range(2); 
    
    intensity_range = [double((median(im(:))-quantile(im(:),qtl1))<=min(im(:)))*min(im(:)) + ...
        double( (median(im(:))-quantile(im(:),qtl1))>min(im(:)) )*(median(im(:))-quantile(im(:),qtl1)) ...
        double( (median(im(:))+quantile(im(:),qtl2))>=max(im(:)) )*max(im(:)) + ...
        double( (median(im(:))+quantile(im(:),qtl2))<max(im(:)) )*(median(im(:))+quantile(im(:),qtl2))];

end