function [outliers] = extract_outliers(varargin)
    
    
    % Get input
    p = inputParser;  
    addParameter(p,'traces',[], @isnumeric);
    addParameter(p,'show_plot',true, @islogical);
    parse(p,varargin{:});
    params = p.Results; 
    
    
    % Reshape trace matrix 
    traces = params.traces'; % rows: variables, columns: values 


    npks = 5; % assuming outlier traces have no more than 5 delta peaks
    nvars = size(traces,1); 
    nvals = size(traces,2); 
    stdnormheight = zeros(nvars,1);
    for i=1:nvars
        tempvec = traces(i,:);
        [~, locs] = findpeaks(tempvec,'NPeaks',npks,'SortStr','descend'); 

        if ~isempty(locs)
            tempvecx = tempvec(setxor(1:nvals,locs)); 
            xxc = find(tempvec(locs) > 10*mad(tempvecx)); 
            locs_trunc = locs(xxc); 
            if ~isempty(locs_trunc) 
                hill = mean(tempvec(locs_trunc));  
                traces_pksrmv = tempvec(setxor(1:nvals,locs_trunc));  
                stdnormheight(i,1) = hill./std(traces_pksrmv);  
            else 
                stdnormheight(i,1) = 0;
            end
        else 
            stdnormheight(i,1) = 0; 
        end
    end

    lcc = find(stdnormheight~=0); 
    outliers = lcc((isoutlier(stdnormheight(lcc)))); 
    

%     if ~isempty(outliers) && params.show_plot
%         include = setxor(1:size(traces,1), outliers);
%         figure;
%         subplot(2,1,1);
%         splot(1:size(traces,2), 1:length(outliers), traces(outliers,:)');
%         title('outliers');
%         subplot(2,1,2);
%         splot(1:size(traces,2), 1:length(include),  traces(include, :)');
%         title('kept');
%     end

end

