function [outliers] = extract_delta(obj, indices)


    traces = obj.C_raw(indices,:);

    npks = 5; % assuming outlier traces have no more than 5 delta peaks
    traces_pksrmv = zeros(size(traces,1),size(traces,2)-npks);
    hill = zeros(size(traces,1),1);
    for i=1:size(traces,1)
        tempvec = traces(i,:);
        [~, locs] = findpeaks(tempvec,'NPeaks',npks,'SortStr','descend');

        if isempty(locs)
            hill(i,1) = 0;
            traces_pksrmv(i,:) = traces(i,size(traces,2)-npks);
        elseif ~isempty(locs) && length(locs)<npks
            hill(i,1) = mean(tempvec(locs));
            tempvec(locs) = [];
            tempvec(end-(npks - length(locs))+1:end) = [];
            traces_pksrmv(i,:) = tempvec;
        else
            hill(i,1) = mean(tempvec(locs));
            tempvec(locs) = [];
            traces_pksrmv(i,:) = tempvec;
        end
    end

    stdnormheight = hill./std(traces_pksrmv,0,2);

    idx = isoutlier(stdnormheight);
    isnan_idx = find(isnan(stdnormheight)); 
    
    if sum(isnan_idx)>0
        exclude_nan_stdnormheight = stdnormheight((idx==0)); 
        exclude_nan_stdnormheight = exclude_nan_stdnormheight(~isnan(exclude_nan_stdnormheight)); 

        if mean(stdnormheight((idx==1))) > mean(exclude_nan_stdnormheight)
            outliers = find(idx==1);
        elseif mean(stdnormheight((idx==1))) < mean(exclude_nan_stdnormheight)
            outliers = find(idx==0);
        else
            if median(stdnormheight((idx==1))) >= median(exclude_nan_stdnormheight)
                outliers = find(idx==1);
            elseif median(stdnormheight((idx==1))) < median(exclude_nan_stdnormheight)
                outliers = find(idx==0);
            end
        end

        outliers = unique([isnan_idx; outliers]);  
    else
        outliers = []; 
    end

    if ~isempty(outliers)
        include = setxor(1:size(traces,1), outliers);
        figure;
        subplot(2,1,1);
        splot(1:size(traces,2), 1:length(outliers), traces(outliers,:)');
        title('outliers');
        subplot(2,1,2);
        splot(1:size(traces,2), 1:length(include),  traces(include, :)');
        title('kept');
    end

end

