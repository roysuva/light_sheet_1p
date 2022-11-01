function [resp , StimComb, stimout] = get_resp_properties_cat(stim, neuronObj, frame2stim, frames_cum, varargin)


    
    global fps;

    p = inputParser;
    addParameter(p,'detrend',true, @islogical);
    addParameter(p,'method','running_percentile', @ischar);
    parse(p,varargin{:});
    params = p.Results;

    StimComb = [];  
    
    fprintf(sprintf('%s Analyzing responses to %s stimulus %s\n',repelem('.',30),stim.type,repelem('.',30)));
    
    switch stim.type
        case 'MB'
            
            for q=1
                
                s = fieldnames(stim.combinations);
                StimComb = zeros(length(stim.combinations),length(fieldnames(stim.params)));
                for i=1:length(s)
                    if strcmpi(s{i}, 'back_rgb') || strcmpi(s{i}, 'rgb')
                        for j=1:length(stim.combinations)
                            StimComb(j,i) = mean(stim.combinations(j).(s{i}));
                        end
                    else
                        StimComb(:,i) = [stim.combinations(:).(s{i})]';
                    end
                end
                
                
                % Set frame range for analysis
                frm_index = frame2stim.frame_rng(1,1):frame2stim.frame_rng(end,2);
                frm_rng_corected = frm_index - frames_cum ;
                
                ncells = numel(neuronObj.ids);
                nrepeats = stim.repetitions;
                resp.t_axis = (1:length(frm_rng_corected))./fps;
                resp.t_axis_trial = (1:diff(frame2stim.frame_rng(1,:))+1)./fps;
                
                
                % Get stimulus parameters
                S = [find(strcmp(fieldnames(stim.combinations),'delta')) find(strcmp(fieldnames(stim.combinations),'DELTA'))];
                unique_speed = unique(StimComb(:,S));
                W = [find(strcmp(fieldnames(stim.combinations),'bar_width')) find(strcmp(fieldnames(stim.combinations),'BAR_WIDTH'))];
                unique_width = unique(StimComb(:,W));
                C = [find(strcmp(fieldnames(stim.combinations),'rgb')) find(strcmp(fieldnames(stim.combinations),'RGB'))];
                unique_contrast = unique(StimComb(:,C));
                D = [find(strcmp(fieldnames(stim.combinations),'direction')) find(strcmp(fieldnames(stim.combinations),'DIRECTION'))];
                unique_dirs = unique(StimComb(:,D));
                
                %                     % correction for missing frames at the end (if)
                %                     if ~isfield(stim.time_stamps,'external_trig_times_corrected')
                %                         external_trig_times = stim.time_stamps.external_trig_times;
                %                     else
                %                         external_trig_times = stim.time_stamps.external_trig_times_corrected;
                %                     end
                %                     if length(external_trig_times)<frame2stim.frame_rng(end,2)
                %                         numfr_missing = abs(length(external_trig_times)-frame2stim.frame_rng(end,2));
                %                         lntrial = mean(diff(frame2stim.frame_rng,[],2));
                %                         fprintf(sprintf('%1.2f fraction of the last trial missing.\n', numfr_missing/lntrial));
                %                         fprintf(sprintf('Reducing repetitions from %d >>>> %d.\n', stim.repetitions, stim.repetitions-1));
                %                         nrepeats = nrepeats -1;
                %                         stim.repetitions_corrected = nrepeats;
                %                     end
                
                
                % correction for missing frames at the end (if)
                if ~isfield(stim.time_stamps,'external_trig_times_corrected')
                    external_trig_times = stim.time_stamps.external_trig_times;
                else
                    external_trig_times = stim.time_stamps.external_trig_times_corrected;
                end
                frames_from_cam = fix(external_trig_times(end)-external_trig_times(1))*fps;
                frames_from_stim = frm_index(end)-frm_index(1)+1;
                if frames_from_stim>frames_from_cam
                    numfr_missing = frames_from_stim-frames_from_cam;
                    lnrep = fix(stim.trial_dur)*length(stim.combinations)*fps;
                    fprintf(sprintf('%1.2f fraction of the last repeat missing.\n', numfr_missing/lnrep));
                    fprintf(sprintf('Reducing repetitions from %d >>>>>> %d.\n', stim.repetitions, stim.repetitions-ceil(numfr_missing/lnrep)));
                    nrepeats = nrepeats -ceil(numfr_missing/lnrep);
                    stim.repetitions_corrected = nrepeats;
                    true_rng = frm_rng_corected(1)+(0:frames_from_cam-1);
                else
                    fprintf(sprintf('All repeats available.\n'));
                    stim.repetitions_corrected = nrepeats;
                    true_rng = frm_rng_corected(1)+(0:frames_from_stim-1);
                end
                
                % set raster for cells
                resp.detrended_trace = zeros(ncells,length(true_rng));
                periods = frame2stim.frame_rng - frames_cum;
                for nc=1:ncells
                    
                    trace = neuronObj.C_raw(nc,true_rng);
                    if numel(resp.t_axis)>numel(trace)
                        resp.t_axis = resp.t_axis(1:length(trace));
                    end
                    
                    % if detrend
                    if params.detrend
                        switch params.method
                            case 'running_percentile'
                                df_prctile = 20; % use low fluorescence values (excludes stimulus induced transients)
                                df_twind = 30; % seconds
                                df_fwind = df_twind*fps;
                                C_raw_baseline_drift = running_percentile(trace, df_fwind, df_prctile)';
                                
                                trace = trace - C_raw_baseline_drift;
%                                 if min(trace)<0
%                                     trace = trace - min(trace); % fluorescence values should be positive
%                                 end
                            case 'spline fit'
                        end
                    end
                    resp.detrended_trace(nc,:) = trace;
                    
                    if sum(trace)~=0
                        % filter trace
                        resp.detrended_filt_trace(nc,:) = get_filtered_trace(trace, fps);
                    else
                        resp.detrended_filt_trace(nc,:) = zeros(1,length(trace));
                    end
                    
                    
                    % extract rasters for stimulus combinations
                    trace2use = resp.detrended_filt_trace(nc,:);
                    traceunfilt = trace;
                    for i=1:length(unique_speed)
                        for j=1:length(unique_width)
                            for k=1:length(unique_contrast)
                                for m=1:length(unique_dirs)
                                    rc = find(StimComb(:,D)==unique_dirs(m));
                                    runid = find(stim.trial_list==rc);
                                    for ntr=1:nrepeats
                                        resp.unique_speed(i).unique_width(j).unique_contrast(k).unique_dirs(m).cell(nc).raster_traces(ntr,:) = ...
                                            trace2use(periods(runid(ntr),1):periods(runid(ntr),2));
                                        
                                        resp.unique_speed(i).unique_width(j).unique_contrast(k).unique_dirs(m).cell(nc).raster_unfilt_traces(ntr,:) = ...
                                            traceunfilt(periods(runid(ntr),1):periods(runid(ntr),2));
                                    end
                                end
                            end
                        end
                    end
                    fprintf(sprintf('%d/%d \n',nc,numel(neuronObj.ids)));
                end
                
                mylist = {};
                for i=1:length(unique_speed)
                    mylist{i} = sprintf('MB speed: %d',unique_speed(i));
                end
                for i=length(mylist)+1:length(unique_width)+length(mylist)
                    mylist{i} = sprintf('MB width: %d',unique_width(i-length(mylist)));
                end
                for i=length(mylist)+1:length(unique_contrast)+length(mylist)
                    mylist{i} = sprintf('MB contrast: %d',unique_contrast(i-length(mylist)));
                end
                
                [indx, ~] = listdlg('PromptString','Select a combination of bar speed, bar width and contrast for DSI calculation',...
                    'SelectionMode','multiple','ListString',mylist);
                [chosenval,paramids] = deal(zeros(1,length(indx)));
                for i=1:3
                    chosenval(i) = str2double(regexprep(mylist{indx(i)}(regexpi(mylist{indx(i)},':')+1:end),' ',''));
                end
                paramids(1) = find(unique_speed==chosenval(1));
                paramids(2) = find(unique_width==chosenval(2));
                paramids(3) = find(unique_contrast==chosenval(3));
                
                
                % calculate dsi
                mylist = cell(1,nrepeats);
                for ntr=1:nrepeats
                    mylist{ntr} = sprintf('%d',ntr);
                end
                mylist{nrepeats+1} = 'All repeats';
                [indx, ~] = listdlg('PromptString','Select which repeats to use for calculating DSI',...
                    'SelectionMode','multiple','ListString',mylist);
                if indx==nrepeats+1
                    repeats2use = 1:nrepeats;
                else
                    repeats2use = indx;
                end
                
                
                str2use = resp.unique_speed(paramids(1)).unique_width(paramids(2)).unique_contrast(paramids(3));
                ncells = length(str2use.unique_dirs(1).cell);
                theta = unique_dirs.*pi./180;
                rho = zeros(ncells,length(unique_dirs));
                [angle,mag] = deal(zeros(1,ncells));
                dF_F0 = cell(ncells,1);
                for nc=1:ncells
                    traparea = zeros(1,length(unique_dirs));  % area under curve (keep DC response: heterogeneous across directions)
                    dF_F0_temp = zeros(length(repeats2use),length(unique_dirs));
                    for d=1:length(unique_dirs)
                        
                        % remove baseline for each filtered trace and calculate area under curve
                        for nr=1:length(repeats2use)
                            ntr = repeats2use(nr);
                            current_trace = str2use.unique_dirs(d).cell(nc).raster_unfilt_traces(ntr,:);
                            traparea(d) = traparea(d) + sum(trapz(resp.t_axis_trial, current_trace));
                            
                            current_trace_ = current_trace - min(current_trace);
                            dF_F0_temp(nr,d) = abs(max(current_trace_)-prctile(current_trace_,10))/prctile(current_trace_,10);
                        end
                        traparea(d) = traparea(d)./length(repeats2use);
                    end
                    rho(nc,:) = traparea./sum(traparea);
                    
                    [X,Y] = pol2cart(theta', rho(nc,:));
                    u = sum(X);
                    v = sum(Y);
                    [angle(nc), mag(nc)] = cart2pol(u,v);
                    dF_F0{nc} = dF_F0_temp;
                end
                
                resp.repeats_used_for_dsi_estimate = repeats2use;
                resp.ds_estimates.rho = rho;
                resp.ds_estimates.theta = theta;
                resp.ds_estimates.mag = mag;
                resp.ds_estimates.angle = angle;
                resp.ds_estimates.dF_F0 = dF_F0;
                
                
                % select ds cells
                hf = figure(); clf;
                if size(mag,1)==1 || size(mag,2)==1
                    plot(1:ncells, mag,'ok','markerfacecolor','k'); hold on;
                    xll = get(gca,'xlim');
                    xlabel('Cell #'); ylabel('DSI'); set(gca,'xlim',[xll(1)-10 xll(2)+10]);
                    dim_ = 1;
                else
                    plot(log10(mag(:,1)), log10(mag(:,2)), 'ok','markerfacecolor','k'); hold on;
                    xlabel('log10(parameter 1)'); ylabel('log10(parameter 2)');
                    xl = min([log10(mag(:,1)); log10(mag(:,2))]);
                    set(gca,'xlim',[xl*1.2 0],'ylim',[xl*1.2 0]);
                    axis square;
                    dim_ = 2;
                end
                [x,y] = deal([]);
                but = 1;
                while ~isempty(but)
                    [x_,y_,but] = ginput(1);
                    x = [x x_];
                    y = [y y_];
                    plot(x,y,'.-r','markersize',8,'linewidth',1);
                end
                if dim_==1
                    IN = inpolygon(1:ncells, mag,x,y);
                    I = find(IN==1);
                    X = 1:ncells;
                    plot(X(I), mag(I),'ko','markerfacecolor','r'); set(gca,'xlim',[-1 inf]);
                elseif dim_==2
                    IN = inpolygon(log10(mag(:,1)), log10(mag(:,2)),x,y);
                    I = find(IN==1);
                    plot(log10(mag(I,1)), log10(mag(I,2)),'ko','markerfacecolor','r');
                end
                resp.ds_ids = I;
                resp.nonds_ids = setxor(1:ncells,I);
            end
            
        case 'MG'
            
            for q=1
                s = fieldnames(stim.combinations);
                StimComb = zeros(length(stim.combinations),length(fieldnames(stim.params)));
                for i=1:length(s)
                    if strcmpi(s{i}, 'back_rgb') || strcmpi(s{i}, 'rgb')
                        for j=1:length(stim.combinations)
                            StimComb(j,i) = mean(stim.combinations(j).(s{i}));
                        end
                    else
                        StimComb(:,i) = [stim.combinations(:).(s{i})]';
                    end
                end
                
                % Set frame range for analysis
                frm_index = frame2stim.frame_rng(1,1):frame2stim.frame_rng(end,2);
                frm_rng_corected = frm_index - frames_cum ;
                
                ncells = numel(neuronObj.ids);
                nrepeats = stim.repetitions;
                resp.t_axis = (1:length(frm_index))./fps;
                resp.t_axis_trial = (1:diff(frame2stim.frame_rng(1,:))+1)./fps;
                
                T = [find(strcmp(fieldnames(stim.combinations),'temporal_period')) find(strcmp(fieldnames(stim.combinations),'TEMPORAL_PERIOD'))];
                unique_tp = unique(StimComb(:,T));
                S = [find(strcmp(fieldnames(stim.combinations),'spatial_period')) find(strcmp(fieldnames(stim.combinations),'SPATIAL_PERIOD'))];
                unique_sp = unique(StimComb(:,S));
                C = [find(strcmp(fieldnames(stim.combinations),'rgb')) find(strcmp(fieldnames(stim.combinations),'RGB'))];
                unique_contrast = unique(StimComb(:,C));
                D = [find(strcmp(fieldnames(stim.combinations),'direction')) find(strcmp(fieldnames(stim.combinations),'DIRECTION'))];
                unique_dirs = unique(StimComb(:,D));
                
                % correction for missing frames at the end (if)
                if ~isfield(stim.time_stamps,'external_trig_times_corrected')
                    external_trig_times = stim.time_stamps.external_trig_times;
                else
                    external_trig_times = stim.time_stamps.external_trig_times_corrected;
                end
                frames_from_cam = fix(external_trig_times(end)-external_trig_times(1))*fps;
                frames_from_stim = frm_index(end)-frm_index(1)+1;
                if frames_from_stim>frames_from_cam
                    numfr_missing = frames_from_stim-frames_from_cam;
                    lnrep = fix(stim.trial_dur)*length(stim.combinations)*fps;
                    fprintf(sprintf('%1.2f fraction of the last repeat missing.\n', numfr_missing/lnrep));
                    fprintf(sprintf('Reducing repetitions from %d >>>>>> %d.\n', stim.repetitions, stim.repetitions-ceil(numfr_missing/lnrep)));
                    nrepeats = nrepeats -ceil(numfr_missing/lnrep);
                    stim.repetitions_corrected = nrepeats;
                    true_rng = frm_rng_corected(1)+(0:frames_from_cam-1);
                else
                    fprintf(sprintf('All repeats available'));
                    stim.repetitions_corrected = nrepeats;
                    true_rng = frm_rng_corected(1)+(0:frames_from_stim-1);
                end
                
                
                % set raster for cells
                resp.detrended_trace = zeros(ncells,length(true_rng));
                periods = frame2stim.frame_rng - frames_cum;
                for nc=1:ncells
                    
                    trace = neuronObj.C_raw(nc,true_rng);
                    if numel(resp.t_axis)>numel(trace)
                        resp.t_axis = resp.t_axis(1:length(trace));
                    end
                    
                    % if detrend
                    if params.detrend
                        switch params.method
                            case 'running_percentile'
                                df_prctile = 20; % use low fluorescence values (excludes stimulus induced transients)
                                df_twind = 30; % seconds
                                df_fwind = df_twind*fps;
                                C_raw_baseline_drift = running_percentile(trace, df_fwind, df_prctile)';
                                
                                trace = trace - C_raw_baseline_drift;
                                %                                     if min(trace)<0
                                %                                         trace = trace - min(trace); % fluorescence values should be positive
                                %                                     end
                            case 'spline fit'
                        end
                    end
                    resp.detrended_trace(nc,:) = trace;
                    
                    if sum(trace)~=0
                        % filter trace
                        resp.detrended_filt_trace(nc,:) = get_filtered_trace(trace, fps);
                    else
                        resp.detrended_filt_trace(nc,:) = zeros(1,length(trace));
                    end
                    
                    
                    % extract rasters for stimulus combinations
                    trace2use = resp.detrended_filt_trace(nc,:);
                    traceunfilt = trace;
                    for i=1:length(unique_tp)
                        for j=1:length(unique_sp)
                            for k=1:length(unique_contrast)
                                for m=1:length(unique_dirs)
                                    rc = find(StimComb(:,D)==unique_dirs(m));
                                    runid = find(stim.trial_list==rc);
                                    for ntr=1:nrepeats
                                        resp.unique_tp(i).unique_sp(j).unique_contrast(k).unique_dirs(m).cell(nc).raster_traces(ntr,:) = ...
                                            trace2use(periods(runid(ntr),1):periods(runid(ntr),2));
                                        
                                        resp.unique_tp(i).unique_sp(j).unique_contrast(k).unique_dirs(m).cell(nc).raster_unfilt_traces(ntr,:) = ...
                                            traceunfilt(periods(runid(ntr),1):periods(runid(ntr),2));
                                    end
                                end
                            end
                        end
                    end
                    fprintf(sprintf('%d/%d \n',nc,numel(neuronObj.ids)));
                end
                
                mylist = {};
                for i=1:length(unique_tp)
                    mylist{i} = sprintf('MG temporal period: %d',unique_tp(i));
                end
                for i=length(mylist)+1:length(unique_sp)+length(mylist)
                    mylist{i} = sprintf('MG spatial period: %d',unique_sp(i-length(mylist)));
                end
                for i=length(mylist)+1:length(unique_contrast)+length(mylist)
                    mylist{i} = sprintf('MB contrast: %d',unique_contrast(i-length(mylist)));
                end
                
                [indx, ~] = listdlg('PromptString','Select a combination of bar speed, bar width and contrast for DSI calculation',...
                    'SelectionMode','multiple','ListString',mylist);
                [chosenval,paramids] = deal(zeros(1,length(indx)));
                for i=1:3
                    chosenval(i) = str2double(regexprep(mylist{indx(i)}(regexpi(mylist{indx(i)},':')+1:end),' ',''));
                end
                paramids(1) = find(unique_tp==chosenval(1));
                paramids(2) = find(unique_sp==chosenval(2));
                paramids(3) = find(unique_contrast==chosenval(3));
                
                
                % calculate dsi
                mylist = cell(1,nrepeats);
                for ntr=1:nrepeats
                    mylist{ntr} = sprintf('%d',ntr);
                end
                mylist{nrepeats+1} = 'All repeats';
                [indx, ~] = listdlg('PromptString','Select which repeats to use for calculating DSI',...
                    'SelectionMode','multiple','ListString',mylist);
                if indx==nrepeats+1
                    repeats2use = 1:nrepeats;
                else
                    repeats2use = indx;
                end
                
                
                
                str2use = resp.unique_tp(paramids(1)).unique_sp(paramids(2)).unique_contrast(paramids(3));
                ncells = length(str2use.unique_dirs(1).cell);
                theta = unique_dirs.*pi./180;
                rho = zeros(ncells,length(unique_dirs));
                [angle_ds,mag_ds,angle_os,mag_os,ds_index,os_index] = deal(zeros(1,ncells));
                dF_F0 = cell(ncells,1);
                for nc=1:ncells
                    traparea = zeros(1,length(unique_dirs));  % area under curve (keep DC response: heterogeneous across directions)
                    dF_F0_temp = zeros(length(repeats2use),length(unique_dirs));
                    for d=1:length(unique_dirs)
                        
                        % remove baseline for each filtered trace and calculate area under curve
                        for nr=1:length(repeats2use)
                            ntr = repeats2use(nr);
                            current_trace = str2use.unique_dirs(d).cell(nc).raster_unfilt_traces(ntr,:);
                            traparea(d) = traparea(d) + sum(trapz(resp.t_axis_trial, current_trace));
                            
                            current_trace_ = current_trace - min(current_trace);
                            dF_F0_temp(nr,d) = abs(max(current_trace_)-prctile(current_trace_,10))/prctile(current_trace_,10);
                        end
                        traparea(d) = traparea(d)./length(repeats2use);
                    end
                    rho(nc,:) = traparea;
                    
                    [X,Y] = pol2cart(theta', rho(nc,:));
                    [XX,YY] =  pol2cart(2.*theta', rho(nc,:));
                    u = sum(X);
                    v = sum(Y);
                    uu = sum(XX);
                    vv = sum(YY);
                    [angle_ds(nc), mag_ds(nc)] = cart2pol(u,v);
                    [angle_os(nc), mag_os(nc)] = cart2pol(uu,vv);
                    dF_F0{nc} = dF_F0_temp;
                    
                    % Calculate orientation selective index
                    po = wrapTo360(angle_os(nc).*180./pi); % preferred orientation
                    max_resp_angs = sort(unique_dirs([find(unique_dirs<po,1,'last') find(unique_dirs>po,1,'first')]));
                    if length(max_resp_angs)==1
                        max_resp_angs = [0 max_resp_angs];
                    end
                    [~,max_resp_inds,~] = intersect(unique_dirs, max_resp_angs,'stable');
                    if po>unique_dirs(end) && po<360
                        po_neighbor_angs = [360-po po-unique_dirs(end)];
                    else
                        po_neighbor_angs = abs(max_resp_angs - po);
                    end
                    
                    no = wrapTo360(po+180);
                    min_resp_angs = sort(wrapTo360(max_resp_angs+180));
                    [~,min_resp_inds,~] = intersect(unique_dirs, min_resp_angs,'stable');
                    if no>unique_dirs(end) && no<360
                        no_neighbor_angs = [360-no no-unique_dirs(end)];
                    else
                        no_neighbor_angs = abs(min_resp_angs - no);
                    end
                    
                    ortho1 = wrapTo360(po+90);
                    ortho1_resp_angs = sort(unique_dirs([find(unique_dirs<ortho1,1,'last') find(unique_dirs>ortho1,1,'first')]));
                    if length(ortho1_resp_angs)==1
                        ortho1_resp_angs = [0 ortho1_resp_angs];
                    end
                    [~,ortho1_resp_inds,~] = intersect(unique_dirs, ortho1_resp_angs,'stable');
                    if ortho1>unique_dirs(end) && ortho1<360
                        ortho1_neighbor_angs = [360-ortho1 ortho1-unique_dirs(end)];
                    else
                        ortho1_neighbor_angs = abs(ortho1_resp_angs - ortho1);
                    end
                    
                    ortho2 = wrapTo360(po-90);
                    ortho2_resp_angs = sort(unique_dirs([find(unique_dirs<ortho2,1,'last') find(unique_dirs>ortho2,1,'first')]));
                    if length(ortho2_resp_angs)==1
                        ortho2_resp_angs = [0 ortho2_resp_angs];
                    end
                    [~,ortho2_resp_inds,~] = intersect(unique_dirs, ortho2_resp_angs,'stable');
                    if ortho2>unique_dirs(end) && ortho2<360
                        ortho2_neighbor_angs = [360-ortho2 ortho2-unique_dirs(end)];
                    else
                        ortho2_neighbor_angs = abs(ortho2_resp_angs - ortho2);
                    end
                    
                    if size(po_neighbor_angs,1)>size(po_neighbor_angs,2); po_neighbor_angs = po_neighbor_angs'; end
                    if size(no_neighbor_angs,1)>size(no_neighbor_angs,2); no_neighbor_angs = no_neighbor_angs'; end
                    if size(ortho1_neighbor_angs,1)>size(ortho1_neighbor_angs,2); ortho1_neighbor_angs = ortho1_neighbor_angs'; end
                    if size(ortho2_neighbor_angs,1)>size(ortho2_neighbor_angs,2); ortho2_neighbor_angs = ortho2_neighbor_angs'; end
                    
                    po_resp = mean(rho(nc,max_resp_inds).*cosd(po_neighbor_angs));
                    no_resp = mean(rho(nc,min_resp_inds).*cosd(no_neighbor_angs));
                    ortho1_resp = mean(rho(nc,ortho1_resp_inds).*cosd(ortho1_neighbor_angs));
                    ortho2_resp = mean(rho(nc,ortho2_resp_inds).*cosd(ortho2_neighbor_angs));
                    
                    os_index(nc) = abs(1 - mean([ortho1_resp ortho2_resp])/po_resp);
                    
                    nume = hypot(sum(rho(nc,:).*cosd(unique_dirs')), sum(rho(nc,:).*sind(unique_dirs')));
                    deno = sum(rho(nc,:));
                    %os_index(nc) = (nume/deno);
                    
                    % Calculate direction selective index
                    pd = wrapTo360(angle_ds(nc).*180./pi); % preferred direction
                    max_resp_angs = sort(unique_dirs([find(unique_dirs<pd,1,'last') find(unique_dirs>pd,1,'first')]));
                    if length(max_resp_angs)==1
                        max_resp_angs = [0 max_resp_angs];
                    end
                    [~,max_resp_inds,~] = intersect(unique_dirs, max_resp_angs,'stable');
                    nd = wrapTo360(pd+180);
                    min_resp_angs = sort(wrapTo360(max_resp_angs+180));
                    [~,min_resp_inds,~] = intersect(unique_dirs, min_resp_angs,'stable');
                    if pd>unique_dirs(end) && pd<360
                        pd_neighbor_angs = [360-pd pd-unique_dirs(end)];
                    else
                        pd_neighbor_angs = abs(max_resp_angs - pd);
                    end
                    if nd>unique_dirs(end) && nd<360
                        nd_neighbor_angs = [360-nd nd-unique_dirs(end)];
                    else
                        nd_neighbor_angs = abs(min_resp_angs - nd);
                    end
                    
                    pd_resp = mean(rho(nc,max_resp_inds).*cosd(pd_neighbor_angs));
                    nd_resp = mean(rho(nc,min_resp_inds).*cosd(nd_neighbor_angs));
                    ds_index(nc) = abs(pd_resp-nd_resp)/(pd_resp);
                end
                
                resp.repeats_used = repeats2use;
                
                resp.ds_estimates.rho = rho;
                resp.ds_estimates.theta = theta;
                resp.ds_estimates.mag = mag_ds;
                resp.ds_estimates.angle = angle_ds;
                resp.ds_estimates.ds_index = ds_index;
                
                resp.os_estimates.rho = rho;
                resp.os_estimates.theta = theta;
                resp.os_estimates.mag = mag_os;
                resp.os_estimates.angle = angle_os;
                resp.os_estimates.os_index = os_index;
                
                
                
                % select ds cells
                ds_os = input('Do you want to select DS cells (1) or OS cells (2)? ','s');
                
                hf = figure(); clf;
                
                if strcmpi(ds_os, '1')
                    mag = ds_index;
                    resp.which_index = 'ds';
                    ylabel('DSI');
                    xlabel('Cell #');
                    if size(mag,1)==1 || size(mag,2)==1
                        plot(1:ncells, mag,'ok','markerfacecolor','k'); hold on;
                        xll = get(gca,'xlim');
                        set(gca,'xlim',[xll(1)-10 xll(2)+10]);
                        dim_ = 1;
                    else
                        plot(log10(mag(:,1)), log10(mag(:,2)), 'ok','markerfacecolor','k'); hold on;
                        xlabel('log10(parameter 1)'); ylabel('log10(parameter 2)');
                        xl = min([log10(mag(:,1)); log10(mag(:,2))]);
                        set(gca,'xlim',[xl*1.2 0],'ylim',[xl*1.2 0]);
                        axis square;
                        dim_ = 2;
                    end
                    [x,y] = deal([]);
                    but = 1;
                    while ~isempty(but)
                        [x_,y_,but] = ginput(1);
                        x = [x x_];
                        y = [y y_];
                        plot(x,y,'.-r','markersize',8,'linewidth',1);
                    end
                    if dim_==1
                        IN = inpolygon(1:ncells, mag,x,y);
                        I = find(IN==1);
                        X = 1:ncells;
                        plot(X(I), mag(I),'ko','markerfacecolor','r'); set(gca,'xlim',[-1 inf]);
                    elseif dim_==2
                        IN = inpolygon(log10(mag(:,1)), log10(mag(:,2)),x,y);
                        I = find(IN==1);
                        plot(log10(mag(I,1)), log10(mag(I,2)),'ko','markerfacecolor','r');
                    end
                    
                elseif strcmpi(ds_os, '2')
                    mag = os_index;
                    resp.which_index = 'os';
                    ylabel('OSI');
                    xlabel('Cell #');
                    if size(mag,1)==1 || size(mag,2)==1
                        plot(1:ncells, mag,'ok','markerfacecolor','k'); hold on;
                        xll = get(gca,'xlim');
                        set(gca,'xlim',[xll(1)-10 xll(2)+10]);
                        dim_ = 1;
                    else
                        plot(log10(mag(:,1)), log10(mag(:,2)), 'ok','markerfacecolor','k'); hold on;
                        xlabel('log10(parameter 1)'); ylabel('log10(parameter 2)');
                        xl = min([log10(mag(:,1)); log10(mag(:,2))]);
                        set(gca,'xlim',[xl*1.2 0],'ylim',[xl*1.2 0]);
                        axis square;
                        dim_ = 2;
                    end
                    [x,y] = deal([]);
                    but = 1;
                    while ~isempty(but)
                        [x_,y_,but] = ginput(1);
                        x = [x x_];
                        y = [y y_];
                        plot(x,y,'.-r','markersize',8,'linewidth',1);
                    end
                    if dim_==1
                        IN = inpolygon(1:ncells, mag,x,y);
                        I = find(IN==1);
                        X = 1:ncells;
                        plot(X(I), mag(I),'ko','markerfacecolor','r'); set(gca,'xlim',[-1 inf]);
                    elseif dim_==2
                        IN = inpolygon(log10(mag(:,1)), log10(mag(:,2)),x,y);
                        I = find(IN==1);
                        plot(log10(mag(I,1)), log10(mag(I,2)),'ko','markerfacecolor','r');
                    end
                    %                 mag1 = os_index;
                    %                 mag2 = ds_index;
                    %                 resp.which_index = 'os';
                    %                 plot(mag1, mag2, 'ok','markerfacecolor','k'); hold on;
                    %                 plot([0 1],[0 0],'-k',[0 0],[0 1],'-k',[0.5 0.5],[0 1],'-k',[0 1],[0.5 0.5],'-k');
                    %                 xll = get(gca,'xlim');
                    %                 yll = get(gca,'ylim');
                    %                 set(gca,'xlim',[xll(1)-0.1 xll(2)+0.1],'ylim',[yll(1)-0.1 yll(2)+0.1]);
                    %                 xlabel('OSI'); ylabel('DSI');
                    %                 [x,y] = deal([]);
                    %                 but = 1;
                    %                 while ~isempty(but)
                    %                     [x_,y_,but] = ginput(1);
                    %                     x = [x x_];
                    %                     y = [y y_];
                    %                     plot(x,y,'.-r','markersize',8,'linewidth',1);
                    %                 end
                    %                 IN = inpolygon(mag1 , mag2, x,y);
                    %                 I = find(IN==1);
                    %                 plot( mag1(I), mag2(I), 'ko','markerfacecolor','r');
                end
                
                
                if strcmpi(ds_os, '1') % ds cells
                    resp.ds_ids = I;
                    resp.nonds_ids = setxor(1:ncells,I);
                elseif strcmpi(ds_os, '2') % os cells
                    resp.os_ids = I;
                    resp.nonos_ids = setxor(1:ncells,I);
                end
                
                
                
            end
            
        case 'CH'
            
            for q=1
                
                StimComb = [];
                
                % Set frame range for analysis
                frm_index = frame2stim.frame_rng_cam(1,1):frame2stim.frame_rng_cam(end,2);
                
                s = fieldnames(stim);
                ncells = numel(neuronObj.ids);
                nrepeats = stim.stimulus.num_trials;
                resp.t_axis = (1:length(frm_index))./fps;
                trial_frame_rng = frame2stim.frame_rng_cam(1,1):frame2stim.frame_rng_cam(find(frame2stim.trial==1,1,'last'),2);
                resp.t_axis_trial = (1:length(trial_frame_rng))/fps;
                
                % Determine the number of frames in a trial
                transits = [0; find(diff(frame2stim.trial)~=0); length(frame2stim.trial)];
                sortfrmnum = zeros(nrepeats,1);
                for ntr=1:nrepeats
                    sortfrmnum(ntr,1) = abs(diff([frame2stim.frame_rng_cam(transits(ntr)+1,1) ...
                        frame2stim.frame_rng_cam(transits(ntr+1),2)]));
                end
                sortfrmnum = min(sortfrmnum);
                resp.t_axis_trial = resp.t_axis_trial(1:sortfrmnum);
                
                % correct for initial offset
                frame_rng_corrected = frame2stim.frame_rng_cam - frame2stim.frame_rng_cam(1,1)+1;
                
                for nc=1:ncells
                    
                    trace = neuronObj.C_raw(nc,frm_index);
                    if numel(resp.t_axis)>numel(trace)
                        resp.t_axis = resp.t_axis(1:length(trace));
                    end
                    
                    % if detrend
                    if params.detrend
                        switch params.method
                            case 'running_percentile'
                                df_prctile = 20; % use low fluorescence values (excludes stimulus induced transients)
                                df_twind = 10; % seconds
                                df_fwind = df_twind*fps;
                                C_raw_baseline_drift = running_percentile(trace, df_fwind, df_prctile)';
                                
                                trace = trace - C_raw_baseline_drift;
                                %                                     if min(trace)<0
                                %                                         trace = trace - min(trace); % fluorescence values should be positive
                                %                                     end
                            case 'spline fit'
                        end
                    end
                    resp.detrended_trace(nc,:) = trace;
                    
                    if sum(trace)~=0
                        
                        % filter trace
                        %resp.detrended_filt_trace(nc,:) = get_filtered_trace(trace, fps);
                        resp.detrended_filt_trace(nc,:) = get_filtered_trace(trace, fps, ...
                            'filter_type','ButterWorth','order',5, 'cutoff_dB',6);
                        
                        % extract rasters for stimulus combinations
                        trace2use = resp.detrended_trace(nc,:);
                        trace2use_filt = resp.detrended_filt_trace(nc,:);
                        
                        
                        for ntr=1:nrepeats
                            tempvec = trace2use(frame_rng_corrected(transits(ntr)+1,1):...
                                frame_rng_corrected(transits(ntr+1),2));
                            resp.cell(nc).raster(ntr,:) = tempvec(1,1:sortfrmnum);
                            tempvec = trace2use_filt(frame_rng_corrected(transits(ntr)+1,1):...
                                frame_rng_corrected(transits(ntr+1),2));
                            resp.cell(nc).raster_filtered(ntr,:) = tempvec(1,1:sortfrmnum);
                        end
                        
                    else
                        resp.detrended_filt_trace(nc,:) = zeros(1,length(frm_index));
                        resp.cell(nc).raster(1:nrepeats,1:sortfrmnum) = zeros(nrepeats,sortfrmnum);
                        resp.cell(nc).raster_filtered(1:nrepeats,1:sortfrmnum) = zeros(nrepeats,sortfrmnum);
                    end
                    
                    fprintf(sprintf('%d/%d \n',nc,numel(neuronObj.ids)));
                end
                
                % get stimulus trace in camera reference frame 
                taxis_stim = linspace(1,frame_rng_corrected(transits(2),2),...
                    frame2stim.frame_rng_dlp(transits(2),2))./fps;
                
                % get stimulus trace in stimulus frame rate
                stim_refframe_cam = interp1(taxis_stim, mean(stim.stimulus.color,1), resp.t_axis_trial,'pchip');  
                
                % generate rasters of 4 randomly selected cells
                permvec = randperm(ncells);
                permvec = permvec(1:4);
                clxx = [linspace(1,0,nrepeats)' zeros(1,nrepeats)' linspace(0,1,nrepeats)' ];
                pxx = get(0,'screensize');
                hf100=figure(100); set(hf100,'color','w','position',[1 1 round(pxx(3)/3)*1.5 round(pxx(3)/3)]);
                for tu=1:length(permvec)
                    subplot(2,2,tu);
                    for tv=1:nrepeats
                        plot(resp.t_axis_trial(1:sortfrmnum), ...
                            resp.cell(permvec(tu)).raster(tv,1:sortfrmnum),'-','color',clxx(tv,:));
                        hold on;
                    end
                    xlabel('Time (s)'); set(gca,'box','off','xlim',[0 taxis_stim(end)],'FontSize',10);
                    title(sprintf('Cell id %d',permvec(tu)));
                    if tu==2
                        cbb = colorbar(gca,'Ticks',[0 1],'TickLabels',[1 nrepeats],'FontSize',10);
                        cbb.Label.String = 'Trials';
                        colormap(clxx);
                    end
                end
                
                
                % generate a figure of stimulus trace and the mean response traces of
                % all cells
                alltraces = zeros(sortfrmnum,ncells);
                ylabelnames = cell(ncells,1);
                for nc=1:ncells
                    %alltraces(:,nc) = mean(resp.cell(nc).raster_filtered,1)';
                    alltraces(:,nc) = mean(resp.cell(nc).raster,1)';
                    ylabelnames{nc} = sprintf('Cell %d',nc);
                end
                numsplots = ceil(ncells/25);
                
                pxx = get(0,'screensize');
                wd = min([numsplots*750 pxx(3)]);
                cls = brewermap(ncells,'Dark2');
                %                     hf101=figure(101);
                %                     set(hf101,'position',[1 1 wd pxx(4)],'color','w');
                %                     for nsp = 1:numsplots
                %                         subplot(26,numsplots,nsp);
                %                         plot(taxis_stim, mean(stim.stimulus.color,1),'-k');
                %                         xlabel('Time (s)'); set(gca,'box','off','xlim',[0 taxis_stim(end)]);
                %                         drawnow;
                %
                %                         rng = (nsp-1)*25+1:min([25*nsp ncells]);
                %                         rngsubplot = numsplots+nsp:numsplots:26*numsplots;
                %                         rngsubplot = rngsubplot(1:length(rng));
                %                         subplot(26,numsplots,rngsubplot);
                %
                %                         ttime = resp.t_axis_trial(1:sortfrmnum);
                %                         sp = stackedplot(ttime, alltraces(:,rng));
                %                         spd = sp.DisplayLabels;
                %                         sp.DisplayLabels = ylabelnames(rng,1);
                %                         for spl=1:length(rng)
                %                             sp.LineProperties(spl).Color = cls(rng(spl),:);
                %                         end
                %                         xlabel('Time (s)');
                %                         drawnow;
                %                     end
                
                resp.alltraces_psth = alltraces;
                
                % --------------- Analysis section -----------------
                
                % Agglomerative hierarchical clustering on the full traces
%                 mxcl = 10; % max number of clusters
%                 clusT = clusterdata(alltraces','linkage','complete',...
%                     'Distance','seuclidean','maxclust',mxcl);
%                 D = pdist(alltraces','seuclidean');
%                 Z = linkage(D,'complete');
%                 clusT_duplicate = cluster(Z,'Maxclust',mxcl);
%                 leafOrder = optimalleaforder(Z,D);
%                 
%                 [clustcontainer,alltraces_clust] = deal(cell(1,mxcl));
%                 for i=1:mxcl
%                     clustcontainer{i} = find(clusT==i);
%                     alltraces_clust{i} = alltraces(:,clustcontainer{i});
%                 end


                % clustering using corrcoeff of traces as distance metric
                mxcl = 10;

                % define distance metric for use in clustering
                cf = corrcoef(alltraces);
                cf(isnan(cf)) = 0;
                cf = cf - diag(diag(cf)) + diag(zeros(ncells,1));

                % redefine cf to positive definite matrix
                cf_pos = cf - min(cf(:));
                cf_pos = cf_pos - diag(diag(cf_pos));
                maxval = max(cf_pos(:));
                for i=1:size(cf_pos,1)
                    for j=1:size(cf_pos,2)
                        if i~=j
                            cf_pos(i,j) = abs(cf_pos(i,j) - maxval);
                        end
                    end
                end

                % linear distance vector
                Y = cf_pos(triu(cf_pos)~=0)';

                % cluster
                Z = linkage(cf_pos,'average');
                clusT = cluster(Z,'MaxClust',10);
                leafOrder = optimalleaforder(Z,cf_pos);

                % assign cluster id
                [clustcontainer,alltraces_clust] = deal(cell(1,mxcl));
                for i=1:mxcl
                    clustcontainer{i} = find(clusT==i);
                    alltraces_clust{i} = alltraces(:,clustcontainer{i});
                end



                % Show figure for traces are clusters
                sstim = mean(stim.stimulus.color,1);
                if ncells < 50
                    clb = cell(1,ncells);
                    for nc=1:ncells
                        clb{nc} = sprintf('Cell %d',nc);
                    end
                    
                    hf102=figure(102); set(hf102,'color','w','position',[1 1 1800 1080]);
                    subplot(ncells+1,3,1);
                    plot(taxis_stim, sstim,'-k','linewidth',1);
                    set(gca,'ytick',[],'xtick',[],'box','off');
                    subplot(ncells+1,3,3);
                    plot(taxis_stim, sstim,'-k','linewidth',1);
                    set(gca,'ytick',[],'xtick',[],'box','off');
                    
                    ix = 4:3:(ncells+1)*3;
                    for i=1:ncells
                        subplot(ncells+1,3,ix(i));
                        plot(alltraces(:,i),'-b','linewidth',1);
                        set(gca,'xtick',[],'ytick',[],'box','off');
                        ylabel(sprintf('Cell %d',i));
                        ylh = get(gca,'ylabel');
                        ylp = get(ylh, 'Position');
                        set(ylh, 'Rotation',0, 'Position',ylp, 'VerticalAlignment','middle',...
                            'HorizontalAlignment','right');
                    end
                    cls = brewermap(mxcl,'Dark2');
                    subplot(ncells+1,3,ix+1);
                    dendlabels = cellstr(num2str((1:size(alltraces',1))', 'Cell %d'));
                    clstree = Z(end-mxcl+2, 3) - eps;
                    [H,T,perm]= dendrogram(Z,0,'Labels',dendlabels, 'Orientation','left', 'Reorder',leafOrder,...
                        'ColorThreshold',clstree);
                    
                    ix = 6:3:(ncells+1)*3;
                    cnt = 1;
                    for cl=1:mxcl
                        for i=1:size(alltraces_clust{cl},2)
                            subplot(ncells+1,3,ix(cnt));
                            yyaxis right;
                            plot(alltraces_clust{cl}(:,i),'-','color',cls(cl,:),'linewidth',1);
                            set(gca,'xtick',[],'ytick',[],'box','off','ycolor',cls(cl,:));
                            ylabel(sprintf('Cell %d (clust %d)',clustcontainer{cl}(i),cl));
                            ylh = get(gca,'ylabel');
                            ylp = get(ylh, 'Position');
                            set(ylh, 'Rotation',0, 'Position',ylp, 'VerticalAlignment','middle',...
                                'HorizontalAlignment','left');
                            yyaxis left;
                            set(gca,'xtick',[],'ytick',[],'box','off','ycolor','w');
                            cnt = cnt+1;
                        end
                    end
                    sgtitle(sprintf('# of cells = %d, # of clusters = %d',ncells,mxcl));
                end
                
                % Gather clustered traces in order of cluster groups 
                alltraces_sort = [];
                draw_white = [];
                for i=1:mxcl
                    alltraces_sort = [alltraces_sort alltraces_clust{i}];
                    draw_white = [draw_white size(alltraces_clust{i},2)];
                end
                draw_white_cumsum = [0 cumsum(draw_white)];
                ytck_lc = round(diff([0 cumsum(draw_white)])/2) + draw_white_cumsum(1:end-1);
                
                cmap = brewermap(100,'RdBu'); 
                hf103=figure(103);
                set(hf103,'position',[1 1 wd*1.2/2 wd*0.5/2],'color','w');
                subplot(1,2,1);
                imagesc(alltraces'); colormap(cmap);
                xlabel('Time (s)'); ylabel('Cell ids');
                title('Unclustered Ca response');
                subplot(1,2,2);
                imagesc(alltraces_sort'); colormap(cmap);
                xlabel('Time (s)'); ylabel('Cluster ids'); hold on;
                xlm = get(gca,'xlim');
                for dw=1:length(draw_white)
                    plot(xlm, repmat(draw_white_cumsum(dw),1,2),'--w','linewidth',2);
                end
                set(gca,'ytick',[ytck_lc],'yticklabel',1:mxcl);
                title('Clustered Ca responses');
                
                
                % append results to resp structure
                resp.cluster.pdist = D;
                resp.cluster.linkage = Z;
                resp.cluster.maxclust = mxcl;
                resp.cluster.leafOrder = leafOrder;
                resp.cluster.linkage = Z;
                resp.cluster.alltraces_psth_clust = alltraces_clust;
                resp.cluster.clustcontainer = clustcontainer;
                resp.stim = sstim; 
                resp.stim_refframe_cam = stim_refframe_cam; 
                
            end
            
        case 'RN'
            
            for q=1
                
                StimComb = [];
                
                % Set frame range for analysis
                frm_index = frame2stim.cam_frame_num_corrected;
                
                % check for initial delay
                cam_frm_num_before_stiminit = fix(frame2stim.stim_frame_time_corrected(1)*fps);
                cam_frm_rng_after_stiminit = frm_index(cam_frm_num_before_stiminit+1:end);
                
                % load stimulus frames (note: bw_images will be
                % returned at refresh interval sampling)
                
                [file, filepath] = uigetfile('*.mat','Select the movie mat file','MultiSelect','off');
                localvar = load(fullfile(filepath,file));
                strname = char(fieldnames(localvar));
                bw_images_ = localvar.(strname);
                if ~isfield(stim.stimulus,'refresh_rate')
                    stim.stimulus.refresh_rate = 60;
                end
                
                %                     try
                %                         [bw_images, stim] = get_bw_images(stim);
                %                         bw_images_ = squeeze(bw_images(:,:,1,:));
                %                     catch ME
                %                         if isa(ME, 'matlab.exception.JavaException')
                %                             [file, filepath] = uigetfile('*.mat','Select the movie mat file','MultiSelect','off');
                %                             localvar = load(fullfile(filepath,file));
                %                             bw_images_ = localvar.img_accum;
                %                             if ~isfield(stim.stimulus,'refresh_rate')
                %                                 stim.stimulus.refresh_rate = 60;
                %                             end
                %                         end
                %                     end
                
                % reshape stim images
                bw_images_ = double(bw_images_);
                nkx1 = size(bw_images_,1); nkx2 = size(bw_images_,2);
                nkx = nkx1*nkx2;
                d3 = size(bw_images_,3);
                bw_images_ = reshape(bw_images_,[nkx,d3]);
                bw_images_ = permute(bw_images_,[2 1]);
                bw_images_ = bw_images_./255; % - 0.5
                
                % length of temporal sta
                stalen = 30; % units of stimulus frames
                frame_refresh_time = (1/stim.stimulus.refresh_rate);
                image_refresh_time = stim.stimulus.refresh * frame_refresh_time; % sec
                nkt = stalen;
                
                % Set t_axis, n cells (note: any initial delay is offset)
                resp.t_axis = (1:length(cam_frm_rng_after_stiminit))./fps;
                
                % Preallocate variables
                [detrended_trace,detrended_filt_trace,sta_all] = deal([]);    
                parfor_progress(numel(neuronObj.ids));
                for nc=1:numel(neuronObj.ids)
                    
                    trace = neuronObj.C_raw(nc,cam_frm_rng_after_stiminit); % discarding responses prior to stim_delay
                    
                    % if detrend
                    if params.detrend
                        switch params.method
                            case 'running_percentile'
                                df_prctile = 20; % use low fluorescence values (excludes stimulus induced transients)
                                df_twind = 10; % seconds
                                df_fwind = df_twind*fps;
                                C_raw_baseline_drift = running_percentile(trace, df_fwind, df_prctile)';
                                
                                trace = trace - C_raw_baseline_drift;
                                %                                     if min(trace)<0
                                %                                         trace = trace - min(trace); % fluorescence values should be positive
                                %                                     end
                            case 'spline fit'
                        end
                    end
                    resp.detrended_trace(nc,:) = trace;
                    
                    if sum(trace)~=0
                        
                        % filter trace
                        resp.detrended_filt_trace(nc,:) = get_filtered_trace(trace, fps);
                        
                        % find peaks in temporal trace
                        use_trace = trace;
                        use_t_axis = resp.t_axis;
                        %[pks, pk_locs] = findpeaks(use_trace, use_t_axis);
                        
                        ck=0;
                        min_psd = 1.0;
                        while sum(ck)<=1
                            [ck, ~, ~] = deconvolveCa(use_trace, 'foopsi','ar1','smin', min_psd,...
                                'optimize_pars',true,'optimize_b',true);
                            min_psd = min_psd-0.05;
                        end
                        [pk_vals, pk_indices] = findpeaks(ck);
                        pk_indices = pk_indices+1;
                        pk_indices(pk_indices>length(use_t_axis)) = [];
                        pk_locs = use_t_axis(pk_indices);
                        
                        
                        %                             figure; set(gcf,'position',[147         285        1661         642]); hold on;
                        %                             plot(use_t_axis, use_trace,'-',use_t_axis(pk_indices), use_trace(pk_indices),'ok')
                        %                             plot(use_t_axis , ck,'-r');
                        
                        
                        
                        % bin spikes in resolution of stim interval
                        lastbin = frame_refresh_time*ceil(use_t_axis(end)/frame_refresh_time);
                        binned_spks = histcounts(pk_locs, 0:frame_refresh_time:lastbin);
                        if size(binned_spks,2)>size(binned_spks,1)
                            binned_spks = binned_spks';
                        end
                        shift = stim.stimulus.refresh; % shift binned spikes by refresh interval
                        binned_spks = [zeros(shift,1); binned_spks(shift+1:end)];
                        
                        % truncate stim if needed to match range of data
                        if length(binned_spks)<=size(bw_images_,1)
                            bw_images_temp = bw_images_(1:length(binned_spks),:);
                        end
                        
                        
                        %                             % calculate spatial and temporal sta based on spike times
                        %                             spkloc = find(binned_spks);
                        %                             weights = pk_vals./max(pk_vals);
                        %                             weights_ = zeros(size(binned_spks));
                        %                             weights_(pk_indices) = weights;
                        %                             weights = weights_; clear weights_;
                        %                             bw_images_temp = double(reshape(bw_images_temp, [28203,96,60]))./255;
                        %                             pk_indices(pk_indices<=30) = [];
                        %                             sta_xyt = zeros(96,60,30);
                        %                             for sp=1:numel(pk_indices)
                        %                                 sta_xyt = sta_xyt + permute(bw_images_temp(pk_indices(sp)-29:pk_indices(sp),:,:), [2 3 1]).*weights(pk_indices(sp));
                        %                             end
                        %                             sta_xyt = sta_xyt./numel(pk_indices);
                        
                        %[usta,ssta,vsta] = svd(permute(reshape(sta_xyt,[96*60 30]),[2 1]));
                        
                        
                        % Calculate STA
                        if sum(binned_spks)==0
                            STA = 0.5.*ones(nkt,nkx);
                        else
                            STA = simpleSTA(bw_images_temp, binned_spks, nkt);
                        end
                        sta = reshape(STA',[nkx1,nkx2,nkt]);
                        
                        % accumulate stats
                        resp.sta_all(nc).events = pk_locs;
                        resp.sta_all(nc).sta = sta;
                        %resp.sta_all(nc).stc = STC;
                        resp.sta_all(nc).dim1 = nkx1;
                        resp.sta_all(nc).dim2 = nkx2;
                        resp.sta_all(nc).dim3 = stalen;
                        
                    else
                        
                        % accumulate stats
                        resp.sta_all(nc).events = [];
                        resp.sta_all(nc).sta = [];
                        %resp.sta_all(nc).stc = STC;
                        resp.sta_all(nc).dim1 = [];
                        resp.sta_all(nc).dim2 = [];
                        resp.sta_all(nc).dim3 = [];
                        continue
                        
                    end
                    
                    
                    
                    %             figure(100); set(gcf,'position',[326 60 1351 909],'color','w');
                    %             clf(gcf);
                    %             subplot(3,2,1); imagesc(1:nkx, 1:nkt, STA); colormap('copper'); colorbar
                    %             title('sta'); xlabel('pixel'); ylabel('time lag');
                    %             subplot(3,2,3); plot(diag(ssta),'ok'); title('singular values of sta');
                    %             subplot(3,2,5); plot(usta(:,1:3)); title('top temporal components'); legend({'C1','C2','C3'},'box','off');
                    %             subplot(3,2,4); imagesc(reshape(vsta(:,1), nkx1, nkx2));
                    %             title('first spatial component'); colormap('copper'); colorbar; axis image;
                    %             subplot(3,2,6); imagesc(reshape(vsta(:,2), nkx1, nkx2));
                    %             title('second spatial component'); colormap('copper'); colorbar; axis image;
                    %             subplot(3,2,2); imagesc(usta(:,1:2)*ssta(1:2,1:2)*vsta(:,1:2)'); colormap('copper'); colorbar;
                    %             title('rank-2 reconstruction'); xlabel('pixel'); ylabel('time lag');
                    %             set(findall(gcf,'-property','FontSize'),'FontSize',12);
                    %             sgtitle(sprintf('cell %d',nc));
                    %             pause();
                    
                    
                    
                    %             % deconvolve trace to extract spike count
                    %             options_new.type = 'ar1';
                    %             options_new.method = 'constrained';
                    %             options_new.smin = neuronObj.options.deconv_options.smin;
                    %             options_new.optimize_pars = true;
                    %             options_new.optimize_b = true;
                    %             options_new.lambda = [];%neuronObj.options.deconv_options.lambda;
                    %             options_new.max_tau = neuronObj.options.deconv_options.max_tau;
                    %             [c, s, options] = deconvolveCa(trace, options_new);
                    %             figure; plot(1:length(trace), trace,'-',1:length(s),s,'-');
                    %
                    %             % adjust length
                    %             t1 = frame2stim.cam_frame_time_corrected';
                    %             y1 = s;
                    %             if length(t1)~=length(y1)
                    %                 ln = min([length(t1) length(y1)]);
                    %                 t1 = t1(1:ln);
                    %                 y1 = y1(1:ln);
                    %             end
                    
                    %             % calculate spatial and temporal sta based on spike times
                    %             stalen = 30; % units of stimulus frames
                    %             stawind = stalen*(1/60);
                    %             spkloc = find(y1~=0);
                    %             spkcnt = y1(spkloc);
                    %             normresp = max(spkcnt);
                    %             spktm = t1(spkloc);
                    %             spktm(spktm<=stawind) = [];
                    %             spkcnt(spktm<=stawind) = [];
                    %             sta_xyt = zeros(stim.stimulus.field_width, stim.stimulus.field_height, stalen, length(spkcnt),'uint8');
                    %             parfor sp=1:numel(spktm)
                    %                 stl = find(frame2stim.stim_frame_time_corrected<spktm(sp),1,'last');
                    %                 weight = spkcnt(sp)/normresp;
                    %                 sta_xyt(:,:,:,sp) = bw_images(:,:,stl-stalen:stl-1).*weight;
                    %             end
                    %             sta_xyt = squeeze(mean(sta_xyt,4));
                    
                    %fprintf(sprintf('%d/%d \n',nc,numel(neuronObj.ids)));
                    parfor_progress();                    
                end
                parfor_progress(0);                             
            end
            
    end
    
    stimout = stim;
    fprintf(sprintf('%s Responses extracted for %s stimulus %s\n',repelem('-',30),stim.type,repelem('-',30)));
    
end
