function [resp , StimComb, stim] = get_resp_properties(stim, neuronObj, frame2stim, varargin)


d3 = size(neuronObj.C,2);   % frames
global fps;

p = inputParser;
addParameter(p,'detrend',true, @islogical);
addParameter(p,'method','running_percentile', @ischar);
addParameter(p,'BWstim_to_camera_img_Tform',[], @ischar);
addParameter(p,'cellids',1:numel(neuronObj.ids), @isnumeric);
addParameter(p,'BW_truestim_img',[], @ischar); % sec 
addParameter(p,'AOI',[], @isstruct);  
addParameter(p,'FOV',[], @isstruct);  
addParameter(p,'medZim',[], @isnumeric);  
addParameter(p,'median_im',[], @isnumeric);  
addParameter(p,'micron_per_pixel',[], @isstruct); 
addParameter(p,'DSI_thresh',0.3, @isnumeric); 
addParameter(p,'OSI_thresh',0.3, @isnumeric); 
addParameter(p,'sta_t_wind',[-1.5 0.5], @isnumeric);  
addParameter(p,'spatial_filt_sd',50, @isnumeric);  % micron (for smoothing spatial RF) 
addParameter(p,'spatial_wind',200, @isnumeric);  % micron (for estimating temporal RF) 
addParameter(p,'fps',10, @isnumeric);  % micron (for estimating temporal RF) 
parse(p,varargin{:});
params = p.Results;

if ~exist('fps','var') || isempty(fps)
    fps = params.fps; 
end


fprintf(sprintf('Analyzing response to %s stimulus .......\n',stim.type));


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


                ncells = numel(params.cellids);
                nrepeats = stim.repetitions;
                resp.t_axis = (1:d3)./fps;
                resp.t_axis_trial = (1:diff(frame2stim.frame_rng(1,:))+1)./fps;

                S = [find(strcmp(fieldnames(stim.combinations),'delta')) find(strcmp(fieldnames(stim.combinations),'DELTA'))];
                unique_speed = unique(StimComb(:,S));
                W = [find(strcmp(fieldnames(stim.combinations),'bar_width')) find(strcmp(fieldnames(stim.combinations),'BAR_WIDTH'))];
                unique_width = unique(StimComb(:,W));
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
                if length(external_trig_times)<frame2stim.frame_rng(end,2)
                    numfr_missing = abs(length(external_trig_times)-frame2stim.frame_rng(end,2));
                    lntrial = mean(diff(frame2stim.frame_rng,[],2));
                    fprintf(sprintf('%1.2f fraction of the last trial missing.\n', numfr_missing/lntrial));
                    fprintf(sprintf('Reducing repetitions from %d >>>> %d.\n', stim.repetitions, stim.repetitions-1));
                    nrepeats = nrepeats -1;
                    stim.repetitions_corrected = nrepeats;
                end


                % set calcium raster for cells
                resp.detrended_trace = zeros(ncells,d3);
                local_t_axis = resp.t_axis;
                if numel(local_t_axis)>numel(neuronObj.C_raw(1,:))
                    local_t_axis = local_t_axis(1:length(neuronObj.C_raw(1,:)));
                end
                
                [detrended_trace_local,detrended_inferred_ca_local,detrended_inferred_spk_local,...
                    detrended_filt_trace_local] = deal(zeros(ncells, length(local_t_axis)));
                
                parfor_progress(length(params.cellids)); 
                parfor nc=params.cellids
           
                    trace = neuronObj.C_raw(nc,:);
                    
                    % if detrend
                    if params.detrend
                        switch params.method
                            case 'running_percentile'
                                df_prctile = 20; % use low fluorescence values (excludes stimulus induced transients)
                                df_twind = 30; % seconds
                                df_fwind = df_twind*fps;
                                C_raw_baseline_drift = running_percentile(trace, df_fwind, df_prctile)';

                                trace = trace - C_raw_baseline_drift;
                                %                             if min(trace)<0
                                %                                 trace = trace - min(trace); % fluorescence values should be positive
                                %                             end
                            case 'spline fit'
                        end
                    end
                    detrended_trace_local(nc,:) = trace;

                    % deconvolve again (based on new parameters)
                    lambda = 2;
                    g = 0.90;
                    [c_oasis, s_oasis, options] = deconvolveCa(detrended_trace_local(nc,:), 'ar1',g,   'foopsi', 'lambda', lambda, ...
                        'optimize_pars');
                    detrended_inferred_ca_local(nc,:) = c_oasis';
                    detrended_inferred_spk_local(nc,:) = s_oasis';



                    % filter trace
                    detrended_filt_trace_local(nc,:) = get_filtered_trace(trace, fps);
                    parfor_progress; 
                end
                parfor_progress(0); 
                
                resp.t_axis = local_t_axis; 
                resp.detrended_trace = detrended_trace_local; 
                resp.detrended_inferred_ca = detrended_inferred_ca_local;
                resp.detrended_inferred_spk = detrended_inferred_spk_local;
                resp.detrended_filt_trace = detrended_filt_trace_local;
                    
                for cc=1:length(params.cellids)
                    nc = params.cellids(cc);  
                    
                    % extract Ca and spike rasters for stimulus combinations
                    trace2use = resp.detrended_filt_trace(nc,:);
                    spktrain2use = resp.detrended_inferred_spk(nc,:);
                    traceunfilt = resp.detrended_trace(nc,:);

                    for i=1:length(unique_speed)
                        for j=1:length(unique_width)
                            for k=1:length(unique_contrast)
                                for m=1:length(unique_dirs)
                                    rc = find(StimComb(:,D)==unique_dirs(m));
                                    runid = find(stim.trial_list==rc);
                                    for ntr=1:nrepeats
                                        resp.unique_speed(i).unique_width(j).unique_contrast(k).unique_dirs(m).cell(nc).raster_traces(ntr,:) = ...
                                            trace2use(frame2stim.frame_rng(runid(ntr),1):frame2stim.frame_rng(runid(ntr),2));

                                        resp.unique_speed(i).unique_width(j).unique_contrast(k).unique_dirs(m).cell(nc).raster_unfilt_traces(ntr,:) = ...
                                            traceunfilt(frame2stim.frame_rng(runid(ntr),1):frame2stim.frame_rng(runid(ntr),2));


                                        resp.unique_speed(i).unique_width(j).unique_contrast(k).unique_dirs(m).cell(nc).raster_inf_spk(ntr,:) = ...
                                            spktrain2use(frame2stim.frame_rng(runid(ntr),1):frame2stim.frame_rng(runid(ntr),2));
                                    end
                                end
                            end
                        end
                    end
                    fprintf(sprintf('%d/%d\n',nc,ncells));
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


                % calculate dsi from Ca traces and inferred spikes
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
                [rho_Ca,rho_spk] = deal(zeros(ncells,length(unique_dirs)));
                [angle_ds_Ca,mag_ds_Ca,angle_ds_spk,mag_ds_spk,...
                    angle_os_Ca,mag_os_Ca,angle_os_spk,mag_os_spk,...
                    ds_index_Ca, os_index_Ca, ds_index_spk, os_index_spk] = deal(zeros(1,ncells));
                [tuning_curve_ds, dirs_centered,tuning_curve_os,ors_centered] = deal(zeros(ncells, length(unique_dirs))); 
                dF_F0 = cell(ncells,1);
                for nc=1:ncells
                    % Case 1: Ca traces
                    traparea = zeros(1,length(unique_dirs));  % area under curve (keep DC response: heterogeneous across directions)
                    dF_F0_temp = zeros(length(repeats2use),length(unique_dirs));
                    for d=1:length(unique_dirs)
                        % remove baseline for each filtered trace and calculate area under curve
                        for nr=1:length(repeats2use)
                            ntr = repeats2use(nr);
                            current_trace = str2use.unique_dirs(d).cell(nc).raster_unfilt_traces(ntr,:);
                            current_trace = current_trace - min(current_trace);
                            traparea(d) = traparea(d) + sum(trapz(resp.t_axis_trial, current_trace));

                            dF_F0_temp(nr,d) = abs(max(current_trace)-prctile(current_trace,10))/prctile(current_trace,10);
                        end
                        traparea(d) = traparea(d)./length(repeats2use);
                    end
                    if sum(traparea)~=0
                        rho_Ca(nc,:) = traparea./sum(traparea);
                    end
                    [X,Y] = pol2cart(theta', rho_Ca(nc,:));
                    [XX,YY] =  pol2cart(2.*theta', rho_Ca(nc,:));
                    u = sum(X);
                    v = sum(Y);
                    uu = sum(XX);
                    vv = sum(YY);
                    [angle_ds_Ca(nc), mag_ds_Ca(nc)] = cart2pol(u,v);
                    [angle_os_Ca(nc), mag_os_Ca(nc)] = cart2pol(uu,vv);

                    dF_F0{nc} = dF_F0_temp;

                     
                    
                    
                    % Case 1: Inferred spikes
                    spkcnt_per_trial = zeros(length(unique_dirs), 1);
                    for d=1:length(unique_dirs)
                        spkcnt_per_trial(d) = numel(find(str2use.unique_dirs(d).cell(nc).raster_inf_spk(:)))/length(repeats2use);
                    end
                    if sum(spkcnt_per_trial)~=0
                        rho_spk(nc,:) = spkcnt_per_trial./sum(spkcnt_per_trial);
                    end
                    [X,Y] = pol2cart(theta', rho_spk(nc,:));
                    [XX,YY] =  pol2cart(2.*theta', rho_spk(nc,:));
                    u = sum(X);
                    v = sum(Y);
                    uu = sum(XX);
                    vv = sum(YY);
                    [angle_ds_spk(nc), mag_ds_spk(nc)] = cart2pol(u,v);
                    [angle_os_spk(nc), mag_os_spk(nc)] = cart2pol(uu,vv);


                    % Calculate direction selective index from Ca trace
                    for kj=1
                        pd = wrapTo360(angle_ds_Ca(nc).*180./pi); % preferred direction
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

                        pd_resp_Ca = mean(rho_Ca(nc,max_resp_inds).*cosd(pd_neighbor_angs));
                        nd_resp_Ca = mean(rho_Ca(nc,min_resp_inds).*cosd(nd_neighbor_angs));
                        ds_index_Ca(nc) = abs(pd_resp_Ca-nd_resp_Ca)/(pd_resp_Ca + nd_resp_Ca);

                        pd_degree_Ca(nc) = pd; nd_degree_Ca(nc) = nd;
                        
                        % Center DS tuning curve 
                        [~,ind_temp] = max(cosd(max_resp_angs - pd)); 
                        dir_near_pd = max_resp_angs(ind_temp); 
                        [~,ind_peak] = intersect(unique_dirs, dir_near_pd); 
                        %norm_resp = rho_Ca(nc,:)./max(rho_Ca(nc,:));
                        norm_resp = (rho_Ca(nc,:) - min(rho_Ca(nc,:)))./max((rho_Ca(nc,:) - min(rho_Ca(nc,:))));
                        [~,lc] = max(norm_resp); 
                        Di = length(unique_dirs); 
                        if ind_peak~=lc
                            if (ind_peak>1 && ind_peak<Di && ind_peak>lc) || (ind_peak==1 && lc==Di) || (ind_peak==Di && ind_peak>lc)
                                norm_resp = circshift(norm_resp,1); 
                            elseif (ind_peak>1 && ind_peak<Di && ind_peak<lc) || (ind_peak==1 && lc<lc) || (ind_peak==Di && lc==1)
                                norm_resp = circshift(norm_resp,-1); 
                            end
                        end
                        if ind_peak~=1
                            idd = circshift(1:Di, -ind_peak+1); 
                            norm_resp = norm_resp(idd);
                        end
                        ticks = 1:Di; 
                        dirs_sorted_centered = wrapTo360(circshift(unique_dirs, round(Di/2)-1)); 
                        ticks_centered = circshift(ticks,round(Di/2)-1); 
                        norm_resp = norm_resp(ticks_centered);  
                        tuning_curve_ds(nc,:) = norm_resp; 
                        temp_var = wrapTo180(dirs_sorted_centered); 
                        dirs_centered(nc,:) = temp_var;  
                    end

                    % Calculate orientation selective index from Ca trace
                    for kj=1

                        po = wrapTo360(angle_os_Ca(nc).*180./pi); % preferred orientation
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

                        po_resp_Ca = mean(rho_Ca(nc,max_resp_inds).*cosd(po_neighbor_angs));
                        no_resp_Ca = mean(rho_Ca(nc,min_resp_inds).*cosd(no_neighbor_angs));
                        ortho1_resp = mean(rho_Ca(nc,ortho1_resp_inds).*cosd(ortho1_neighbor_angs));
                        ortho2_resp = mean(rho_Ca(nc,ortho2_resp_inds).*cosd(ortho2_neighbor_angs));

                        %os_index_Ca(nc) = abs(1 - mean([ortho1_resp ortho2_resp])/po_resp_Ca);
                        os_index_Ca(nc) = abs(po_resp_Ca -  mean([ortho1_resp ortho2_resp])) / (po_resp_Ca + mean([ortho1_resp ortho2_resp]));

                        nume = hypot(sum(rho_Ca(nc,:).*cosd(unique_dirs')), sum(rho_Ca(nc,:).*sind(unique_dirs')));
                        deno = sum(rho_Ca(nc,:));
                        %os_index(nc) = (nume/deno);

                        po_degree_Ca(nc) = po; no_degree_Ca(nc) = no;
                        
                        
                        
                        % Center OS tuning curve 
                        [~,ind_temp] = max(cosd(max_resp_angs - po)); 
                        dir_near_po = max_resp_angs(ind_temp); 
                        [~,ind_peak] = intersect(unique_dirs, dir_near_po); 
                        %norm_resp = rho_Ca(nc,:)./max(rho_Ca(nc,:));
                        norm_resp = (rho_Ca(nc,:) - min(rho_Ca(nc,:)))./max((rho_Ca(nc,:) - min(rho_Ca(nc,:))));
                        [~,lc] = max(norm_resp); 
                        Di = length(unique_dirs); 
                        if ind_peak~=lc
                            if (ind_peak>1 && ind_peak<Di && ind_peak>lc) || (ind_peak==1 && lc==Di) || (ind_peak==Di && ind_peak>lc)
                                norm_resp = circshift(norm_resp,1); 
                            elseif (ind_peak>1 && ind_peak<Di && ind_peak<lc) || (ind_peak==1 && lc<lc) || (ind_peak==Di && lc==1)
                                norm_resp = circshift(norm_resp,-1); 
                            end
                        end
                        if ind_peak~=1
                            idd = circshift(1:Di, -ind_peak+1); 
                            norm_resp = norm_resp(idd);
                        end
                        ticks = 1:Di; 
                        ors_sorted_centered = wrapTo360(circshift(unique_dirs, round(Di/2)-1)); 
                        ticks_centered = circshift(ticks,round(Di/2)-1); 
                        norm_resp = norm_resp(ticks_centered);  
                        tuning_curve_os(nc,:) = norm_resp; 
                        temp_var = wrapTo180(ors_sorted_centered); 
                        ors_centered(nc,:) = temp_var; 
                        
                    end

                    % Calculate direction selective index from Spk train
                    for kj=1
                        pd = wrapTo360(angle_ds_spk(nc).*180./pi); % preferred direction
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

                        pd_resp_spk = mean(rho_spk(nc,max_resp_inds).*cosd(pd_neighbor_angs));
                        nd_resp_spk = mean(rho_spk(nc,min_resp_inds).*cosd(nd_neighbor_angs));
                        ds_index_spk(nc) = abs(pd_resp_spk-nd_resp_spk)/(pd_resp_spk + nd_resp_spk);


                        pd_degree_spk(nc) = pd; nd_degree_spk(nc) = nd;
                    end

                    % Calculate orientation selective index from Spk train
                    for kj=1
                        po = wrapTo360(angle_os_spk(nc).*180./pi); % preferred orientation
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

                        po_resp_spk = mean(rho_spk(nc,max_resp_inds).*cosd(po_neighbor_angs));
                        no_resp_spk = mean(rho_spk(nc,min_resp_inds).*cosd(no_neighbor_angs));
                        ortho1_resp = mean(rho_spk(nc,ortho1_resp_inds).*cosd(ortho1_neighbor_angs));
                        ortho2_resp = mean(rho_spk(nc,ortho2_resp_inds).*cosd(ortho2_neighbor_angs));

                        %os_index_spk(nc) = abs(1 - mean([ortho1_resp ortho2_resp])/po_resp_spk);
                        os_index_spk(nc) = abs(po_resp_spk -  mean([ortho1_resp ortho2_resp])) / (po_resp_spk + mean([ortho1_resp ortho2_resp]));

                        nume = hypot(sum(rho_spk(nc,:).*cosd(unique_dirs')), sum(rho_spk(nc,:).*sind(unique_dirs')));
                        deno = sum(rho_spk(nc,:));

                        po_degree_spk(nc) = po; no_degree_spk(nc) = no;
                    end
                end


                resp.repeats_used_for_dsi_estimate = repeats2use;

                resp.ds_estimates_Ca.rho = rho_Ca;
                resp.ds_estimates_Ca.theta = theta;
                resp.ds_estimates_Ca.mag = mag_ds_Ca;
                resp.ds_estimates_Ca.angle = angle_ds_Ca;
                resp.ds_estimates_Ca.ds_index = ds_index_Ca;
                resp.ds_estimates_Ca.pd_degree = pd_degree_Ca;
                resp.ds_estimates_Ca.nd_degree = nd_degree_Ca;            
                resp.ds_estimates_Ca.tuning_curve = tuning_curve_ds;
                resp.ds_estimates_Ca.dirs_centered = dirs_centered;
                


                resp.os_estimates_Ca.rho = rho_Ca;
                resp.os_estimates_Ca.theta = theta;
                resp.os_estimates_Ca.mag = mag_os_Ca;
                resp.os_estimates_Ca.angle = angle_os_Ca;
                resp.os_estimates_Ca.os_index = os_index_Ca;
                resp.os_estimates_Ca.po_degree = po_degree_Ca;
                resp.os_estimates_Ca.no_degree = no_degree_Ca;
                resp.os_estimates_Ca.tuning_curve = tuning_curve_os;
                resp.os_estimates_Ca.ors_centered = ors_centered;

                resp.ds_estimates_spk.rho = rho_spk;
                resp.ds_estimates_spk.theta = theta;
                resp.ds_estimates_spk.mag = mag_ds_spk;
                resp.ds_estimates_spk.angle = angle_ds_spk;
                resp.ds_estimates_spk.ds_index = ds_index_spk;
                resp.ds_estimates_spk.pd_degree = pd_degree_spk;
                resp.ds_estimates_spk.nd_degree = nd_degree_spk;

                resp.os_estimates_spk.rho = rho_spk;
                resp.os_estimates_spk.theta = theta;
                resp.os_estimates_spk.mag = mag_os_spk;
                resp.os_estimates_spk.angle = angle_os_spk;
                resp.os_estimates_spk.os_index = os_index_spk;
                resp.os_estimates_spk.po_degree = po_degree_spk;
                resp.os_estimates_spk.no_degree = no_degree_spk;



                % select ds cells 
                hf1 = figure(); clf;
                title('Select DS cells');
                if size(ds_index_Ca,1)==1 || size(ds_index_Ca,2)==1
                    plot(1:ncells, ds_index_Ca,'ok','markerfacecolor','k'); hold on;
                    xll = get(gca,'xlim');
                    xlabel('Cell #'); ylabel('DSI'); set(gca,'xlim',[xll(1)-10 xll(2)+10]);
                    dim_ = 1;
                else
                    plot(log10(ds_index_Ca(:,1)), log10(ds_index_Ca(:,2)), 'ok','markerfacecolor','k'); hold on;
                    xlabel('log10(parameter 1)'); ylabel('log10(parameter 2)');
                    xl = min([log10(ds_index_Ca(:,1)); log10(ds_index_Ca(:,2))]);
                    set(gca,'xlim',[xl*1.2 0],'ylim',[xl*1.2 0]);
                    axis square;
                    dim_ = 2;
                end
                xi = get(gca,'xlim'); yi = get(gca,'ylim'); 
                if isempty(params.DSI_thresh)
                    params.DSI_thresh = 0.3; 
                end
                x = [xi fliplr(xi) xi(1)]; y = [params.DSI_thresh params.DSI_thresh yi(2) yi(2) params.DSI_thresh]; 
                plot(x,y,'.-r','markersize',8,'linewidth',1);
%                 [x,y] = deal([]);
%                 but = 1;
%                 while ~isempty(but)
%                     [x_,y_,but] = ginput(1);
%                     x = [x x_];
%                     y = [y y_];
%                     plot(x,y,'.-r','markersize',8,'linewidth',1);
%                 end
                if dim_==1
                    IN = inpolygon(1:ncells, ds_index_Ca,x,y);
                    I = find(IN==1);
                    X = 1:ncells;
                    plot(X(I), ds_index_Ca(I),'ko','markerfacecolor','r'); set(gca,'xlim',[-1 inf]);
                elseif dim_==2
                    IN = inpolygon(log10(ds_index_Ca(:,1)), log10(ds_index_Ca(:,2)),x,y);
                    I = find(IN==1);
                    plot(log10(ds_index_Ca(I,1)), log10(v(I,2)),'ko','markerfacecolor','r');
                end
                resp.ds_ids = I;



                % select os cells
                hf2 = figure(); clf;
                title('Select DS cells');
                if size(os_index_Ca,1)==1 || size(os_index_Ca,2)==1
                    plot(1:ncells, os_index_Ca,'ok','markerfacecolor','k'); hold on;
                    xll = get(gca,'xlim');
                    xlabel('Cell #'); ylabel('DSI'); set(gca,'xlim',[xll(1)-10 xll(2)+10]);
                    dim_ = 1;
                else
                    plot(log10(os_index_Ca(:,1)), log10(os_index_Ca(:,2)), 'ok','markerfacecolor','k'); hold on;
                    xlabel('log10(parameter 1)'); ylabel('log10(parameter 2)');
                    xl = min([log10(os_index_Ca(:,1)); log10(os_index_Ca(:,2))]);
                    set(gca,'xlim',[xl*1.2 0],'ylim',[xl*1.2 0]);
                    axis square;
                    dim_ = 2;
                end
                xi = get(gca,'xlim'); yi = get(gca,'ylim');                 
                if isempty(params.OSI_thresh)
                    params.OSI_thresh = 0.3; 
                end
                x = [xi fliplr(xi) xi(1)]; y = [params.OSI_thresh params.OSI_thresh yi(2) yi(2) params.OSI_thresh]; 
                plot(x,y,'.-r','markersize',8,'linewidth',1);
%                 [x,y] = deal([]);
%                 but = 1;
%                 while ~isempty(but)
%                     [x_,y_,but] = ginput(1);
%                     x = [x x_];
%                     y = [y y_];
%                     plot(x,y,'.-r','markersize',8,'linewidth',1);
%                 end
                if dim_==1
                    IN = inpolygon(1:ncells, os_index_Ca,x,y);
                    I = find(IN==1);
                    X = 1:ncells;
                    plot(X(I), os_index_Ca(I),'ko','markerfacecolor','r'); set(gca,'xlim',[-1 inf]);
                elseif dim_==2
                    IN = inpolygon(log10(os_index_Ca(:,1)), log10(os_index_Ca(:,2)),x,y);
                    I = find(IN==1);
                    plot(log10(os_index_Ca(I,1)), log10(os_index_Ca(I,2)),'ko','markerfacecolor','r');
                end
                resp.os_ids = I;




            end

        case 'CH'

            for q=1

                StimComb = [];

                s = fieldnames(stim);
                ncells = numel(params.cellids);
                resp.t_axis = (1:d3)./fps;
                resp.t_axis_trial = resp.t_axis(1:frame2stim.frame_rng_cam(find(frame2stim.trial==1,1,'last'),2));
                

                % correction for nrepeats (if image acquisition stopped before stim ended)
                transits = [0; find(diff(frame2stim.trial)~=0)];
                transits = [[1; transits(2:end)+1] [transits(2:end); length(frame2stim.trial)]]; 
                fr_transits = [frame2stim.frame_rng_cam(transits(:,1),1) frame2stim.frame_rng_cam(transits(:,2),2)]; 
                numfrs_per_trial = diff(fr_transits,[],2)+1; 
                cumfr_trials = cumsum(numfrs_per_trial);
                
                transits_local = find(diff(frame2stim.trial)~=0); 
                numframes_lastphase = diff(frame2stim.frame_rng_cam(transits_local(2),1:2)); % last phase of 1st trial 
                buffer = ceil(numframes_lastphase.*0.5); 
                
                if numfrs_per_trial(end)<numfrs_per_trial(1)-buffer 
                    missing_frnum = numfrs_per_trial(end) - numfrs_per_trial(1); 
                    fprintf(sprintf('More than 0.5 of last phase of last trial missing.\n'));
                    fprintf(sprintf('Reducing repetitions from %d >>>> %d.\n', numel(numfrs_per_trial), numel(numfrs_per_trial)-1));
                    nrepeats = numel(numfrs_per_trial)-1; 
                else 
                    missing_frnum = numfrs_per_trial(end) - numfrs_per_trial(1); 
                    fprintf(sprintf('%d frames of last phase of last trial missing. Retaining %d trials.\n',missing_frnum,numel(numfrs_per_trial)));
                    nrepeats = numel(numfrs_per_trial); 
                end
                
                    
                % correction for nrepeats (if neuron trace is shorter than frame triggers - not sure why that would happen)
                if length(neuronObj.C_raw(1,:))<cumfr_trials(end)
                    missing_frnum = cumfr_trials(end) - length(neuronObj.C_raw(1,:)); 
                    if missing_frnum>=buffer 
                        fprintf(sprintf('More than 0.5 of last phase of last trial missing.\n'));
                        fprintf(sprintf('Reducing repetitions from %d >>>> %d.\n', numel(numfrs_per_trial), numel(numfrs_per_trial)-1));
                        nrepeats = numel(numfrs_per_trial)-1;  
                    else 
                        missing_frnum = numfrs_per_trial(end) - numfrs_per_trial(1); 
                        fprintf(sprintf('%d frames of last phase of last trial missing. Retaining %d trials.\n',missing_frnum,numel(numfrs_per_trial)));
                        nrepeats = numel(numfrs_per_trial); 
                    end
                end

                % restructure trial and frame range based on updated
                % nrepeats 
                numfrs_per_trial = numfrs_per_trial(1:nrepeats); 
                cumfr_trials = cumfr_trials(1:nrepeats); 
                fr_transits = fr_transits(1:nrepeats,:); 
                sortfrmnum = min(numfrs_per_trial);
                diff_fr = diff(fr_transits,[],2)+1 - sortfrmnum; 
                fr_transits_touse = [fr_transits(:,1) fr_transits(:,2)-diff_fr]; 
                

                % Extract trial responses of cells 
                for cc=1:ncells
                    
                    nc = params.cellids(cc); 

                    trace = neuronObj.C_raw(nc,:);
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
                                %                             if min(trace)<0
                                %                                 trace = trace - min(trace); % fluorescence values should be positive
                                %                             end
                            case 'spline fit'
                        end
                    end
                    resp.detrended_trace(nc,:) = trace;


                    % deconvolve again (based on new parameters)
                    lambda = 2.4;
                    [c_oasis, s_oasis, options] = deconvolveCa(resp.detrended_trace(nc,:), 'ar1', 'foopsi', 'lambda', lambda, ...
                        'optimize_pars');
                    resp.detrended_inferred_ca(nc,:) = c_oasis';
                    resp.detrended_inferred_spk(nc,:) = s_oasis';



                    % filter trace
                    if sum(trace)~=0
                        resp.detrended_filt_trace(nc,:) = get_filtered_trace(trace, fps);
                    else
                        resp.detrended_filt_trace(nc,:) = trace;
                    end
                    
                    % extract rasters for stimulus combinations
                    trace2use = resp.detrended_trace(nc,:); 
                    trace2use_filt = resp.detrended_filt_trace(nc,:); 
                    for ntr=1:nrepeats
                        resp.cell(nc).raster(ntr,:) = trace2use(fr_transits_touse(ntr,1):fr_transits_touse(ntr,2));
                        resp.cell(nc).raster_filtered(ntr,:) = trace2use_filt(fr_transits_touse(ntr,1):fr_transits_touse(ntr,2));
                    end

                end


                % get psth from raw traces and filtered raw traces
                [tracefilt_psth, trace_psth, trace_psth_std, tracefilt_psth_std] = deal(zeros(sortfrmnum,ncells));
                for nc=1:ncells
                    trace_psth(:,nc) = mean(resp.cell(nc).raster,1)';
                    trace_psth_std(:,nc) = std(resp.cell(nc).raster,[],1)'./2;
                    tracefilt_psth(:,nc) = mean(resp.cell(nc).raster_filtered,1)';
                    tracefilt_psth_std(:,nc) = std(resp.cell(nc).raster_filtered,[],1)'./2;
                end

                % get corresponding stimulus trace in camera rate
                taxis_stim = linspace(1,length(1:frame2stim.frame_rng_cam(transits(2),2)),...
                    length(1:frame2stim.frame_rng_dlp(transits(2),2)))./fps;

                % Extract rasters for deconvolved spikes
                for nc=1:ncells
                    trace2use = resp.detrended_inferred_spk(nc,:);
                    for ntr=1:nrepeats
                        tempvec = trace2use(frame2stim.frame_rng_cam(transits(ntr)+1,1):...
                            frame2stim.frame_rng_cam(transits(ntr+1),2));
                        resp.cell(nc).spk_raster{ntr} = find(tempvec)/fps;
                    end
                end

                % Mean stimulus
                sstim = mean(stim.stimulus.color,1);


                % append results to resp structure
                resp.t_stim = taxis_stim;
                resp.y_stim = sstim;
                resp.nrepeats = nrepeats; 
                resp.trace_psth = trace_psth;
                resp.trace_psth_std = trace_psth_std;
                resp.tracefilt_psth = tracefilt_psth;
                resp.tracefilt_psth_std = tracefilt_psth_std;
                if exist('D','var'); resp.cluster.pdist = D; end
                if exist('Z','var'); resp.cluster.linkage = Z; end
                if exist('mxcl','var'); resp.cluster.maxclust = mxcl; end
                if exist('leafOrder','var'); resp.cluster.leafOrder = leafOrder; end
                if exist('alltraces_clust','var'); resp.cluster.alltraces_psth_clust = alltraces_clust; end
                if exist('clustcontainer','var'); resp.cluster.clustcontainer = clustcontainer; end

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


                ncells = numel(params.cell_ids);
                nrepeats = stim.repetitions;
                resp.t_axis = (1:d3)./fps;
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
                if length(external_trig_times)<frame2stim.frame_rng(end,2)
                    numfr_missing = abs(length(external_trig_times)-frame2stim.frame_rng(end,2));
                    lntrial = mean(diff(frame2stim.frame_rng,[],2));
                    fprintf(sprintf('%1.2f fraction of the last trial missing.\n', numfr_missing/lntrial));
                    fprintf(sprintf('Reducing repetitions from %d >>>> %d.\n', stim.repetitions, stim.repetitions-1));
                    nrepeats = nrepeats -1;
                    stim.repetitions_corrected = nrepeats;
                end

                % set raster for cells
                resp.detrended_trace = zeros(ncells,d3);
                for cc=1:ncells
                    
                    nc = params.cellids(cc); 

                    trace = neuronObj.C_raw(nc,:);
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
                                if min(trace)<0
                                    trace = trace - min(trace); % fluorescence values should be positive
                                end
                            case 'spline fit'
                        end
                    end
                    resp.detrended_trace(nc,:) = trace;

                    % filter trace
                    resp.detrended_filt_trace(nc,:) = get_filtered_trace(trace, fps);


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
                                            trace2use(frame2stim.frame_rng(runid(ntr),1):frame2stim.frame_rng(runid(ntr),2));

                                        resp.unique_tp(i).unique_sp(j).unique_contrast(k).unique_dirs(m).cell(nc).raster_unfilt_traces(ntr,:) = ...
                                            traceunfilt(frame2stim.frame_rng(runid(ntr),1):frame2stim.frame_rng(runid(ntr),2));
                                    end
                                end
                            end
                        end
                    end
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
                for nc=1:ncells
                    traparea = zeros(1,length(unique_dirs));  % area under curve (keep DC response: heterogeneous across directions)
                    for d=1:length(unique_dirs)

                        % remove baseline for each filtered trace and calculate area under curve
                        for nr=2:length(repeats2use)
                            ntr = repeats2use(nr);
                            current_trace = str2use.unique_dirs(d).cell(nc).raster_unfilt_traces(ntr,:);
                            current_trace = current_trace;
                            traparea(d) = traparea(d) + sum(trapz(resp.t_axis_trial, current_trace));
                        end
                        traparea(d) = traparea(d)./(nrepeats-1);
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

        case 'RN'

            for q=1

                StimComb = [];

                s = fieldnames(stim);
                ncells = numel(params.cellids);
                nfr = frame2stim.cam_frame_num_corrected(end);
                resp.t_axis = (1:nfr)./fps;

                % set refresh time
                frame_refresh_time = (1/stim.stimulus.refresh_rate);
                image_refresh_time = stim.stimulus.refresh * frame_refresh_time; % sec
                
                % sta time window 
                twind = abs(diff(params.sta_t_wind)); 
                nimage_sta = floor(twind./image_refresh_time); 
                
                % set micron per pixel
                if size(params.median_im,1)==576
                    pix2mu = params.micron_per_pixel.at_576; 
                elseif size(params.median_im,1)==1152
                    pix2mu = params.micron_per_pixel.at_1152;
                elseif size(params.median_im,1)==2304
                    pix2mu = params.micron_per_pixel.at_1152;
                end
                 

                % load stimulus frames
                try
                    [bw_images, stim] = get_bw_images(stim);
                    bw_images_ = squeeze(bw_images(:,:,1,:));
                catch ME
                    warning('Error getting images from movie xml file.');
                    [file, filepath] = uigetfile('*.mat','Select the movie mat file','MultiSelect','off');
                    localvar = load(fullfile(filepath,file));
                    strname = char(fieldnames(localvar));
                    bw_images_ = localvar.(strname);
                end


                % get transformation matrix if mapping image to stimulus
                mapreply = questdlg('Would you like to map BW stim to imaging FOV?', ...
                    'Mapping','Yes','No','No');
                if strcmpi(mapreply,'Yes')
                    mapexistreply = questdlg('Does the mapping file exist?', ...
                        'Mapping file','Yes','No','No');
                    if strcmpi(mapexistreply,'Yes')
                        tform_struct = struct;
                        if exist(params.BWstim_to_camera_img_Tform,'file')
                            load(params.BWstim_to_camera_img_Tform);
                        else
                            [file, filepath] = uigetfile('*.mat','Select the transformation file','MultiSelect','off');
                            if isempty(file) || isempty(filepath)
                                warning('No transformation file found. Will default to unmapped BW stim.');
                            else
                                load(fullfile(filepath,file));
                            end
                        end
                    else
                        [file, filepath] = uigetfile('*.tif','Select the photographic mapping file','MultiSelect','off');
                        photographic_mapping = imread([filepath,file]);
                        if size(photographic_mapping,1)>576 || size(photographic_mapping,2)>576 % conserving memory 
                            photographic_mapping = imresize(photographic_mapping, [576 576]);
                        end   
                        if isempty(params.BW_truestim_img)
                            stimimg = bw_images_(:,:,1); 
                        else 
                            stimimg = load(params.BW_truestim_img); 
                            fldname = fieldnames(stimimg); 
                            stimimg = stimimg.(fldname{1}); 
                        end
                        [movingPoints,fixedPoints] = cpselect(stimimg, photographic_mapping, 'Wait', true);
                        tform = fitgeotrans(movingPoints,fixedPoints,'affine');
                        
                        tform_struct.Ifixed = photographic_mapping;
                        tform_struct.Imoving = stimimg;
                        tform_struct.moving_points = movingPoints;
                        tform_struct.fixedPoints = fixedPoints;
                        tform_struct.tform = tform;
                        uisave('tform_struct','tform.mat');
                        
%                         % Show mapping figures 
%                         temp_loc = [ceil(size(stimimg,2)/2) ceil(size(stimimg,1)/2)]; % x, y
%                         [xm, ym] = transformPointsForward(tform, temp_loc(1), temp_loc(2));
%                         mapped_im = imwarp(stimimg, tform_struct.tform,'nearest',...
%                                 'OutputView',imref2d(size(photographic_mapping)),'SmoothEdges',true);
                      
                    end
                end
                
                % Transform the first image 
                bw_image_mapped = imwarp(bw_images_(:,:,1), tform_struct.tform,'nearest',...
                    'OutputView',imref2d(size(tform_struct.Ifixed)),'SmoothEdges',true);

                

                % reshape stim images
                nkx1 = size(bw_images_,1); nkx2 = size(bw_images_,2);
                nkx = nkx1*nkx2;
                d3 = size(bw_images_,3);
                bw_images_ = permute(reshape(bw_images_,[nkx,d3]),[2 1]);
                bw_images_ = bw_images_-0.5; 
               
                % check for initial delay
                init_delay = stim.stimulus.delay_frames/stim.stimulus.refresh_rate; % sec
%                 if init_delay~=0
%                     stim_fr_time = frame2stim.stim_frame_time_corrected(2:end)-init_delay;
%                 end
                
                 

                % set a limit on how long stimulus to use for calculating STA
                prompt = {sprintf('%s\n%s%3.2f%s%3.2f','Enter stimulus duration in sec you want to use for STA analysis.',...
                    'Stim last frame: ',frame2stim.stim_frame_time_corrected(end),'. Last Ca image: ',resp.t_axis(end))};
                dlgtitle = 'Duration in sec (value set to default)';
                definput = {num2str(resp.t_axis(end))};
                dims = [1 70];
                opts.Interpreter = 'tex';
                answer = inputdlg(prompt,dlgtitle,dims,definput,opts);
                stimdur2use = str2double(answer);
                
                
                resp.sta_all = [];
                for cc=1:ncells
                    
                    nc = params.cellids(cc); 
                    
                    % get com of cell (All COMs are : [Y, X])
                    com_AOI = com(neuronObj.A(:,nc), neuronObj.options.d1, neuronObj.options.d2); % [y, x]
                    com_FOV = [com_AOI(1)+min(params.AOI.y) com_AOI(2)+min(params.AOI.x)]; % [y, x]
                    [x_bw, y_bw] =  transformPointsInverse(tform_struct.tform, com_FOV(2), com_FOV(1)); % [xm, ym]
                    com_bw_frame = [y_bw x_bw];
                        
                    % get trace 
                    trace = neuronObj.C_raw(nc,1:nfr);
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
                                %                         if min(trace)<0
                                %                             trace = trace - min(trace); % fluorescence values should be positive
                                %                         end
                            case 'spline fit'
                        end
                    end
                    resp.detrended_trace(nc,:) = trace;
                    
                    % filter trace
                    if sum(trace)~=0
                        resp.detrended_filt_trace(nc,:) = get_filtered_trace(trace, fps);
                    else
                        resp.detrended_filt_trace(nc,:) = trace;
                    end
                    
                    % find peaks in temporal trace
                    use_trace = trace; % resp.detrended_filt_trace(nc,:);
                    use_trace = use_trace(resp.t_axis>init_delay);
                    if ~isempty(resp.t_axis(find(resp.t_axis<=init_delay,1,'last')))
                        use_t_axis = resp.t_axis(resp.t_axis>init_delay) - resp.t_axis(find(resp.t_axis<=init_delay,1,'last'));
                    else
                        use_t_axis = resp.t_axis;
                    end
                    %[pks, pk_locs] = findpeaks(use_trace, use_t_axis, 'MinPeakProminence',2);
                    
                    % redefine limits for response and stimulus (note, we are only
                    % using non-repeating bw images)
                    use_trace_redef = use_trace(1:find(use_t_axis<=stimdur2use,1,'last'));
                    use_t_axis_redef = use_t_axis(1:find(use_t_axis<=stimdur2use,1,'last'));
                    %use_bw_images_redef = bw_images_(1:floor(stimdur2use/image_refresh_time),:);
                    
                    
                    if sum(use_trace_redef)~=0

                        ck=0;
                        min_psd = 0.7;
                        while sum(ck)<=1 
                            [ck, sk, ~] = deconvolveCa(use_trace_redef, 'foopsi','ar1','smin', min_psd,...
                                'optimize_pars',true,'optimize_b',true);
                            min_psd = min_psd-0.05;
                        end
                        spkt = use_t_axis_redef((sk~=0)); % sec 
                        spkwt = sk(sk~=0); % weight of spike response 
                        %spkwt_norm = spkwt; % normalized weight 
                        
                        % get slope at calcium transients at upsampled trace (10x)
                        dt = diff(use_t_axis_redef(1:2))./10; 
                        use_t_axis_redef_upsamp = use_t_axis_redef(1):dt:use_t_axis_redef(end);
                        use_trace_redef_upsamp = interp1(use_t_axis_redef, use_trace_redef, use_t_axis_redef_upsamp);
                        slope = gradient(use_trace_redef_upsamp)./dt; 
                        slope_at_transient = zeros(1,length(spkt)); 
                        for nk=1:length(spkt)
                            [~,lc] = find(use_t_axis_redef_upsamp<spkt(nk),1,'last'); 
                            slope_at_transient(nk) = slope(lc); 
                        end
                        spkwt_norm = slope_at_transient; 
                        
                         
                        % adjust spk time for STA calculation 
                        lcc = find(spkt>nimage_sta.*image_refresh_time,1,'first'); 
                        spkt = spkt(lcc:end);  
                        spkwt_norm = spkwt_norm(lcc:end);  
                        
                        % calculate spatiotemporal STA 
                        STA = zeros(nkx1,nkx2,nimage_sta);
                        for nspk=1:length(spkt)
                            image_rng = floor(spkt(nspk)./image_refresh_time)-nimage_sta:floor(spkt(nspk)./image_refresh_time)-1; 
                            for k=1:nimage_sta
                                temp = reshape(bw_images_(image_rng(k),:),[nkx1 nkx2]);
                                temp = imgaussfilt(temp,1); 
                                STA(:,:,k) = STA(:,:,k)+temp.*spkwt_norm(nspk); 
                            end
                        end
                        STA = STA./length(spkt);        
                        
                        % transform STA to mapping image reference 
                        [STA_mapped,STA_mapped_filt] = deal(zeros(size(params.median_im))); 
                        for k=1:size(STA,3)
                            STA_mapped(:,:,k) = imwarp(STA(:,:,k), tform_struct.tform,'nearest',...
                                'OutputView',imref2d(size(params.median_im)),'SmoothEdges',true);
                            %STA_mapped_filt(:,:,k) = imgaussfilt(STA_mapped(:,:,k), floor(params.spatial_filt_sd/pix2mu)); 
                        end
                        STA_mapped_filt = STA_mapped; 
                        
                        
%                         % calculate spatiotemporal STA 
%                         STA = zeros(nimage_sta,nkx1*nkx2);
%                         for nspk=1:length(spkt)
%                             image_rng = floor(spkt(nspk)./image_refresh_time)-nimage_sta:floor(spkt(nspk)./image_refresh_time)-1; 
%                             STA = STA + bw_images_(image_rng,:).*spkwt_norm(nspk); 
%                         end
%                         STA = STA./length(spkt); 
%                         STA = permute(reshape(STA,[nimage_sta nkx1 nkx2]), [2 3 1]);
%                         
%                         % transform STA to mapping image reference 
%                         [STA_mapped,STA_mapped_filt] = deal(zeros(size(params.median_im))); 
%                         for k=1:size(STA,3)
%                             STA_mapped(:,:,k) = imwarp(STA(:,:,k), tform_struct.tform,'nearest',...
%                                 'OutputView',imref2d(size(params.median_im)),'SmoothEdges',true);
%                             STA_mapped_filt(:,:,k) = imgaussfilt(STA_mapped(:,:,k), floor(params.spatial_filt_sd/pix2mu)); 
%                         end
                        
                    

                        % get temporal RF 
                            % set local spatial window 
                        numpix_tRF = ceil(params.spatial_wind./pix2mu); 
                        xl = ceil(max([com_AOI(2)-ceil(numpix_tRF/2) 1])); 
                        xu = floor(min([com_AOI(2)+ceil(numpix_tRF/2) size(params.medZim,2)])); 
                        yl = ceil(max([com_AOI(1)-ceil(numpix_tRF/2) 1])); 
                        yu = floor(min([com_AOI(1)+ceil(numpix_tRF/2) size(params.medZim,1)])); 
                        tRF_spatwind_AOI.x = [xl xu xu xl xl];
                        tRF_spatwind_AOI.y = [yl yl yu yu yl];
                        tRF_spatwind_FOV.x = tRF_spatwind_AOI.x + min(params.AOI.x);
                        tRF_spatwind_FOV.y = tRF_spatwind_AOI.y + min(params.AOI.y);
                            % calculate temporal RF from spatially filtered STA     
                        rng_y = min(tRF_spatwind_FOV.y):max(tRF_spatwind_FOV.y);
                        rng_x = min(tRF_spatwind_FOV.x):max(tRF_spatwind_FOV.x);
                        for k=1:size(STA_mapped_filt,3) 
                            resp.sta_all(nc).RF_tcourse(k) = mean(mean(STA_mapped_filt(rng_y,rng_x,k))); 
                            resp.sta_all(nc).RF_taxis = params.sta_t_wind(1):image_refresh_time:floor(params.sta_t_wind(2)/image_refresh_time)*image_refresh_time; 
                        end
                        STA_mapped_local = STA_mapped(rng_y,rng_x,:); 
                        STA_mapped_filt_local = STA_mapped_filt(rng_y,rng_x,:); 
                         
                        
                        if cc==1
                            
                            [FOV_inv_map.x,FOV_inv_map.y] = transformPointsInverse(tform_struct.tform, params.FOV.x,params.FOV.y); 
                            rng_rec = [prctile(params.median_im(:),2) prctile(params.median_im(:),99)];
                           
                            hf=figure; 
                            set(hf,'Position',[133 55 1135 742],'renderer','painters'); 
                            subplot(2,4,1);
                            imshow(tform_struct.Imoving); hold on; 
                            plot(FOV_inv_map.x, FOV_inv_map.y,'-m','linewidth',1); 
                            title('BW image'); 
                            subplot(2,4,2);
                            imshow(tform_struct.Ifixed); hold on; 
                            plot(params.FOV.x, params.FOV.y,'-m','linewidth',1); 
                            title('Photograhic mapping image'); 
                            subplot(2,4,3);
                            imshow(bw_image_mapped); hold on; 
                            plot(params.FOV.x, params.FOV.y,'-m','linewidth',1); 
                            title('Mapped BW image'); 
                            subplot(2,4,5);
                            rng_rec = [prctile(params.median_im(:),2) prctile(params.median_im(:),99)];
                            imshow(params.median_im,rng_rec); hold on
                            plot(params.AOI.x, params.AOI.y,'-y','linewidth',1);
                            plot(com_FOV(2), com_FOV(1),'*m');
                            plot(tRF_spatwind_FOV.x,tRF_spatwind_FOV.y,'-r'); 
                            title('Full Ca projection image');
                            ax=subplot(2,4,6);  
                            ln_cl = [0 1 0]; %'g'
                            ln_wid = 1;
                            txt_cl = [0 1 0]; % 'g'
                            fill_alpha = 0;
                            clrmap = [linspace(0,1,100)' linspace(0,1,100)' linspace(0,1,100)'];
                            show_rois = 'none'; % 'all', or indices, or 'none' [8 10 11 13 24]
                            contour_plot_simple(neuronObj.A,params.medZim,neuronObj,...
                                'display_numbers',false,'max_thresh',0.2,'line_clr',ln_cl,'line_wid',ln_wid,...
                                'text_clr',txt_cl,'fill_alpha',fill_alpha,'clrmap',clrmap,'ha',ax,...
                                'ind_show',false,'which_rois',nc,'clust_info',[]);
                            caxis(rng_rec);
                            plot(params.AOI.x-min(params.AOI.x)+1, params.AOI.y-min(params.AOI.y)+1,'-y','linewidth',1); 
                            plot(com_AOI(2), com_AOI(1),'o','markerfacecolor','g'); 
                            plot(tRF_spatwind_AOI.x,tRF_spatwind_AOI.y,'-r');
                            title('Truncated Ca projection image');    
                            subplot(2,4,7);
                            imshow(bw_image_mapped); hold on; 
                            plot(params.AOI.x, params.AOI.y,'-y','linewidth',1);
                            plot(params.FOV.x, params.FOV.y,'-m','linewidth',1);
                            plot(com_FOV(2), com_FOV(1),'o','markerfacecolor','g');
                            plot(tRF_spatwind_FOV.x,tRF_spatwind_FOV.y,'-r'); 
                            title('Location of cell on BW image');
                            subplot(2,4,8);
                            imshow(bw_image_mapped(min(params.AOI.y):max(params.AOI.y),...
                                min(params.AOI.x):max(params.AOI.x))); hold on; 
                            plot(com_AOI(2), com_AOI(1),'o','markerfacecolor','g');
                            plot(tRF_spatwind_AOI.x,tRF_spatwind_AOI.y,'-r'); 
                            title('Location of cell on BW image');
                            
                        end
                        
                        
                    else
                        
                        STA = [];  
                        STA_mapped = [];
                        STA_mapped_local = [];
                        STA_mapped_filt = [];
                        STA_mapped_filt_local = [];
                         
                    end
                    
                    fprintf(sprintf('Cell %d/%d\n',cc,ncells));
                    
                    
                    % Store stats 
                    resp.sta_all(nc).sta = STA;  
                    resp.sta_all(nc).sta_mapped = STA_mapped;
                    resp.sta_all(nc).sta_mapped_filt = STA_mapped_filt;
                    resp.sta_all(nc).sta_mapped_local = STA_mapped_local; 
                    resp.sta_all(nc).sta_mapped_filt_local = STA_mapped_filt_local;
                    resp.sta_all(nc).com_AOI = com_AOI; % [y, x]
                    resp.sta_all(nc).com_FOV = com_FOV; % [y, x]
                    resp.sta_all(nc).com_bw_frame = com_bw_frame; % [y, x]
                    resp.sta_all(nc).sta_local_size_micron = params.spatial_wind; 
                    %resp.sta_all(nc).corr_img = corr_img;
                    
                end
            end

        case 'LS'

            for q=1
                

                StimComb = [];
                ncells = numel(params.cellids);         
                resp.t_axis = (1:d3)./fps;
                resp.t_axis_trial = (1:frame2stim.minfrpertrial)./fps;

               
                % Extract trial responses of cells 
                [gather_traces, gather_filt_traces]  = deal(zeros(ncells, size(neuronObj.C_raw,2)));  
                parfor cc=1:ncells
                    nc = params.cellids(cc);                    
                    trace = neuronObj.C_raw(nc,:);                    

                    % if detrend
                    if params.detrend
                        switch params.method
                            case 'running_percentile'
                                df_prctile = 20; % use low fluorescence values (excludes stimulus induced transients)
                                df_twind = 10; % seconds
                                df_fwind = df_twind*fps;
                                C_raw_baseline_drift = running_percentile(trace, df_fwind, df_prctile)';

                                trace = trace - C_raw_baseline_drift;
                                %                             if min(trace)<0
                                %                                 trace = trace - min(trace); % fluorescence values should be positive
                                %                             end
                            case 'spline fit'
                        end
                    end
                    gather_traces(cc,:) = trace; 

                    % filter trace
                    if sum(trace)~=0
                        gather_filt_traces(cc,:) = get_filtered_trace(trace, fps);
                    else
                        gather_filt_traces(cc,:) = trace;
                    end
                    
                end
 
                
                % Gather data and store in resp structure  
                if numel(resp.t_axis)>size(neuronObj.C_raw,2)
                    resp.t_axis = resp.t_axis(1,1:size(neuronObj.C_raw,2));
                end
                [resp.detrended_trace,resp.detrended_filt_trace] = deal(zeros(size(neuronObj.C_raw,1), size(neuronObj.C_raw,2)));  
                for ci=1:ncells
                    nci = params.cellids(ci);            
                    resp.detrended_trace(nci,:) = gather_traces(nci,:); 
                    resp.detrended_filt_trace(nci,:) = gather_filt_traces(ci,:);                     
                end
                 
                
            end
            
        case 'None'  % spontaneous response 
            
            for q=1
                
                StimComb = [];
                ncells = numel(params.cellids);
                resp.t_axis = (1:d3)./fps;
                
                
                % Extract trial responses of cells 
                [gather_traces, gather_filt_traces]  = deal(zeros(ncells, size(neuronObj.C_raw,2)));  
                parfor cc=1:ncells
                    nc = params.cellids(cc);
                    trace = neuronObj.C_raw(nc,:);
                    
                    % if detrend
                    if params.detrend
                        switch params.method
                            case 'running_percentile'
                                df_prctile = 20; % use low fluorescence values (excludes stimulus induced transients)
                                df_twind = 10; % seconds
                                df_fwind = df_twind*fps;
                                C_raw_baseline_drift = running_percentile(trace, df_fwind, df_prctile)';
                                
                                trace = trace - C_raw_baseline_drift;
                                %                             if min(trace)<0
                                %                                 trace = trace - min(trace); % fluorescence values should be positive
                                %                             end
                            case 'spline fit'
                        end
                    end
                    gather_traces(cc,:) = trace;
                    
                    % filter trace
                    if sum(trace)~=0
                        gather_filt_traces(cc,:) = get_filtered_trace(trace, fps);
                    else
                        gather_filt_traces(cc,:) = trace;
                    end
                end
                
                % Gather data and store in resp structure
                if numel(resp.t_axis)>size(neuronObj.C_raw,2)
                    resp.t_axis = resp.t_axis(1,1:size(neuronObj.C_raw,2));
                end
                [resp.detrended_trace,resp.detrended_filt_trace] = deal(zeros(size(neuronObj.C_raw,1), size(neuronObj.C_raw,2)));
                for ci=1:ncells
                    nci = params.cellids(ci);
                    resp.detrended_trace(nci,:) = gather_traces(nci,:);
                    resp.detrended_filt_trace(nci,:) = gather_filt_traces(ci,:);
                end
                   
            end

        fprintf(sprintf('Response extracted.\n'));
    end

end




