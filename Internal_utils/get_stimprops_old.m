function stim = get_stimprops(stim, dataname, initpath)

    
    exs_us = regexpi(dataname,'_');
    if numel(exs_us)>1
        exs_us = [exs_us(1)+1 exs_us(2)-1];
    else
        exs_us = [exs_us(1)+1 numel(dataname)]; 
    end     
    if numel(dataname(exs_us(1):exs_us(2)))>6
        lookfor_n = dataname(exs_us(1):exs_us(1)+1);
    elseif numel(dataname(exs_us(1):exs_us(2)))==6
        lookfor_n = dataname(exs_us(1));   
    end
    
    
    % Distinguish files containing size(movie_wnr,1) parameters and camera triggers 
    fnames = dir(stim.filepath);
    if isempty(fnames)
        ix = strfind(stim.filepath,'Experiments')-1; 
        tempstr = stim.filepath; 
        tempstr = strrep(tempstr, tempstr(1:ix), initpath); 
        fnames = dir(tempstr); 
    else
        tempstr = stim.filepath; 
    end
    lcc = find(cell2mat(cellfun(@isempty, cellfun(@(x) find(x(:)==lookfor_n), arrayfun(@(x) x(:).name,fnames,'uniformoutput', false),...
        'uniformoutput', false), 'uniformoutput', false))==0);
    stimfullfname = cell(1,1); 
    lookfor_fulln = strcat('s0',lookfor_n); 
    cnt = 1; 
    for ii=1:length(lcc)
        if ~isempty(regexpi(fnames(lcc(ii)).name,'.m')) && ~isempty(regexpi(fnames(lcc(ii)).name,'.mat')) && ...
                ~isempty(regexpi(fnames(lcc(ii)).name,lookfor_fulln))
            stimfullfname{cnt,1} = strcat(tempstr, fnames(lcc(ii)).name); 
            cnt = cnt+1; 
        end
    end
    for ii=1:length(stimfullfname)
        if ~isempty(regexpi(stimfullfname{ii},'tstamp'))
            tstid = ii; 
        end
    end
    
    % Extract stimulus specific properties 
    if strcmpi(stim.type,'MB') || strcmpi(stim.type,'MG')
       
            
            tempstr = load(stimfullfname{1}); 
            if isfield(tempstr,'stimulus')   
                if isobject(tempstr.stimulus)
                    curpwd = pwd;
                    setlocalpath = uigetdir('/Volumes/','Select stimulus class directory'); 
                    cd(setlocalpath);
                    tempstr = load(stimfullfname{1}); 
                end
            else
                warning('''stimulus'' field does not exist in stim file. Cannot proceed !\n');
                return;
            end
            pfields = fieldnames(tempstr.stimulus.shuff_parameters); 
            ntrials = length(tempstr.stimulus.shuff_parameters); 
            for nf=1:length(pfields)
                for ntr=1:ntrials
                    stim.trials(ntr).(pfields{nf}) = tempstr.stimulus.shuff_parameters(ntr).(pfields{nf}); 
                end
            end
            
            
            time_stamps = tempstr.time_stamps;
            initoffset_camtrig_stim.dt = abs(time_stamps.stim_time_stamps(1,1) - time_stamps.external_trig_times(1)); 
            if time_stamps.stim_time_stamps(1,1)<time_stamps.external_trig_times(1)
                initoffset_camtrig_stim.who_leads = 'stim';
                initoffset_camtrig_stim.who_lags = 'cam';
                time_stamps.stim_time_stamps(1,2) = time_stamps.stim_time_stamps(1,2) - initoffset_camtrig_stim.dt; 
                time_stamps.stim_time_stamps(2:end,:) = time_stamps.stim_time_stamps(2:end,:) - initoffset_camtrig_stim.dt; 
                time_stamps.stim_time_stamps(1,1) = 0; 
            else
                initoffset_camtrig_stim.who_leads = 'cam';
                initoffset_camtrig_stim.who_lags = 'stim';
            end
            trial_dur = mean(diff(time_stamps.stim_time_stamps(2:end,2:end),[],2));  
            
            % Get all stimulus combinations
            stim.combinations(1) = stim.trials(1);  
            stim.trial_list(1)=1;
            for j=2:length(stim.trials)
                t=0;
                for jj=1:length(stim.combinations)
                    if isequal(stim.combinations(jj),stim.trials(j))
                        stim.trial_list(j)=jj;
                        t=jj;
                    end    
                end
                if ~t
                    stim.combinations(length(stim.combinations)+1)=stim.trials(j);
                    stim.trial_list(j)=length(stim.combinations);
                end
            end
                              
            t=fieldnames(stim.combinations);
            for j=1:length(t)
                if length(stim.combinations(1).(t{j}))==1
                    stim.params.(t{j})=sort(unique([stim.combinations.(t{j})]));
                else
                    for cnd = 1:length(stim.combinations)
                        temp(cnd,:) = stim.combinations(cnd).(t{j});
                    end
                    stim.params.(t{j}) = unique(temp, 'rows');
                    stim.params.(t{j}) = mat2cell(stim.params.(t{j}), ones(size(stim.params.(t{j}), 1), 1), size(stim.params.(t{j}), 2))';
                 end
            end
            stim.repetitions=length(stim.trials)/length(stim.combinations);
            
            % Assign trial duration
            stim.trial_dur = trial_dur; 
            
            % Add the field values for time_stamps into stim structure 
            pfields = fieldnames(time_stamps); 
            for nf=1:length(pfields)
                stim.time_stamps.(pfields{nf}) = time_stamps.(pfields{nf}); 
            end
            
            % check for missing camera triggers at the end (NEED TO
            % IMPLEMENT) 
            
            % Return status 
            fprintf(sprintf('\n Stimulus and image frame trigger info loaded.\n')); 
            
    elseif strcmpi(stim.type,'CH')    
        
            
            tempstr = load(stimfullfname{1}); 
            if isfield(tempstr,'stimulus')   
                if isobject(tempstr.stimulus)
                    curpwd = pwd;
                    setlocalpath = uigetdir('/Volumes/','Select stimulus class directory'); 
                    cd(setlocalpath);
                    tempstr = load(stimfullfname{1}); 
                end
            else
                warning('''stimulus'' field does not exist in stim file. Cannot proceed !\n');
                return;
            end
            numframes = length(tempstr.stimulus.color); 
            
            time_stamps = tempstr.time_stamps;
            initoffset_camtrig_stim.dt = abs(time_stamps.stim_time_stamps(1,1) - time_stamps.external_trig_times(1)); 
            % Set camera as the master and adjust offset in stim time_stamps 
            if time_stamps.stim_time_stamps(1,1)<time_stamps.external_trig_times(1) 
                initoffset_camtrig_stim.who_leads = 'stim';
                initoffset_camtrig_stim.who_lags = 'cam';
            else
                initoffset_camtrig_stim.who_leads = 'cam';
                initoffset_camtrig_stim.who_lags = 'stim';
            end
             
            
            % Add the field values for time_stamps into stim structure 
            pfields = fieldnames(time_stamps); 
            for nf=1:length(pfields)
                stim.time_stamps.(pfields{nf}) = time_stamps.(pfields{nf}); 
            end
            
            % Add all fields to stim
            if isobject(tempstr.stimulus)                
                local_props = properties(tempstr.stimulus);
                for p=1:numel(local_props)
                    stim.(local_props{p}) = tempstr.stimulus.(local_props{p});
                end
                cd(curpwd); 
            else
                pfields = fieldnames(tempstr.stimulus);
                for nf=1:length(pfields)
                    stim.stimulus.(pfields{nf}) = tempstr.stimulus.(pfields{nf});
                end
            end
            
            
            
    elseif strcmpi(stim.type,'RN')
        
        tempstr = load(stimfullfname{1});
        if isfield(tempstr,'stimulus')
            if isobject(tempstr.stimulus)
                curpwd = pwd;
                setlocalpath = uigetdir('/Volumes/','Select stimulus class directory');
                cd(setlocalpath);
                tempstr = load(stimfullfname{1});
            end
        else
            warning('''stimulus'' field does not exist in stim file. Cannot proceed !\n');
            return;
        end
        
        time_stamps = tempstr.time_stamps;
        
        % Add the field values for time_stamps into stim structure
        pfields = fieldnames(time_stamps);
        for nf=1:length(pfields)
            stim.time_stamps.(pfields{nf}) = time_stamps.(pfields{nf});
        end
        
        % Add all fields to stim
        if isobject(tempstr.stimulus)
            local_props = properties(tempstr.stimulus);
            for p=1:numel(local_props)
                stim.(local_props{p}) = tempstr.stimulus.(local_props{p});
            end
            cd(curpwd);
        else
            pfields = fieldnames(tempstr.stimulus);
            for nf=1:length(pfields)
                stim.stimulus.(pfields{nf}) = tempstr.stimulus.(pfields{nf});
            end
        end
        
    elseif strcmpi(stim.type,'FP')
    end
    
end
    
    