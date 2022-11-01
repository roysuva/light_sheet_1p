function [stim2framecollect, stimcollect] = get_stim2frame_association_cat(stiminp, fps, varargin)

 

    p = inputParser;
    addParameter(p,'fps_gather',[], @isnumeric);
    addParameter(p,'upsamp_factor',1, @isnumeric); 
    addParameter(p,'nframes_gather',[], @isnumeric);
    parse(p,varargin{:});
    params = p.Results; 

    [stimcollect,stim2framecollect] = deal(cell(1,length(stiminp))); 
    cum_offset = [0; cumsum(params.nframes_gather(1:end-1))];
    fprintf(sprintf('%s Initiating frame adjustment %s\n\n',repelem('.',30),repelem('.',30))); 
    for dn=1:length(stiminp)
        stim = stiminp{dn}; 
        
        [stimcollect{dn}, stim2framecollect{dn}] = ...
            gen_association(stim, cum_offset(dn),params.fps_gather(dn),params.upsamp_factor(dn));
        fprintf(sprintf('%s Frames adjusted for %s stimulus %s\n',repelem('-',30),stim.type,repelem('-',30))); 
    end
    
    function [stim, stim2frame] = gen_association(stim, offset_fr,fps_gather,upsamp_f)
    
        if strcmpi(stim.type,'MB') || strcmpi(stim.type,'MG')

            % Check for missing triggers from camera (threshold = 1/2 of the
            % inter-frame interval)
            if abs(max(diff(stim.time_stamps.external_trig_times))-(1/fps_gather))>0.5*(1/fps_gather)
                lcc = find(diff(stim.time_stamps.external_trig_times)>1.5*(1/fps_gather));
                fprintf('\nA few frame triggers are missing. Adjusting frame numbers.\n');
                if ~isempty(lcc) 
                    lcc = [0; lcc]; 
                end

                stim.time_stamps.external_trig_times_corrected = [];
                for lc=1:length(lcc)-1
                    stim.time_stamps.external_trig_times_corrected = [stim.time_stamps.external_trig_times_corrected; ...
                        stim.time_stamps.external_trig_times(lcc(lc)+1:lcc(lc+1))]; 
                    stim.time_stamps.external_trig_times_corrected = [stim.time_stamps.external_trig_times_corrected; ...                    
                        diff(stim.time_stamps.external_trig_times(lcc(lc+1):lcc(lc+1)+1))/2 + stim.time_stamps.external_trig_times(lcc(lc+1))]; 
                end
                stim.time_stamps.external_trig_times_corrected = [stim.time_stamps.external_trig_times_corrected; ...
                    stim.time_stamps.external_trig_times(lcc(lc+1)+1:end)];
            end


            % Check if initial offset between frame acquisition and stimulus delivery
            % is equivalent to more than 1 frame
            if stim.time_stamps.stim_time_stamps(2,2) > (1/fps_gather)            
                fprintf('\nInital offset requires correction.\n');
            else
                fprintf('\nInital offset doesn''t require correction.\n');
            end

            % Assign frame range for all stimulus combinations 
            frame2stim = zeros(length(stim.trials),2); 
            %frametimes = stim.time_stamps.external_trig_times_corrected - stim.time_stamps.external_trig_times_corrected(1); 
            fps_upsamp = fps_gather*upsamp_f; 
            for ntr=1:length(stim.trials)
                frame2stim(ntr,:) = [ceil(stim.time_stamps.stim_time_stamps(ntr+1,2).*fps_upsamp) ...
                    fix(stim.time_stamps.stim_time_stamps(ntr+1,3).*fps_upsamp)]; 
                if ntr>1 
                    if frame2stim(ntr-1,2)==frame2stim(ntr,1)
                        frame2stim(ntr,1)= frame2stim(ntr,1)+1;
                    end
                end
            end
            if frame2stim(1,1)==0
                frame2stim(1,1) = 1;
            end

            inp = input('Are all trials equal length ? (y/n)   ','s'); 
            if strcmpi(inp, 'y')
                minnumfr = min(diff(frame2stim,[],2)+1); % equalize trials to this many frames 
                frame2stim(:,2) = frame2stim(:,1) + minnumfr-1; 
            end

            if frame2stim(end,2)>length(stim.time_stamps.external_trig_times_corrected)*upsamp_f
                fprintf('\nStimulus run time longer than image acquisition time. Missing frames at the end.\n');
            elseif frame2stim(end,2)<length(stim.time_stamps.external_trig_times_corrected)*upsamp_f
                fprintf('\nStimulus run time shorter than image acquisition time. Extra frames at the end.\n');
            else
                fprintf('\nStimulus run time matches image acquision time.\n');
            end

            stim2frame.trial = [1:length(stim.trials)]';
            stim2frame.frame_rng = frame2stim+offset_fr;


        elseif strcmpi(stim.type,'CH')



            % Check for missing triggers from camera (threshold = 1/2 of the
            % inter-frame interval)
            if abs(max(diff(stim.time_stamps.external_trig_times))-(1/fps_gather))>0.5*(1/fps_gather)
                lcc = find(diff(stim.time_stamps.external_trig_times)>1.5*(1/fps_gather));
                fprintf(sprintf('\n%d frames missing.\n',length(lcc)));
            else
                fprintf('\nNo missing frames.\n');
            end


            % Check if initial offset between frame acquisition and stimulus delivery
            % is equivalent to more than 1 frame
            if abs(stim.time_stamps.stim_time_stamps(1,1) - stim.time_stamps.external_trig_times(1)) > (1/fps_gather)
                fprintf('Inital offset requires correction.\n');
            else
                fprintf('Inital offset doesn''t require correction.\n');
            end

            numphase = max(stim.stimulus.phase);

            state = stim.stimulus.transition_states'; 
            trial  = stim.stimulus.transition_trial_num';
            phase = stim.stimulus.phase'; 
            frame_rng_dlp = zeros(length(stim.stimulus.transition_states),2);
            totalnumcycles = length(stim.stimulus.transition_states); 
            for tnc=1:totalnumcycles
                if tnc==1
                    frame_rng_dlp(tnc,1:2) = [1 stim.stimulus.transition_frames(tnc)-1];
                else
                    frame_rng_dlp(tnc,1:2) = [stim.stimulus.transition_frames(tnc-1) stim.stimulus.transition_frames(tnc)-1]; 
                end
            end
            mould = stim.stimulus.phase(1:length(stim.stimulus.phase)/stim.stimulus.num_trials); 
            indx = cell(1,numphase); 
            for i=1:numphase 
                indx{i} = find(mould(1,:)==i); 
            end
            for i=1:numphase 
                mould(indx{i}) = 1:length(indx{i}); 
            end
            cycle = repmat(mould, 1, stim.stimulus.num_trials)'; 

            frame_rng_cam = zeros(length(stim.stimulus.transition_states),2);
            fps_local = fps_gather*upsamp_f; 
            for fr=1:size(frame_rng_cam,1)            
                frame_rng_cam(fr,:) = fix(stim.time_stamps.stim_time_stamps(fr+1,2:3).*fps_local); 
            end
            frame_rng_cam(:,1) = frame_rng_cam(:,1)+1; 
            if frame_rng_cam(end,2)>length(stim.time_stamps.external_trig_times)*upsamp_f
                frame_rng_cam(end,2) = length(stim.time_stamps.external_trig_times)*upsamp_f; 
            end
            
            frame_rng_cam = frame_rng_cam + offset_fr;

            stim2frame = table(frame_rng_dlp,frame_rng_cam,trial,state,phase,cycle, ...
                'VariableNames',{'frame_rng_dlp','frame_rng_cam','trial','state','phase','cycle'});
            
            


        elseif strcmpi(stim.type,'RN') 

            % Determine offset between stim triggers and camera triggers 
            initoffset = abs(stim.time_stamps.external_trig_times(1)-stim.time_stamps.stim_time_stamps(1));
            if stim.time_stamps.external_trig_times(1)>stim.time_stamps.stim_time_stamps(1)
                leadname = 'stim';
            else
                leadname = 'cam';
            end
            fprintf(sprintf('Inital offset = %1.4f s.\n',initoffset));
            
            % Display refresh rate 
            if isfield(stim.stimulus,'refresh_rate')
                display_ifi = 1/stim.stimulus.refresh_rate; 
            else
                display_ifi = 1/60; 
                stim.stimulus.refresh_rate = 60; 
            end

            % Check for missing triggers from camera 
            ideal_ifi = median(diff(stim.time_stamps.external_trig_times)); 
            lcc = find(diff(stim.time_stamps.external_trig_times)>1.5*ideal_ifi);
            lcc = [lcc; length(stim.time_stamps.external_trig_times)]; 
            stim2frame.cam_frame_num_corrected = [1:lcc(1)]';
            stim2frame.cam_frame_time_corrected = stim.time_stamps.external_trig_times(1:lcc(1));
            for lc=1:length(lcc)-1
                insert = round(diff(stim.time_stamps.external_trig_times(lcc(lc):lcc(lc)+1))/(1/fps_gather))-1; 
                stim2frame.cam_frame_num_corrected = [stim2frame.cam_frame_num_corrected; length(stim2frame.cam_frame_num_corrected)+(1:insert)']; 
                stim2frame.cam_frame_time_corrected = [stim2frame.cam_frame_time_corrected; stim2frame.cam_frame_time_corrected(end)+(1:insert)'.*ideal_ifi]; 

                stim2frame.cam_frame_num_corrected = [stim2frame.cam_frame_num_corrected; ...
                    length(stim2frame.cam_frame_num_corrected)+(1:length(lcc(lc)+1:lcc(lc+1)))']; 
                stim2frame.cam_frame_time_corrected = [stim2frame.cam_frame_time_corrected; ...
                    stim.time_stamps.external_trig_times(lcc(lc)+1:lcc(lc+1))]; 
            end
            stim2frame.cam_frame_time_corrected = stim2frame.cam_frame_time_corrected - stim2frame.cam_frame_time_corrected(1); % assuming 'vsync' mode
            stim2frame.cam_frame_num_corrected = stim2frame.cam_frame_num_corrected-1; % matching frame num with frame time

            % Check for missing triggers for stimulus, adjust initial
            % offset between stim and camera trigs, and initial delay in
            % stimulus display  
            stim_time_stamps = stim.time_stamps.stim_time_stamps(2:end); 
            lcc1 = find(diff(stim_time_stamps)>display_ifi*1.5); 
            lcc2 = find(diff(stim_time_stamps)<display_ifi*0.5); 
            rng = min([lcc1; lcc2]):max([lcc1; lcc2]);
            if abs((sum(diff(stim_time_stamps(rng)))-(length(rng)*display_ifi)))/display_ifi < 1
                fprintf(sprintf('No missing stimulus triggers.\n')); 
            end

            if stim.stimulus.delay_frames~=0 
                stim2frame.stim_frame_time_corrected = initoffset+[stim_time_stamps(1)-display_ifi; stim_time_stamps];
            else
                stim2frame.stim_frame_time_corrected = initoffset+[0; stim_time_stamps]; 
            end
            stim2frame.stim_frame_num_corrected = (0:length(stim2frame.stim_frame_time_corrected)-1)'; 
            
            % adjust for upsampling 
            yt = stim2frame.cam_frame_time_corrected; 
            x = 1:upsamp_f:upsamp_f*length(yt);
            x_ = 1:x(end); 
            stim2frame.cam_frame_time_corrected = interp1(x,yt,x_)'; 
            stim2frame.cam_frame_num_corrected = [0:1:length(stim2frame.cam_frame_time_corrected)-1]';
            

            % Align end times 
            if stim2frame.stim_frame_time_corrected(end)>=stim2frame.cam_frame_time_corrected(end)
                final_stim_frame_needed = find(stim2frame.stim_frame_time_corrected>stim2frame.cam_frame_time_corrected(end),1,'first'); 
                stim2frame.stim_frame_time_corrected = stim2frame.stim_frame_time_corrected(1:final_stim_frame_needed); 
                stim2frame.stim_frame_num_corrected = stim2frame.stim_frame_num_corrected(1:final_stim_frame_needed); 
            else 
                final_cam_frame_needed = find(stim2frame.cam_frame_time_corrected<stim2frame.stim_frame_time_corrected(end),1,'last');
                stim2frame.cam_frame_time_corrected = stim2frame.cam_frame_time_corrected(1:final_cam_frame_needed); 
                stim2frame.cam_frame_num_corrected = stim2frame.cam_frame_num_corrected(1:final_cam_frame_needed); 
            end
            
      
            % Adjust for concatenation offset
            stim2frame.cam_frame_time_corrected = stim2frame.cam_frame_time_corrected + (offset_fr/fps);
            stim2frame.cam_frame_num_corrected = stim2frame.cam_frame_num_corrected + offset_fr; 

            % Get other relevant parameter values 
            stim2frame.display_ifi = display_ifi; 
            stim2frame.initoffset_cam_stim = initoffset; 
            stim2frame.who_leads = leadname; 
            
            % 

        end
    end
        
end




