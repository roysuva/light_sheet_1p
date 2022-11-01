function stim = get_stim_info(notebookpath, dataname, datadate, varargin)

    % Specify list of optional parameters
    p = inputParser; 
    addParameter(p,'initpath','/Volumes/dusom_fieldlab/All_Staff/lab/', @ischar); 
    addParameter(p,'stim_options',  {'Pulse','Square','Sine','Noise'}, @iscell);
    addParameter(p,'num_images', 2400, @isnumeric);
    addParameter(p,'fps', 30, @isnumeric);
    p.parse(varargin{:});
    prms = p.Results;
    
    fnames = dir(notebookpath);
    for i=1:length(fnames)
        if ~isempty(regexpi(fnames(i).name,datadate))
            indx = i;
        end
    end

    fullfpath = [notebookpath, fnames(indx).name]; 
    filetext = fileread(fullfpath);
    exprS = dataname(1:end-5);
    exprE = strcat(exprS(1:end-1),num2str(str2double(exprS(end))+1)); 
    filetextmodf = regexprep(filetext,' ','_'); 
    [startindxS,~] = regexp(filetextmodf,exprS);
    [startindxE,~] = regexp(filetextmodf,exprE);
    if ~isempty(startindxE)
        phrase = (filetextmodf(startindxS(1):startindxE(1)-1));
    else
        phrase = filetextmodf(startindxS(1):end);
    end
    phrase = regexprep(phrase,',',''); 
    bolpat = (regexpi(phrase,prms.stim_options)); 
    stimlc = find(cellfun(@(x) ~isempty(x), bolpat)); 
    waveform = prms.stim_options{stimlc}; 

    %-----------

    
    %-----------
    t_axis = (1:prms.num_images).*1/prms.fps; % ms
    
    if strcmpi(waveform,'pulse')
        
        [startindx,endindx] = regexpi(phrase, 'pulse_width');
        if sum(isstrprop(phrase(startindx-7:startindx),'digit'))==0 && sum(isstrprop(phrase(endindx:endindx+7),'digit'))~=0
            phrasechopped = phrase(endindx:endindx+7);
        elseif sum(isstrprop(phrase(startindx-7:startindx),'digit'))~=0 && sum(isstrprop(phrase(endindx:endindx+7),'digit'))==0 
            phrasechopped = phrase(startindx-7:startindx);
        else
            if isempty(sscanf(phrase(startindx-2),'%d')) && ~isempty(sscanf(phrase(endindx+2),'%d'))
                phrasechopped = phrase(endindx:endindx+7);
            else   
                phrasechopped = phrase(startindx-7:startindx);
            end
        end   
        pulse_width = str2double(phrasechopped(isstrprop(phrasechopped,'digit'))); 
        if pulse_width>2 % pulse widths are typically > 2ms, and smaller than 2s
            pulse_width = pulse_width/1000; % sec
        end
        
        [startindx,~] = regexpi(phrase, 'period'); 
        phrasechopped = phrase(startindx+1:end); 
        indentloc = find(phrasechopped=='_',2,'first'); 
        period = sscanf(phrasechopped(indentloc(1)+1:indentloc(2)-1),'%g',1); 
        if ~isempty(regexpi(phrasechopped(indentloc(1)+1:indentloc(2)-1),'ms'))
            period = period/1000; 
        end
        
        [startindx,~] = regexpi(phrase, 'Amp'); 
        phrasechopped = phrase(startindx+1:end); 
        target = [regexpi(phrasechopped,'mV') regexpi(phrasechopped,'V')];
        val = regexpi(phrasechopped(1:target(1)+1),'\d*','Match');
        amplitude = str2double(val{1});
        
        
        numcycles = floor(t_axis(end)/period);
        ytrace = repmat([0 1 1 0],1,numcycles).*amplitude;
        xtrace = [];
        for nc=1:numcycles
            xtrace = [xtrace [0 0 pulse_width pulse_width]+(nc-1)*period];
        end
        
    elseif strcmpi(waveform,'square')
        
        pulse_width = []; 
        
        [startindx,~] = regexpi(phrase, 'period'); 
        phrasechopped = phrase(startindx+1:end); 
        indentloc = find(phrasechopped=='_',2,'first'); 
        period = sscanf(phrasechopped(indentloc(1)+1:indentloc(2)-1),'%g',1); 
        if ~isempty(regexpi(phrasechopped(indentloc(1)+1:indentloc(2)-1),'ms'))
            period = period/1000; 
        end
        
        [startindx,~] = regexpi(phrase, 'Amp'); 
        phrasechopped = phrase(startindx+1:end); 
        target = [regexpi(phrasechopped,'mV') regexpi(phrasechopped,'V')];
        val = regexpi(phrasechopped(1:target(1)+1),'\d*','Match');
        amplitude = str2double(val{1});
        
        numcycles = floor(t_axis(end)/period);
        ytrace = repmat([0 1 1 0],1,numcycles).*amplitude;
        xtrace = [];
        for nc=1:numcycles
            xtrace = [xtrace [0 0 period/2 period/2]+(nc-1)*period];
        end
        
    elseif strcmpi(waveform,'sine')
        
        pulse_width = []; 
        
        [startindx,~] = regexpi(phrase, 'period'); 
        phrasechopped = phrase(startindx+1:end); 
        indentloc = find(phrasechopped=='_',2,'first'); 
        period = sscanf(phrasechopped(indentloc(1)+1:indentloc(2)-1),'%g',1); 
        if ~isempty(regexpi(phrasechopped(indentloc(1)+1:indentloc(2)-1),'ms'))
            period = period/1000; 
        end
        
        [startindx,~] = regexpi(phrase, 'Amp'); 
        phrasechopped = phrase(startindx+1:end); 
        target = [regexpi(phrasechopped,'mV') regexpi(phrasechopped,'V')];
        val = regexpi(phrasechopped(1:target(1)+1),'\d*','Match');
        amplitude = str2double(val{1}); 
        
        
        numcycles = floor(t_axis(end)/period);
        ytrace = repmat(sin(2.*pi.*(0:0.01:period-0.01)./period).*amplitude,1,numcycles);
        xtrace = linspace(0,numcycles.*period, length(ytrace));
        
    end
    
   
    
    % Generate stimulus info 
    stim.waveform = waveform; 
    stim.period = period; % sec 
    stim.amplitude = amplitude; % sec 
    stim.pulse_width = pulse_width; % sec 
    stim.xtrace = xtrace; 
    stim.ytrace = ytrace; 

    % Throw
    fprintf('Retrieved stimulus info\n');
    
end