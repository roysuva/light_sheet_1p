function ind_del = viewNeurons(obj, ind, C2, folder_nm)
%% view all components and delete components manually. it shows spatial
%   components in the full-frame and zoomed-in view. It also shows temporal
%   components
%% input:
%   ind: vector, indices of components to be displayed, no bigger than the maximum
%       number of neurons
%   C2:  K*T matrix, another temporal component to be displayed together
%       with the esitmated C. usually it is C without deconvolution.

%% Author: Pengcheng Zhou, Carnegie Mellon University, 2016

if ~exist('ind', 'var') || isempty(ind)
    % display all neurons if ind is not specified.
    ind = 1:size(obj.A, 2);
elseif ind==-1
    ind = size(obj.A,2):-1:1;
end
if ~exist('C2', 'var'); C2=[]; end

if exist('folder_nm', 'var')&&(~isempty(folder_nm))
    % create a folder to save images
    save_img = true;
    cur_cd = cd();
    if ~exist(folder_nm, 'dir'); mkdir(folder_nm);
    else
        fprintf('The folder has been created and old results will be overwritten. \n');
    end
    cd(folder_nm);
else
    save_img = false;
end

% obj.delete(sum(obj.A>0, 1)<max(obj.options.min_pixel, 1));

Amask = (obj.A~=0);
ind_trim = false(size(ind));    % indicator of trimming neurons
ind_del = false(size(ind));     % indicator of deleting neurons
ctr = obj.estCenter();      %neuron's center
gSiz = obj.options.gSiz;        % maximum size of a neuron

% time
T = size(obj.C, 2);
t = 1:T;
if ~isnan(obj.Fs)
    t = t/obj.Fs;
    str_xlabel = 'Time (Sec.)';
else
    str_xlabel = 'Frame';
end

%% keep the log information
if ~save_img
    try
        log_file =  obj.P.log_file;
        flog = fopen(log_file, 'a');
        log_data = matfile(obj.P.log_data, 'Writable', true); %#ok<NASGU>
        manual_intervention.before = obj.obj2struct();
        
        fprintf(flog, '[%s]\b', get_minute());
        fprintf(flog, 'Start manual interventions:\n');
    end
end

%% reconstruct image for viewing 

% set max array time bins to 6000 samples 
tlim = min([6000 size(obj.C_raw,2)]); 
reconstructed_im = reshape(obj.A*obj.C_raw(:,1:tlim), obj.options.d1, obj.options.d2, tlim);
reconstructed_im = std(reconstructed_im,[],3) + obj.b0_new; % median is costly 
im_rng = [prctile(reconstructed_im(:),0.5) prctile(reconstructed_im(:),99.5)];  

%% extract outliers with delta function like response -false positives 

outlier_idx = extract_delta(obj, ind); 


%% start viewing neurons

% look for preexisting file 
llc = find(obj.file=='/',1,'last'); 
fpath = obj.file(1:llc-1);
if exist(fullfile(fpath,'temp_manual_selection.mat'),'file')
    answer = questdlg('Do you want to start from where you left off?',...
        'Manual sorting',...
        'Yes','No','No');
    switch answer 
        case 'Yes'
            load(fullfile(fpath,'temp_manual_selection.mat')); 
        case 'No'
            m=1; 
            cntr=1; 
            delete(fullfile(fpath,'temp_manual_selection.mat'));
    end
else 
    m=1; 
    cntr=1; 
end


figure('position', [100, 100, 1024, 512]);
while and(m>=1, m<=length(ind))
    
    % save location to pick up later 
    save(fullfile(fpath,'temp_manual_selection.mat'),'m','ind_del','cntr');
        
    %if isempty(outlier_idx) || (outlier_idx(cntr)==ind(m) && cntr<numel(outlier_idx))
    if ~isempty(outlier_idx) && outlier_idx(cntr)==ind(m) && cntr<numel(outlier_idx)
        ind_del(m) = true;
        m = m+1;
        cntr=cntr+1; 
         
    else
        
        %% center of ROI
        x0 = ctr(ind(m), 2);
        y0 = ctr(ind(m), 1);
        aa = Amask(:, ind(m));
        aa = reshape(full(aa),obj.options.d1,obj.options.d2,[]);
        
        %% reconstruction of original image 
        ax1=subplot(2,3,1); cla; 
        obj.image(reconstructed_im, im_rng); 
        colormap(ax1,gray); 
        axis(ax1,'equal','off'); 
        if ~isnan(x0)
            hold on; 
            xotl = x0 + [-gSiz gSiz gSiz -gSiz -gSiz]*4;
            yotl = y0 + [-gSiz -gSiz gSiz gSiz -gSiz]*4;
            line(ax1,xotl,yotl,'Color','r','linewidth',1); 
        end
        %% full-frame view
        ax2=subplot(2,3,2); cla;
        axis equal;
        obj.image(reconstructed_im, im_rng);
        colormap(ax2,gray);
        %obj.image(obj.A(:, ind(m)).*Amask(:, ind(m))); %
        axis(ax2,'equal','off');
        if ~isnan(x0)
            xlim(x0+[-gSiz, gSiz]*4);
            ylim(y0+[-gSiz, gSiz]*4);
        end
        hold on;
        contour(ax2,reshape(full(Amask(:, ind(m))),obj.options.d1,obj.options.d2,[]),'y'); 
        line(ax2,xotl,yotl,'Color','r','linewidth',1); 
        if ind_del(m)
            title(sprintf('Neuron %d', ind(m)), 'color', 'r');
        else
            title(sprintf('Neuron %d', ind(m)));
        end
        %linkaxes([ax1 ax2]); 
        %% zoomed-in view
        subplot(2,3,3); cla;
        obj.image(obj.A(:, ind(m)).*Amask(:, ind(m))); %
        %     imagesc(reshape(obj.A(:, ind(m)).*Amask(:,ind(m))), obj.options.d1, obj.options.d2));
        axis equal; axis off;
        x0 = ctr(ind(m), 2);
        y0 = ctr(ind(m), 1);
        if ~isnan(x0)
            xlim(x0+[-gSiz, gSiz]*2);
            ylim(y0+[-gSiz, gSiz]*2);
        end

        %% temporal components
        subplot(2,3,4:6);cla;
        if ~isempty(C2)
            plot(t, C2(ind(m), :)*max(obj.A(:, ind(m))), 'linewidth', 2); hold on;
            plot(t, obj.C(ind(m), :)*max(obj.A(:, ind(m))), 'r');
        else

            plot(t, obj.C(ind(m), :)*max(obj.A(:, ind(m))));
        end
        xlim([t(1), t(end)]);
        xlabel(str_xlabel);
         

        %% save images
        if save_img
            drawnow();
            saveas(gcf, sprintf('neuron_%d.png', ind(m)));
            m = m+1;
        else
            fprintf('Neuron %d/%d, keep(k, default)/delete(d)/split(s)/trim(t)\n\t/trim cancel(tc)/delete all(da)/backward(b)/end(e):    ', ind(m),length(ind));

            temp = input('', 's');
            if temp=='d'
                ind_del(m) = true;
                m = m+1;
            elseif strcmpi(temp, 'b')
                m = m-1;
            elseif strcmpi(temp, 'da')
                ind_del(m:end) = true;
                break;
            elseif strcmpi(temp, 'k')
                ind_del(m) = false;
                m= m+1;
            elseif strcmpi(temp, 's')
                try
                    subplot(2,3,3);
                    temp = imfreehand();
                    tmp_ind = temp.createMask();
                    tmpA = obj.A(:, ind(m));
                    obj.A(:, end+1) = tmpA.*tmp_ind(:);
                    obj.C(end+1, :) = obj.C(ind(m), :);
                    obj.A(:, ind(m)) = tmpA.*(1-tmp_ind(:));
                    obj.S(end+1, :) = obj.S(ind(m), :);
                    obj.C_raw(end+1, :) = obj.C_raw(ind(m), :);
                    obj.P.kernel_pars(end+1, :) = obj.P.kernel_pars(ind(m), :);
                    k_ids = obj.P.k_ids;
                    obj.ids(end+1) = k_ids+1;   % assign an neuron id
                    obj.tags(end+1) = obj.tags(ind(m));
                    obj.P.k_ids = k_ids+1;
                    fprintf(flog, '\tSplit %d --> %d + %d\n', obj.ids(ind(m)),obj.ids(ind(m)), k_ids);
                catch
                    fprintf('the neuron was not split\n');
                end
            elseif strcmpi(temp, 't')
                try
                    subplot(2,3,3);
                    title('Freehand select region you would like to keep'); 
                    temp = imfreehand();
                    tmp_ind = temp.createMask();
                    Amask(:, ind(m)) = tmp_ind(:);
                    ind_trim(m) = true;
                catch
                    fprintf('the neuron was not trimmed\n');
                end
            elseif strcmpi(temp, 'tc')
                Amask(:, ind(m)) = (obj.A(:, ind(m)) > 0);
                ind_trim(m) = false;
            elseif strcmpi(temp, 'e')
                break;
            elseif ~isnan(str2double(temp))
                m = m + floor(str2double(temp));
                m = max(m, 1);
                m = min(m, length(ind));
                fprintf('jump to neuron %d / %d\n', m, length(ind));
            else
                m = m+1;
            end
            
            
        end
    end   
end

if exist(fullfile(fpath,'temp_manual_selection.mat'),'file')
    delete(fullfile(fpath,'temp_manual_selection.mat')); 
end

if save_img
    cd(cur_cd);
else
    if ~isempty(ind(ind_trim))
        obj.A(:, ind(ind_trim)) = obj.A(:,ind(ind_trim)).*Amask(:, ind(ind_trim));
        try
            fprintf(flog, '\n\tFollowing neurons were trimmed:\n');
            ids_trimmed = ind(ind_trim);
            for m=1:length(ids_trimmed)
                fprintf(flog, '%2d, ', ids_trimmed(m));
            end
            fprintf(flog, '\n');
        end
    end
    
    if ~isempty(ind(ind_del))
        try
            fprintf(flog, '\tDeleting manually selected neurons:\n');
        end
        obj.delete(ind(ind_del));
    end
    %     obj.Coor = obj.get_contours(0.9);
    
    
    return;
end
try
    fprintf(flog, '[%s]\b', get_minute());
    fprintf(flog, 'Finished the manual intervention.\n');
    fprintf(flog, '[%s]\b', get_minute());
    if obj.options.save_intermediate
        manual_intervention.after = obj.obj2struct(); %#ok<STRNU>
        tmp_str = get_date();
        tmp_str=strrep(tmp_str, '-', '_');
        eval(sprintf('log_data.manual_%s = manual_intervention;', tmp_str));
        
        fprintf(flog, '\tThe results were saved as intermediate_results.manual%s\n\n', tmp_str);
    end
    fclose(flog);
end
end

