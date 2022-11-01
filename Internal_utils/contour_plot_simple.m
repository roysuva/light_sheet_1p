function [CC, ha_return] = contour_plot_simple(Aor,Cn,options,varargin)


% save and plot the contour traces of the found spatial components against
% a specified background image. The contour can be determined in two ways:
% options.thr_method = 'max': For every component pixel values below
%       options.thr_method are discarded
% options.thr_method = 'nrg': The contour is drawn around the value above
%        which a specified fraction of energy is explained (default 99%)

% INPUTS:
% Aor:              set of spatial components (matrix d x K)
% Cn:               background image (matrix d1 x d2)
% options:          options structure (optional)
% Varargin:    
%       display_numbers: true or false for showing cell indices
%       max_thresh:      threshold for binarizing   
%       line_clr:        line color  
%       line_wid:        line thickness 
%       text_clr:        text color  
%       fill_clr:        patch color    
%       fill_alpha:      patch transparency (1=opaque, 1=full transparent)
%       clrmap:          pseudo colormap for spatial image 
%       ind_show:        true or false
%       which_rois:      which cell indices to show 
%       clust_info:      if there is cluster info 
%       ranked_rois:     rois ranked by order of cluster rank 
%       ha:              graphics handle for axis 
%       hf:              graphics handle for figure 

% OUTPUTS:
% CC:               contour plots coordinates


% Author: Suva Roy
% Date: 04/15/2020

p = inputParser;
addRequired(p,'Aor');
addRequired(p,'Cn');
addRequired(p,'options');
addParameter(p,'display_numbers',true, @islogical);
addParameter(p,'max_thresh',0.25, @isnumeric); 
addParameter(p,'line_clr',[1 0 1], @(x)(isnumeric(x) || ischar(x))); 
addParameter(p,'line_wid',2, @isnumeric); 
addParameter(p,'text_clr',[1 1 0], @isnumeric); 
addParameter(p,'fill_clr',[1 1 1], @(x)(isnumeric(x) || ischar(x))); 
addParameter(p,'fill_alpha',1, @isnumeric);
addParameter(p,'clrmap',parula(size(Aor,2)+5), @isnumeric);
addParameter(p,'ind_show',true, @islogical);
addParameter(p,'which_rois',[],  @(x)(isnumeric(x) || ischar(x)));
addParameter(p,'clust_info',[],  @isnumeric);
addParameter(p,'ranked_rois',[],  @isnumeric);
addParameter(p,'ha',[],  @isgraphics);
addParameter(p,'hf',[],  @isgraphics);
parse(p,Aor,Cn,options,varargin{:});
params = p.Results; 

% Assign local variables   


intrang = get_intensity_range(params.Cn,'range', [0.5 0.995]);
imshow(params.Cn,intrang,'InitialMagnification','fit','Parent',params.ha); hold on; 
colormap(params.clrmap); 
set(params.ha,'xtick',[],'ytick',[]); hold on; 

if ~isnan(params.ind_show) || params.ind_show
    
    d1 = options.options.d1;
    d2 = options.options.d2;
    if strcmpi(params.which_rois, 'all')
        params.which_rois = 1:size(options.C,1); 
    end
    if ~isempty(params.clust_info)
        clust_avail = true;
    else
        clust_avail = false;
        params.fill_clr = params.fill_clr(1,:); 
        params.text_clr = [1 1 0];
    end

     
    if size(params.line_clr,1)>1 && ~isempty(params.ranked_rois)
        lincols = params.line_clr(params.ranked_rois,:); 
    else
        lincols = repmat(params.line_clr(1,:),length(params.which_rois),1); 
    end
    if size(params.text_clr,1)>1 && isfield('params','ranked_rois') 
        textcols = params.text_clr(params.ranked_rois,:); 
    elseif size(params.text_clr,1)>1 && ~isfield('params','ranked_rois')
        textcols = params.text_clr; 
    elseif (size(params.text_clr,1)==1 && isfield('params','ranked_rois')) || (size(params.text_clr,1)==1 && ~isfield('params','ranked_rois')) 
        textcols = repmat(params.text_clr,length(params.which_rois),1); 
    end
    

    for k = 1:length(params.which_rois)
        i = params.which_rois(k);
        
        A_temp = full(reshape(params.Aor(:,i),d1,d2));
        A_temp = medfilt2(A_temp,[3,3]);
        if sum(A_temp(:)==0)
            A_temp = full(reshape(params.Aor(:,i),d1,d2));     %-------------- modified by Suva Roy. 04/15/2020 --------------$=%
        end
        A_temp(A_temp<params.max_thresh*max(A_temp(:))) = 0;
        BW = bwareafilt(A_temp>0,1);
        BW2 = bwboundaries(BW);
        if ~isempty(BW2) 

            for ii = 1:length(BW2)
                BW2{ii} = fliplr(BW2{ii});

                plot(BW2{ii}(:,1),BW2{ii}(:,2),'Color',lincols(k,:),'linewidth',params.line_wid,'Parent',params.ha); hold on;  
                if ~clust_avail
                    if strcmpi(params.fill_clr,'none')
                        local_color = 'r'; 
                        params.fill_alpha = 1; 
                    else 
                        local_color = params.fill_alpha; 
                    end
                    patch(BW2{ii}(:,1),BW2{ii}(:,2),local_color,'FaceColor',params.fill_clr,'FaceAlpha',params.fill_alpha,...
                        'EdgeColor','none','Parent',params.ha); hold on; 
                else
                    if ~isnan(params.clust_info(k))
                        patch(BW2{ii}(:,1),BW2{ii}(:,2),params.fill_clr(params.clust_info(k),:),...
                            'FaceColor',params.fill_clr(params.clust_info(k),:),'FaceAlpha',params.fill_alpha,...
                            'EdgeColor',params.fill_clr(params.clust_info(k),:),'Parent',params.ha);  hold on;  
                    end
                end
            end
            CC{i} = BW2{1}';
            fp = find(BW);
            [ii,jj] = ind2sub([d1,d2],fp);
            CR{i,1} = [ii,jj]';
            CR{i,2} = A_temp(fp)';
        end
        hold on;
    end
    
    cm=com(params.Aor(:,params.which_rois),d1,d2); % center of mass
    if params.display_numbers
        for k = 1:length(params.which_rois)
            i=params.which_rois(k);
            j=params.clust_info(k); 
            if ~clust_avail
                text(round(cm(k,2)),round(cm(k,1)),strtrim(cellstr(num2str(i))),...
                    'color',textcols(j,:),'fontsize',16,'fontweight','normal','Parent',params.ha); hold on; 
            else
                text(round(cm(k,2))+5,round(cm(k,1))+5,strtrim(cellstr(num2str(params.which_rois(k)))),...
                    'color',textcols(j,:),'fontsize',16,'fontweight','normal','Parent',params.ha); hold on; 
            end
        end
    end
end

ha_return = params.ha;

