function hf = plot_pref_dirs(resp, StimComb , varargin)


if strcmpi(varargin{1},'Ca')
    primaryvar = resp.ds_estimates_Ca;
elseif strcmpi(varargin{1},'spk')
    primaryvar = resp.ds_estimates_spk;
elseif isempty(varargin{1})
    primaryvar = resp.ds_estimates;
end

ncols = 2;
numcells = length(primaryvar.mag); 
 

hf = figure; 
set(hf,'position',[357 386 1178 502],'color','w'); 

for nc=1:ncols
    ha_cent = subplot(1,ncols, nc, polaraxes);
    if nc==1
        al = zeros(numcells,1); cnt = 1; 
        for cc=1:length(resp.non_dsos_ids)
            arrow_length = primaryvar.mag(resp.non_dsos_ids(cc)); 
            arrow_theta = primaryvar.angle(resp.non_dsos_ids(cc));
            polarplot(ha_cent, [arrow_theta arrow_theta],[0 arrow_length],'-k','linewidth',1); hold on;
            al(cnt) = arrow_length;
            cnt = cnt+1; 
        end
        for cc=1:length(resp.ds_ids)
            arrow_length = primaryvar.mag(resp.ds_ids(cc)); 
            arrow_theta = primaryvar.angle(resp.ds_ids(cc));
            polarplot(ha_cent, [arrow_theta arrow_theta],[0 arrow_length],'-r','linewidth',1); hold on;
            al(cnt) = arrow_length;
            cnt = cnt+1;
        end
        ha_cent.FontSize = 10; 
        ha_cent.RLim = [0 max(al)]; 
        legend_curate({'All RGCs', 'DS RGCs'},'markershow',false,...
            'lineshow',false,'text_color',{'k','r'},'box',false,'foa',ha_cent,'font_size',10); 
    else 
        al = zeros(length(resp.ds_ids),1); 
        for cc=1:length(resp.ds_ids)
            arrow_length = primaryvar.mag(resp.ds_ids(cc)); 
            arrow_theta = primaryvar.angle(resp.ds_ids(cc));
            polarplot(ha_cent, [arrow_theta arrow_theta],[0 arrow_length],'-r','linewidth',2); hold on;
            al(cc) = arrow_length;
        end
        ha_cent.FontSize = 10; 
        ha_cent.RLim = [0 max(al)];  
        legend_curate({'DS RGCs'},'markershow',false,...
            'lineshow',false,'text_color','r','box',false,'foa',ha_cent,'font_size',10); 
    end
end
sgtitle('Preferred direction'); 






