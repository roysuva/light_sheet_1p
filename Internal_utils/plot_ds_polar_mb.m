function plot_ds_polar_mb(stim, StimComb, resp, savepath,  print_out ,varargin)

% generate subplot indices for DS raster plot with multiple dataruns.
% So far assume that there are 8 directions.
% DS: DS structure with all datasets.
% raster: rasters of all datasets.
% print: 1:print as pdf and close
%        0:don't print

% Author: Suva Roy (2020)) 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 12 Directions : 4x4 subplot
%     Degree              Plot location        Plot index
%  120  90  60  30          5  4  3  2          1  2  3  4
%  150           0   ==     6        1    ==    5  6  7  8
%  180         330          7       12          9 10 11 12
%  210 240 270 300          8  9 10 11         13 14 15 16
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 8 Directions : 3x3 subplot
%     Degree             Plot location        Plot index
%  135  90   45          4  3  2          1  2  3
%  180        0   ==     5     1    ==    4  5  6
%  225  270 305          6  7  8          7  8  9
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmpi(varargin{1},'Ca')
    primaryvar = resp.ds_estimates_Ca;
elseif strcmpi(varargin{1},'spk')
    primaryvar = resp.ds_estimates_spk;
elseif isempty(varargin{1})
    primaryvar = resp.ds_estimates;
end


D = [find(strcmp(fieldnames(stim.combinations),'direction')) find(strcmp(fieldnames(stim.combinations),'DIRECTION'))];
unique_dirs = unique(StimComb(:,D));
num_dirs = length(unique_dirs);
cell_ids = resp.ds_ids;



% check for existence of file
ixx = find(savepath=='/',1,'last');
datname = savepath(ixx+1:end);
name{1} = strcat(savepath,datname,'MB_DSpolar');
name{2} = strcat(savepath,datname,'MB_DSpolar_shaded');

if print_out
    for in=1:length(name)
        if exist(strcat(name{in},'.pdf'),'file')
            s = input('File exists!  Delete (1) or Append to existing (2) ?   : ','s');
            if strcmpi(s,'1')
                delete(strcat(name{in},'.pdf'));
                appendname = '';
            elseif strcmpi(s,'2')
                appendname = '';
            else
                fprintf('Input unrecognized. Creating new file to append...');
                appendname = replace(replace(datestr(now),' ','-'),':','-');
            end
        else
            appendname = '';
        end
        name{in} = strcat(name{in},appendname);
    end
    fprintf('Files are being saved...\n');
else 
    fprintf('Files are not being saved...\n');
end



switch num_dirs
    case 8
        nR = 4; nC = 4;
        subplt_lc = [6 3 2 1 4 7 8 9]';
    case 12
        nR = 4; nC = 4;
        subplt_lc = [8 4 3 2 1 5 9 13 14 15 16 12 6 7 10 11];
        
        
        hf200 = figure(200); set(hf200,'position',[447   254   873   699],'color','w'); 
        for cc=1:length(cell_ids)
            clf(hf200); set(hf200, 'color','w');
            for dr=1:num_dirs
                subplot(nR, nC, subplt_lc(dr));
                
                for ntr=1:stim.repetitions_corrected
                    localtrace = resp.unique_speed(1).unique_width(1).unique_contrast(1).unique_dirs(dr).cell(cell_ids(cc)).raster_unfilt_traces(ntr,:);
                    xx = resp.t_axis_trial;
                    plot(xx, localtrace, '-k'); hold on;
                end
                if dr == 8
                    ylabel('Detrended Calcium signal'); 
                    xlabel('Time (sec)');
                    set(gca,'xtick',[0 ceil(stim.trial_dur)],'xlim',[0 ceil(stim.trial_dur)]); 
                else
                    set(gca,'xtick',[],'ytick',[],'xlim',[0 ceil(stim.trial_dur)]);
                    
                end
                set(gca,'fontsize',12,'box','off');
            end
            
            tt = [unique_dirs; 360].*pi./180;
            tt_rho = primaryvar.rho(cell_ids(cc),:);
            tt_rho = [tt_rho tt_rho(1)];
            arrow_ang = primaryvar.angle(cell_ids(cc));
            arrow_length = max(tt_rho);
            
            ha_cent = subplot(nR, nC, subplt_lc(num_dirs+1:end), polaraxes);
            polarplot(ha_cent,tt, tt_rho,'-','color','b','linewidth',2); hold on;
            polarplot(ha_cent,[arrow_ang arrow_ang],[0 arrow_length],'-b','linewidth',2);
            ha_cent.FontSize = 12;
            
            dsi = primaryvar.ds_index(cell_ids(cc));
            sgtitle(sprintf('%s stim. DSI = %1.2f. Cell id %d',stim.type, dsi, cell_ids(cc)));
            
            if print_out
                export_fig(name{1},'-pdf','-nocrop','-append',hf200);
                %export_fig([pathforfig,figname],'-transparent','-append',gcf);
            else
                pause(1);
            end
        end
        
        hf201 = figure(201); set(hf201,'position',[447   254   873   699],'color','w'); 
        for cc=1:length(cell_ids)
            clf(hf201); set(hf201, 'color','w');
            yllcollect = zeros(num_dirs,2); 
            hacollect = []; 
            for dr=1:num_dirs
                hacollect(dr) = subplot(nR, nC, subplt_lc(dr));
                
                xx = resp.t_axis_trial;
                localtrace = [];
                for ntr=1:stim.repetitions_corrected
                    localtrace = [localtrace; ...
                        resp.unique_speed(1).unique_width(1).unique_contrast(1).unique_dirs(dr).cell(cell_ids(cc)).raster_unfilt_traces(ntr,:)];
                end
                tracestd = std(localtrace,[],1)./2; 
                tracemean = mean(localtrace,1); 
                fillx = [xx xx(end:-1:1)];
                filly = [tracemean+tracestd tracemean(end:-1:1)-tracestd(end:-1:1)];
                plot(xx, tracemean, '-k','linewidth',2); hold on;
                h1 = fill(fillx,filly,'k'); 
                set(h1,'facealpha',0.3,'edgecolor','none'); 
                
                if dr == 8
                    ylabel('Detrended Calcium signal'); 
                    xlabel('Time (sec)');
                    set(gca,'xtick',[0 ceil(stim.trial_dur)],'xlim',[0 ceil(stim.trial_dur)]); 
                else
                    set(gca,'xtick',[],'ytick',[],'xlim',[0 ceil(stim.trial_dur)]);
                end
                set(gca,'fontsize',12,'box','off');
                yllcollect(dr,:) = get(gca,'ylim'); 
            end
            for dr=1:num_dirs
                set(hacollect(dr),'ylim',[0 max(yllcollect(:,2))]); 
            end
            tt = [unique_dirs; 360].*pi./180;
            tt_rho = primaryvar.rho(cell_ids(cc),:);
            tt_rho = [tt_rho tt_rho(1)];
            arrow_ang = primaryvar.angle(cell_ids(cc));
            arrow_length = max(tt_rho);
            
            ha_cent = subplot(nR, nC, subplt_lc(num_dirs+1:end), polaraxes);
            polarplot(ha_cent,tt, tt_rho,'-','color','b','linewidth',2); hold on;
            polarplot(ha_cent,[arrow_ang arrow_ang],[0 arrow_length],'-b','linewidth',2);
            ha_cent.FontSize = 15;
            
            dsi = resp.ds_estimates.ds_index(cell_ids(cc));
            sgtitle(sprintf('%s stim. DSI = %1.2f. Cell id %d',stim.type, dsi, cell_ids(cc)));
            
            if print_out
                export_fig(name{2},'-pdf','-nocrop','-append',hf201);
                %export_fig([pathforfig,figname],'-transparent','-append',gcf);
            else
                pause(1);
            end
        end
        
end


