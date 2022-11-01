% Temporary script for generating stimulus, calcium and dF/F0 traces for
% specific ROIs

figure;
cids = [61:64];
for cc=1:length(cids)
    xx = (1:size(neuron.C_raw,2)).*(1/fps);
    yy = neuron.C_raw(cids(cc),:);
    span = ceil(1*fps);
    yy_smooth = smooth(yy,span,'sgolay',2);
    plot(xx,yy_smooth,'-');
    title(sprintf('Comp %d',cids(cc)));
    pause()
end


% Select frames
cids = [8 10 11 13 14 18 21 22 24 25 29 30 31 33]; 
dT = [20 40]; % seconds
frame_rng = max([ceil(dT(1).*fps) 1]):ceil(dT(2).*fps);
d3_temp = numel(frame_rng);
ncycles = floor(d3 * (1/fps) / stim.period);
transition_fr = round(stim.period*(1:ncycles) * fps);
transition_fr = [1 transition_fr(1:end-1) d3];
[~,ia,ib] = intersect(frame_rng, transition_fr);
transition_fr_within = frame_rng(ia); 
ncycles_within = length(transition_fr_within); 
fr_per_cycle = diff(transition_fr_within(1:2));

[trace_smooth,trace,dF_F,dF_F_smooth] = deal(zeros(length(cids),d3_temp));
[sig, base, sig_smooth, base_smooth] = deal(zeros(length(cids),ncycles_within-1));
span = ceil(0.1*fps);
for cc=1:length(cids)
    trace_smooth_temp = smooth(neuron.C_raw(cids(cc),:),span,'sgolay',1); % (1) *C_raw_dedrifted (2) *neuron.C_raw (3) *dF_F0_raw_dedrifted
    if size(trace_smooth_temp,1)>size(trace_smooth_temp,2)
        trace_smooth_temp = trace_smooth_temp';
    end
    trace_smooth(cc,:) = trace_smooth_temp(1,frame_rng);
    trace(cc,:) = neuron.C_raw(cids(cc),frame_rng); 
    for nc=2:ncycles_within
        sig(cc,nc-1) = max(neuron.C_raw(cids(cc),transition_fr_within(nc-1):transition_fr_within(nc)-1));
        %base(cc,nc-1) = mean(neuron.C_raw(cids(cc),transition_fr_within(nc)-3:transition_fr_within(nc)));
        base(cc,nc-1) = min(neuron.C_raw(cids(cc),transition_fr_within(nc-1):transition_fr_within(nc)-1));
        dF_F(cc,(nc-2)*fr_per_cycle+1:(nc-1)*fr_per_cycle) = ...
            (neuron.C_raw(cids(cc),transition_fr_within(nc-1):transition_fr_within(nc)-1)-base(cc,nc-1))./abs(base(cc,nc-1));
        
        sig_smooth(cc,nc-1) = max(trace_smooth_temp(1,transition_fr_within(nc-1):transition_fr_within(nc)-1));
        %base(cc,nc-1) = mean(trace_smooth_temp(1,transition_fr_within(nc)-3:transition_fr_within(nc)));
        base_smooth(cc,nc-1) = min(trace_smooth_temp(1,transition_fr_within(nc-1):transition_fr_within(nc)-1));
        dF_F_smooth(cc,(nc-2)*fr_per_cycle+1:(nc-1)*fr_per_cycle) = ...
            (trace_smooth_temp(1,transition_fr_within(nc-1):transition_fr_within(nc)-1)-base_smooth(cc,nc-1))./abs(base_smooth(cc,nc-1));
    end
end


% Get stim trace 
xx = frame_rng.*(1/fps);
stim_temp = get_stim_info(notebookpath, datanam, datadat, 'num_images', d3_temp, 'fps', fps);

n=7;
figure; set(gcf,'position',[561   617   651   302],'color','w'); 
subplot(3,1,1); 
plot(stim_temp.xtrace, stim_temp.ytrace./max(stim_temp.ytrace),'-r','linewidth',0.5); 
set(gca,'fontsize',15);
subplot(3,1,2); 
plot(xx,trace(n,:),'-k');
set(gca,'fontsize',15);
subplot(3,1,3); 
plot(xx,dF_F(n,:),'-k');
set(gca,'fontsize',15);


figure; set(gcf,'position',[562   727   683   135],'color','w'); 
plot(stim_temp.xtrace, stim_temp.ytrace./max(stim_temp.ytrace),'-r','linewidth',0.5);
set(gca,'fontsize',15,'box', 'off');
n=7;
figure; set(gcf,'position',[562   727   683   135],'color','w');
plot(xx,dF_F(n,:),'-k','linewidth',1);
ylabel('dF/F0'); 
set(gca,'fontsize',15,'box', 'off');

figure; 
for cc=1:length(cids)
    subplot(length(cids),1,cc); 
    plot(neuron.C_raw(cids(cc),:)); 
end



