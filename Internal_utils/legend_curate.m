function legend_curate(labels, varargin)


% This function will set the appropriate color, text color, box positions, 
% location, background etc. in a legend based on user input.  
%
% Required parameter:
%       labels          :   legend names
%
% Optional parameters: 
%       markertype      :   the marker type ('*', 'o' etc.)
%       linetype        :   line type ('-', '--', etc.)
%       markercolor     :   if marker, then marker color 
%       linecolor       :   if line, then line color 
%       markershow      :   logical true/false to show/not show marker  (default: false)
%       lineshow        :   logical true/false to show/not show line  (default: true)
%       text_color      :   text color
%       font_size       :   fontsize (dafault: 9)
%       font_weight     :   'normal' or 'bold' (default: 'normal')
%       box             :   logical true/false for showing/not showing bounding box (default: true)
%       transparency    :   [0 - 1]. 0: opaque. 1: transparent (default: 0)
%       background_color:   background color of box  (default: 'w') 
%       box_pos         :   position of box   (default normalized units)
%       location        :   location of legend  (default: 'NE')
%       foa             :   figure or axis handle  
%       
%
% Author: Suva Roy, Duke University, 2020.



% SET UP OPTIONAL ARGUMENTS
p = inputParser;
p.addParameter('markertype', [], @ischar);
p.addParameter('linetype', '-', @ischar);
p.addParameter('markercolor', [], @(x)(ischar(x)||iscell(x)||isnumeric(x)));
p.addParameter('linecolor', 'k', @(x)(ischar(x)||iscell(x)||isnumeric(x))); 
p.addParameter('markershow', false, @islogical);
p.addParameter('lineshow', true, @islogical);
p.addParameter('text_color', 'k', @(x)(ischar(x)||iscell(x)||isnumeric(x)));
p.addParameter('font_size', 9, @isnumeric);
p.addParameter('font_weight', 'normal', @ischar)
p.addParameter('box', true, @islogical);
p.addParameter('transparency', 0, @isnumeric);
p.addParameter('background_color', [1 1 1], @(x)(ischar(x)||iscell(x)||isnumeric(x)));
p.addParameter('box_pos',[0.7875 0.8702 0.1036 0.0357], @isnumeric); 
p.addParameter('location', 'northeast', @ischar);
p.addParameter('foa', gca); 

p.parse(varargin{:});
prms = p.Results; 


if ~exist('labels','var') 
    error('You need to provide legends!') 
elseif isempty(labels)
    error('Legend is empty!');
end

if strcmpi(class(prms.foa),'matlab.ui.Figure')
    hf = prms.foa;
    ha = gca; 
elseif strcmpi(class(prms.foa),'matlab.graphics.axis.Axes')
    ha = prms.foa; 
    hf = gcf; 
end

% Set legend 
[leg_h,obj_h,plt_h,~] = legend(labels,'color',prms.background_color);

% Set location
set(leg_h,'location',prms.location); 


% Set line/marker properties (line precedes marker) 
disptype = findobj(plt_h,'type','line');
if ~isempty(prms.markertype) % marker takes precedence 
    if ~prms.markershow
        legline = findobj(obj_h,'type','line'); 
        for lg=1:length(legline)
            set(legline(lg),'visible','off'); 
        end
    else
        clrval = order_color(prms.markercolor, labels); 
        for dt=1:length(disptype)
            set(disptype(dt), 'Color', clrval{dt}); 
        end
    end
else                         % line takes precedence 
    if ~prms.lineshow
        legline = findobj(obj_h,'type','line'); 
        for lg=1:length(legline)
            set(legline(lg),'visible','off'); 
        end
    else
        clrval = order_color(prms.linecolor, labels); 
        for dt=1:length(disptype)
            set(disptype(dt), 'Color', clrval{dt}); 
        end
    end
end

    

% Set text color, weight, fontweight, fontsize 
clrval = order_color(prms.text_color, labels); 
legtxt = findobj(obj_h,'type','text'); 
for lt=1:length(legtxt)
    pause(0.000001); 
    set(legtxt(lt),'color', clrval{lt});
    set(legtxt(lt),'fontsize', prms.font_size); 
    set(legtxt(lt),'fontweight', prms.font_weight); 
end



% Set box properties
if ischar(prms.background_color)
    boxclr = search_colors(prms.background_color);
else
    boxclr = prms.background_color; 
end
set(leg_h.BoxFace, 'ColorType','truecoloralpha', 'ColorData',uint8([boxclr.*255 (1-prms.transparency)*255]'));
if ~prms.box
    set(leg_h.BoxEdge,'Visible','off'); % Turn off visibility 
end
%set(leg_h,'Position',prms.box_pos,'Units','normalized');





function clrval = order_color(input_colorarr, labels)
if ischar(input_colorarr) || iscell(input_colorarr)
    clrlen = length(input_colorarr);
    clrval = cell(1,length(labels));
    for u=1:length(labels)
        if clrlen==1
            clrval{u} = input_colorarr;
        else
            clrval{u} = input_colorarr{u};
        end
    end
else
    clrlen = size(input_colorarr,1);
    clrval = cell(1,length(labels));
    for u=1:length(labels)
        if clrlen==1
            clrval{u} = input_colorarr(1,:);
        else
            clrval{u} = input_colorarr(u,:);
        end
    end
end


function retclrval = search_colors(assignment)
charassign = {'k','w','r','g','b','c','m','y'};
rgbassign = [0 0 0; 1 1 1; 1 0 0; 0 1 0; 0 0 1; 0 1 1; 1 0 1; 1 1 0]; 

if iscell(assignment)
    [~,lc] = intersect(charassign, assignment,'stable');
elseif ischar(assignment)
    [~,lc] = intersect(charassign, cellstr(assignment),'stable');
end
retclrval = rgbassign(lc,:); 

    

    

