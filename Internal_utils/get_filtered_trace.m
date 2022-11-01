function [yy_filtered] = get_filtered_trace(yy,fs,varargin)


p = inputParser;
addRequired(p,'yy',@isnumeric);
addRequired(p,'fs',@isnumeric); 
addParameter(p,'filter_type','LowPass', @ischar); 
addParameter(p,'cutoff_dB',-10, @isnumeric);  % default value (a soft cutoff of -6dB in the PSD)
parse(p,yy,fs,varargin{:});
params = p.Results; 

yy = params.yy; 
fs = params.fs; 
%cutoff_dB = params.cutoff_dB+0.1;


switch params.filter_type
    case 'LowPass' 
        
        nyquist_f = fs/2; 
        cutoff_freq = nyquist_f/5; % reasonable choice 
        yy_filtered = lowpass(yy,cutoff_freq,fs); 


%         nyquist_f = fs/2; lc = 1; 
%         [psdx,ff] = pwelch(yy,[],[],[],fs);
%         interp_ff = linspace(ff(1),ff(end),10000);
%         interp_pow = interp1(ff,psdx,interp_ff,'spline'); 
%         while lc<=1 
%             cutoff_dB = cutoff_dB-0.1; % decrement in cutoff frequency to accomodate range of sample frequencies
%             cutoff_freq = 10^(cutoff_dB/10);
%             [~,lc] = min(abs(interp_ff-cutoff_freq));
%         end 
%         yy_filtered = lowpass(yy,cutoff_freq,fs);
        
    case 'FIR' 
end


