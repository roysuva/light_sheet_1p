function bw_frames = get_bw_frames(stim)

    curpwd = pwd;
    duke_devo_path = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Shared/sroy/mgl-devo/DLP-devo/';
    if ~exist(duke_devo_path,'dir')
        warning('Duke devo directory does not exist. Exiting! \n');
        return;
    end
    addpath(genpath(duke_devo_path));
    cd(duke_devo_path);


    width = 912; height = 1140; refresh_rate = 60;
    def_params.x_start =  stim.stimulus.x_start;
    def_params.x_end = stim.stimulus.span_width;
    def_params.y_start =  stim.stimulus.y_start;
    def_params.y_end = stim.stimulus.span_height;
    def_params.back_rgb = stim.stimulus.back_rgb;
    def_params.delay_frames = stim.stimulus.delay_frames;
    def_params.tail_frames = stim.stimulus.tail_frames;
    def_params.refresh_rate = stim.stimulus.tail_frames;


    parameters.class = 'RN';
    parameters.back_rgb = stim.stimulus.parameters.back_rgb;
    parameters.rgb = stim.stimulus.parameters.rgb;
    parameters.seed = stim.stimulus.parameters.seed;
    parameters.binary = stim.stimulus.parameters.binary;
    parameters.probability = stim.stimulus.parameters.probability;
    parameters.jitter = stim.stimulus.parameters.jitter;
    parameters.delay_frames = stim.stimulus.parameters.delay_frames; % onset latency in # of frames (60*10 == 10sec)

    %%%%%%%%%%%%% OLED %%%%%%%%%%%%%%
    parameters.x_start = stim.stimulus.parameters.x_start;
    parameters.x_end = stim.stimulus.parameters.x_end;
    parameters.y_start = stim.stimulus.parameters.y_start;
    parameters.y_end = stim.stimulus.parameters.y_end;


    parameters.independent = stim.stimulus.parameters.independent; % (0: rgb values vary together, 1: rgb values vary independently)
    parameters.interval = stim.stimulus.parameters.interval;  % # of frames before the image changes
    parameters.stixel_height = stim.stimulus.parameters.stixel_height; % (76 default for testing)
    parameters.frames = stim.stimulus.parameters.frames;% % (fps x stimulus time) 60x1200 == 20min stimulus

    parameters.stixel_width = stim.stimulus.parameters.stixel_width;
    parameters.field_width = stim.stimulus.parameters.field_width;
    parameters.field_height = stim.stimulus.parameters.field_height;


    parameters_send = {};
    tmp = [fieldnames(parameters) struct2cell(parameters)]';
    parameters_send = [parameters_send tmp(:)'];


    stimulus = Random_Noise(def_params, parameters_send);
    stimulus.rng_init.state = Init_RNG_JavaStyle(stimulus.rng_init.seed);


    bw_frames = zeros(stimulus.field_width, stimulus.field_height, stimulus.frames,'uint8');
    countdown = 1;
    for i=1:stimulus.frames
        if countdown == 1
            countdown = stimulus.refresh;  % reset to the number of frames specified by "interval"
            eval(stimulus.make_frame_script);  % get new WN pattern
        else
            countdown = countdown - 1;    % it wasn't time to get a new frame, so instead we just decrement the count down
        end
        bw_frames(:,:,i) = squeeze(img_frame(1,:,:));
    end
    
    cd(curpwd); 

end
 