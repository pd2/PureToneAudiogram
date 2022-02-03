% This script implements hearing test using train of pure tones
%   and outputs Pure Tone Audiogram (PTA) plots
%
% Author:- Pradeep D, Uhlhaas lab,
% Institute of Neuroscience & Psychology, University of Glasgow, Glasgow, UK
% (C) copyright reserved 2022

% Version History:-
%  ver - DD MMM YYYY - Feature added
%  1.0 - 31 Jan 2022 - first implementation

clear; close all; clc

%% user defined configuration

% specify the path to store the results
out_path = './results/';

%% Configure experiment

% set the list of frequencies to be tested in kHz
list_freq =  [ 0.25 0.5 1 2 4 6 8]; % Standard PTA frequencies

% set the list of dB(HL) hearing levels to systematically test hearing
list_dB_HL = 0 : 5 : 85;

% set flag to enable GUI for user response
flag_enable_gui = 1;


%% default configurations

% set the number of repetitions
num_repetitions = 3;

% set the ear sides that needs to be tested
ear_sides = {'Left', 'Right'};

% set sampling frequency
samp_freq = 48000;

% set max hearing level that can be tested
max_dB_HL = 85;

% ensure that the state of the pause command is set to pause
pause_state = 'on';

% set response time in seconds at end of trial
resp_window_duration = inf; % sec

% set the question to request a response from user through a GUI window
resp_question = 'did you hear these tones?';

% set the title of the GUI window to request a response
resp_title = 'Please respond';

% set the strings to display for both the options
str_choice_1 = 'Yes';
str_choice_2 = 'No';

% set the default choice
str_def_choice = str_choice_2;

% set colours for the audiogram plot
col = {'b' 'r'};

% set the statement to display between runs
inter_ear_statement = 'Press any key to continue testing on other ear';

% set the names for header columns
header_names = ...
    {
    'Ear';
    'Freq';
    'dB SPL'
    'Repeat';
    'Choice';
    'Latency';
    };

% set flag to enable debug code
flag_debug = 0;


%% initialization

flag_0 = 0;
flag_1 = 1;

% obtain participant name
if flag_debug==0
    participant = input('Please enter your name: ', 's');
    if isempty(participant)
        % assign default name upon no response
        participant = 'Subject';
    end
else
    participant = 'testing';
end

% compute the number of factors excluding the test factor
num_factors = 2;

% number of ears for this hearing test
num_ears = length(ear_sides);

% number of frequencies for this hearing test
num_freqs = length(list_freq);

% number of hearing levels for this hearing test
num_dB_HL = length(list_dB_HL);

% number of levels in each factor
list_num_levels = [num_ears num_freqs];

% initialize the counter for threshold for entire experiment
thresholds = zeros(list_num_levels);

list_heard_levels = zeros([list_num_levels num_repetitions]);

% initialize the results tabulation
num_cols = length(header_names);

% results_table = zeros(num_loops, num_cols);
result_cur_trial = zeros(1,num_cols);

kk = 1;

%% file processing

% note time for file name generation
log_time = clock;
time_suffix = sprintf('%s_%.2d_%.2d_%.2d',date,log_time(4),log_time(5),round(log_time(6)));

% create output directory if not existent
if ~(exist(out_path,'dir'))
    mkdir(out_path);
end

% formulate the file names to log and store the results
if flag_debug==0
    log_file_name = [participant '_' time_suffix '_log.txt'];
else
    log_file_name = [participant '_log.txt'];
end
out_file_name = [participant '_results'];

% open file for writing log
hlog_file = fopen([out_path '/' log_file_name],'w');

if flag_debug==0
    out_path_exp = [out_path '/' time_suffix '/' ];
    
    % create output directory if not existent
    if ~(exist(out_path_exp,'dir'))
        mkdir(out_path_exp);
    end
else
    out_path_exp = out_path;
end

% concatenate output file name with path
hout_file = fopen([out_path_exp '/' out_file_name '.txt'],'w');

% results file header
for ii = 1:num_cols
    fprintf(hout_file,'%s \t', header_names{ii});
end
fprintf(hout_file,'\n\n');

%% processing

% set the state of pause funcionality
pause(pause_state);

% iterate over each ear
for ear_indx = 1:num_ears
    
    ear = ear_sides{ear_indx};
    
    % iterate over all frequencies
    for fr = 1:num_freqs
        
        % get current frequency
        cur_freq = list_freq(fr)*1000;
        
        % generate the acoustic stimulus to be tested
        stimulus_orig = generate_tones(cur_freq, ear, samp_freq);
        
        if flag_debug==0

            % load stimulus
            player = audioplayer(stimulus_orig, samp_freq);

            % play the audio
            playblocking(player);

            % inform subject through a popup window
            if flag_enable_gui
                str_resp = questdlg(['In your ' ear ' ear, sample tones were played that you need to detect at low volume'], ...
                    resp_title, 'OK','Cancel','OK');
                
                if strcmp(str_resp, 'Cancel')
                    response = 0;
                else
                    response = 1;
                end
            else
               % inform sbject through command line
                resp = input(['In your ' ear ' ear, sample tones were played that you need to detect at low volume: '],'s');
                
                if (strcmpi(resp,'N') || strcmpi(resp,'No'))
                    response = 0;
                else
                    response = 1;
                end
            end
            
            %% if user does not make a choice, give opportunity to exit
            if response==0
                disp ' ';
                % disp('To continue, Press F5. To exit, press Shift+F5'); keyboard;
                intent = input('Do you wish to quit? Y or [N]: ', 's');
                if isempty(intent)
                    % assign default answer
                    intent = 'n';
                end
                if strcmpi('Y',intent)
                    disp('Quitting..');
                    fclose(hlog_file);
                    fclose(hout_file);
                    return; % exit the program
                else
                    disp('Continuing..');
                    disp ' ';
                end
            end
        end
        
        % reset repeat counter
        repeat = 0;
        
        % heard levels
        heard_levels = repmat(max(list_dB_HL),num_repetitions, 1);
        
        %% start loop counters
        jj = 1; ii = 1;
        
        % iterate over all hearing levels
        while (jj <= num_dB_HL)
            
            % get current hearing level
            cur_dB_HL = list_dB_HL(jj);
            
            if cur_dB_HL > max_dB_HL
                cur_dB_HL = max_dB_HL;
            end
            
            disp(['Ear: ' ear ' Freq: ' num2str(cur_freq) ' Trial: ' num2str(ii)]);
            
            % store parameters of the current trial
            result_cur_trial(1) = ear_indx;
            result_cur_trial(2) = cur_freq;
            result_cur_trial(3) = cur_dB_HL;
            result_cur_trial(4) = repeat;
            
            % log the trial parameters
            fprintf(hlog_file, 'Ear: %s Freq: %d dB HL: %d\n', ear, cur_freq, cur_dB_HL);
            
            %% stimulus presentation
            
            % compute the gain level
            gain = 10^((cur_dB_HL - max_dB_HL)/20);
            
            % compute the stimulus with appropriate gain level
            stimulus = stimulus_orig * gain;
            
            if flag_debug==0
                
                % load stimulus
                player = audioplayer(stimulus, samp_freq);
                
                % play the audio
                playblocking(player);
                
            end
            
            %% obtain user response
            
            % start the counter for measuring response latency
            t = tic;
            
            if flag_debug==0
                
                % obtain response through a popup window
                if flag_enable_gui
                    str_resp = questdlg(['In your ' ear ' ear, ' resp_question], ...
                        ['Trial: ' num2str(ii) ' ' resp_title], str_choice_1,str_choice_2,str_def_choice);
                    
                    % convert the recorded response to value
                    if strcmp(str_resp, str_choice_1)
                        response = 1;
                    elseif strcmp(str_resp, str_choice_2)
                        response = 2;
                    else
                        response = 0; % empty string when user does not respond
                    end
                else
                    resp = input([resp_question '  Y (1) or N (2)?: '],'s');
                    
                    % convert user response to standard variable
                    if (strcmpi(resp,'Y') || strcmpi(resp,'1'))
                        str_resp = str_choice_1;
                        response = 1;
                    elseif (strcmpi(resp,'N') || strcmpi(resp,'2'))
                        str_resp = str_choice_2;
                        response = 2;
                    elseif isempty(resp)
                        str_resp = ''; % empty string when user does not respond
                        response = 0;
                    else
                        str_resp = '?';
                        response = -1;% some response but neither choice nor empty
                    end
                end
                
            else
                % Avoid UI during debug. Assign a random choice as output
                rnd_choice = (rand(1) > (0.1 + repeat/4) ) + 1;
                if rnd_choice(1)==1
                    str_resp = str_choice_1;
                else
                    str_resp = str_choice_2;
                end
                response = rnd_choice(1);
            end
            
            % stop the counter measuring response latency
            resp_lat = toc(t);
            
            % store parameters of the current trial
            result_cur_trial(end-1) = response;
            
            % store parameters of the current trial
            result_cur_trial(end) = resp_lat;
            
            % obtain response through a popup window - with or w/o time limits
            if flag_debug
                disp(['Response: ' str_resp '   Latency: ' num2str(resp_lat,'%.3f') ' s']);
            end
            
            % log the user response
            fprintf(hlog_file,'User choice: %s \nLatency: %.3f s\n',str_resp,resp_lat);
            
            %% record the trial details
            
            % update the overall results table
            results_table(kk,:) = result_cur_trial;
            
            % store intermediate results in a text file
            fprintf(hout_file,'%.4g\t',result_cur_trial);
            
            fprintf(hout_file,'\n');
            fprintf(hlog_file,'\n');
            
            %% Repeat assessment logic
            
            if response == 1
                
                heard_levels(repeat+1) = cur_dB_HL;
                
                % has all repeat assessments complete
                if repeat == (num_repetitions-1)
                    % stop assessment for this frequency
                    break;
                else
                    % perform a repeat trial
                    repeat = repeat + 1;
                    
                    % increase gain by 9 dB HL
                    jj = jj - 3;
                    
                    % ensure index is positive
                    if (jj < 1)
                        jj = 1;
                    end
                end
            elseif response >= 1
                % increment counter
                jj = jj + 1;
            end
            
            % increment counters
            ii = ii + 1;
            kk = kk + 1;
            
            %% if user does not make a choice, give opportunity to exit
            if response==0
                disp ' ';
                % disp('To continue, Press F5. To exit, press Shift+F5'); keyboard;
                intent = input('Do you wish to quit? Y or [N]: ', 's');
                if isempty(intent)
                    % assign default answer
                    intent = 'n';
                end
                if strcmpi('Y',intent)
                    disp('Quitting..');
                    fclose(hlog_file);
                    fclose(hout_file);
                    return; % exit the program
                else
                    disp('Continuing..');
                    disp ' ';
                end
            end
        end
        
        %% threshold computation
        
        % record the list of heard levels at this frequency
        list_heard_levels(ear_indx, fr, :) = heard_levels;
        
        % compute hearing threshold at current frequency as median of three heard levels
        thresholds(ear_indx, fr) = median(heard_levels);
        
    end
    
    % update the accuracy table across entire experiment
    % n_corr_resp_overall(:,run) = n_corr_resp;
    
    %% Processing plots
    if mod(ear_indx,2)
        h_fig = figure('visible','off');
        h_fig.WindowState = 'maximized';
    end
    
    subplot(1,2,ear_indx);
    semilogx(list_freq, thresholds(ear_indx,:), [col{ear_indx} 'o-']);
    hold on;
    semilogx(list_freq, squeeze(list_heard_levels(ear_indx,:,:)), [col{ear_indx} 'x']);
    xlabel('Pure Tone Frequency (kHz)'); % ' ' factor_suffix{fac} ]);
    ylabel('Hearing Threshold (dB-HL)');
    title([ear ' ear']);
    axis tight;
    grid on;
    ylim([-10 100]);
    set(gca,'XTick',list_freq);
    xlim([min(list_freq)-0.1 max(list_freq)+1]);
    set(gca,'XTickLabel',string(list_freq));
    
    if mod(ear_indx,2)==0
        saveas(h_fig, [out_path_exp '/' out_file_name '_audiogram'], 'png');
        close(h_fig);
        disp(['Writing plots at: ' out_path_exp]);
    end
    
    
    %% display inter run message
    disp ' ';
    fprintf(hlog_file,'\n\n');
    
    % display and pause message except for last run
    if mod(ear_indx,2) && flag_debug==0
        % display message
        disp(inter_ear_statement);
        
        % pause the program until the user is ready
        pause;
        
    end
end

% close the log file and results file at end of experiment
fclose(hlog_file);
fclose(hout_file);

%% This function generates acoustic stimulus using tones for hearing test
function stimulus = generate_tones(tone_freq, ear, samp_rate)

% set the duration of the tone % in sec
tone_dur = 0.2;

% set inter tone interval % in sec
int_tone_interval = 0.100;

% set the RMS level of the tones
rms_level = 0.7;

% set the number of tones in a train
num_tones = 3;

% set default sampling frequency in Hz
if nargin < 3
    samp_rate = 48000;
end

% set default ear 
if nargin < 2
    ear = 'Left';
end

% set default tone frequency
if nargin < 1
    tone_freq = 1000;
end

% set the onset/offset duration of each tone 
onset_duration = 20; % in ms

%% Tone construction

% compute the time indices
t_vec = (0 : 1/samp_rate : tone_dur - 1/samp_rate);

% compute the standard tone
std_tone = sin(2*pi * tone_freq * t_vec);

% set the level and temporal characteristics of the sound

% set the RMS level of standard
std_tone = rms_level * std_tone/std(std_tone);

% shape the onset and offset of the standard tone
std_tone = window_adsr(std_tone, samp_rate, onset_duration);

%% Stimulus construction

int_tone_gap = zeros(round(int_tone_interval * samp_rate),1);

stimulus = [repmat([std_tone; int_tone_gap], num_tones-1, 1); std_tone];

% convert to stereo
if strcmpi(ear,'left')
    stimulus = [stimulus zeros(length(stimulus),1)];
elseif strcmpi(ear,'right')
    stimulus = [zeros(length(stimulus),1) stimulus];
end

end

function out_data = window_adsr(data, samp_rate, window_ms)

% out_data = window_adsr(data, samp_rate, window_ms)
%
% This function applies a raised cosine window at onset and offset of the
% data. This is used to shape attack and decay characterisitics of stimuli
%
% Notes:
% Hanning window a.k.a raised cosine window 
%  h(theta) = 0.5*(1-cos(theta))    theta = 0 to 2*pi
%  h(n) = 0.5*(1-cos(2*pi*n/(N-1)))     n = 0 to N-1
%
% Inputs:
%   data      - input data vector assumes a single column
%   samp_rate - sampling rate in Hz
%   window_ms - attack or decay duration in ms
%
% Outputs:
%   out_data  - output data always with row-time & column-ch format
%

% Version History:-
%  ver - DD MMM YYYY - Feature added
%  1.0 - 31 Jan 2022 - first implementation

%% default values, error checks and handling

if nargin < 2
    % set sampling rate in Hz
    samp_rate = 44100;
end

if nargin < 3
    % set attack, decay duration in ms
    window_ms = 10;
end

if nargin < 1
    % error handling
    disp('Error: no valid input');
    return;
end


%% initializations

% extract the number of channels and length of the input data
if size(data, 1) < size(data, 2)
    n_ch = size(data, 1);
    
    dat_len = size(data, 2);
    
    % intialize output data
    out_data = data.';
else
    n_ch = size(data, 2);
    
    dat_len = size(data, 1);
    
    % intialize output data
    out_data = data;
end

% warn user if in case they are using this function incorrectly
if n_ch > 2
    disp('Warning: unsupported data')
end

% compute length of the window needed
win_len = round(window_ms / 1000 * samp_rate);

% compute the Hann or raised cosine window function
hann_win = ( 1 - cos(linspace(0, 2*pi, 2*win_len).') )/2;


%% processing


% iterate for each channel
for ch = 1:n_ch
    
    % onset window is first half of the Hann window
    onset_win = hann_win(1:win_len);
    
    % apply onset window - attack characteristics
    out_data(1:win_len, ch) = out_data(1:win_len, ch) .* onset_win;
    
    % compute the start index of the window
    start_indx = dat_len - win_len + 1;
    
    % offset window is second half of Hann window
    offset_win = hann_win(win_len+1:end);
    
    % apply offset window - decay characteristics
    out_data(start_indx:dat_len, ch) = out_data(start_indx:dat_len, ch) .* offset_win;
    
end

end