% Xiaoqin May 2019
% Auditory discrimination task for moving binaural sounds

clear
close all
clc
sca

%% Define and save current working directory
expcodeDir = cd;
cd ..
expDir = cd;
resultsDir = [expDir filesep 'Results'];
cd(expcodeDir)

%% Initialisation of experiment 

% Do a series of checks to ensure that data do not get overwritten
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

param.subj = input('Subject Number:');
param.runno = input('Run No (Enter 0 for practice run):');
param.seqno = input('Sequence Order:');

% Regardless of what the actual numeric values were to be associated with the digits, 
% the agreement is that digits were to be arranged in consecutive order, '0', '1', ... up to '9' .
%
% Because this ordering property has been guaranteed for decades, 
% what we are doing here is to convert the character for a digit into the numeric value for the digit 
% We then find the relative position of the digit character to the character for '0', 
% which is easily done by subtracting off the representation of '0' from the representation of the digit character.
param.seqno_array = num2str(param.seqno) - '0';

% Both TEMPORAL and SPATIAL runs will have the same trials presented although
% in randomised order
% Information is just to facilitate analysis and instructions.
param.exp_type = input('Experiment Type (1 = Temporal, 2 = Spatial):');
if param.exp_type == 1
    exp_runtype = 'Temporal';
elseif param.exp_type == 2
    exp_runtype = 'Spatial';
end

subjDir = [resultsDir filesep num2str(param.subj, '%03d')];
subjrunDir = [subjDir filesep 'Run_' num2str(param.runno, '%02d') '_' exp_runtype];
paths.subjresDir = subjrunDir;

% Check for existence of current results directory
while exist(paths.subjresDir, 'dir')
    
    fprintf('\n\n%s exist\n', paths.subjresDir)
    param.ow_resultsDir = input('Overwrite results directory (y = yes, n = no)?', 's');
    
    if strcmp(param.ow_resultsDir, 'n')
        fprintf('\n\nResults directory for Subject %d will not be overwritten... \n\n', param.subj)
        param.subj = input('Re-enter Subject Number:');
        param.runno = input('Re-enter Run No (Enter 0 for practice run):');
        param.exp_type = input('Re-enter Experiment Type (1 = Temporal, 2 = Spatial):');
        if param.exp_type == 1
            exp_runtype = 'Temporal';
        elseif param.exp_type == 2
            exp_runtype = 'Spatial';
        end
        
        subjDir = [resultsDir filesep num2str(param.subj, '%03d')];
        subjrunDir = [subjDir filesep 'Run_' num2str(param.runno, '%02d') '_' exp_runtype];
        paths.subjresDir = subjrunDir;
        
    elseif strcmp(param.ow_resultsDir, 'y')
        fprintf('\n\nResults directory for Subject %d will be overwritten... \n\n', param.subj)
        break
    end
    
end

param.resphand = input('Enter response hand for participant (r = Right, l = Left):', 's');

% Check for response hand entry
while ~strcmp(param.resphand, 'r') && ~strcmp(param.resphand, 'l')
    param.resphand = input('Unrecognised input!!!\nRe-enter response hand for participant (r = Right, l = Left):', 's');
    if strcmp(param.resphand, 'r') || strcmp(param.resphand, 'l')
        break
    end
end

if strcmp(param.resphand, 'r')
    param.resphandtext = 'RIGHT';
elseif strcmp(param.resphand, 'l')
    param.resphandtext = 'LEFT';
end

paths.subjtrialDir = [paths.subjresDir filesep 'Trial_Data']; % Directory for trial data dumping after completion of each trial
if (~exist(paths.subjtrialDir, 'dir'))
    mkdir(paths.subjtrialDir)
end

param.exp_runtype = exp_runtype;

%% Set trial number

param.nprac_trials = 20; 
param.ntrialrep = 1;
param.trialbreak = 24;

%% Set screen and viewing parameters

param.view_dist = 60;                % Viewing distance in cm
param.display_size_x = 53;           % Actual length of screen in cm
param.display_size_y = 30;           % Actual width of screen in cm
param.x_res = 1920;                  % Screen resolution (x)
param.y_res = 1080;                  % Screen resolution (y)
param.screen_no = 0;                 % Screen to present the stimuli on
param.backgrd_col = [0 0 0];         % Screen background color

%% Set response parameters

param.indxresp = 'b';
param.middleresp = 'n';
param.ringresp = 'm';

%% Set Stimuli and Timing paramters

% Set colour and size settings
param.cross_size = 1.00;                        % In visual angle
param.cross_thickness = 0.10;                   % Line thickness for fixation cross (in deg)
param.cross_bg_col = [0 0 0];                   % Background colour that cross should be superimposed on
param.stim_col = [255 255 255];                 % Stimuli colour
param.text_size = 30;                           % Font size for text
param.resp_conf_text_size = 50;                 % Font size for response configuration
param.linespacing = 1.5;                        % Spacing between lines when presenting text

% Sound parameters for defining sound directories
param.sndfreq = 800;                            % Frequency of sound to play
param.audiodist = [400 800 1200 1600];          % Auditory distance array (cm)
param.audiodur = [800 1600 2400 3200];          % Audio duration (ms)
param.startarray = [-1600 -800 0];              % Starting point

% Sound parameters for PsychportAudio
param.sndmode = 2;                              % Device id for opening up PsychPortAudio
param.snd_nchannels = 2;                        % No of sound channels for sound playback
param.snd_sampfreq = 48000;                     % Sampling frequency for sound playback
param.snd_startcue = 0;                         % Defines the starttime of device, 0 = start immediately, > 0, PTB will start device at requested time
param.snd_waitfordevice = 1;                    % 1 = wait until device has really started
param.snd_vol = 1.0;                            % Set sound playback volume, value of 1.0 will pass samples unmodified, 0.5 would reduce intensity by 50%
param.repetitions = 1;                          % Defines how often the sound playback should be repeated. 0 = infinite repeatitions, 1 = play exactly once and stop, 1.5 = 1.5 repetitions

% Set durations
param.fix_cross_dur_range = [0.2 0.5];          % Interval range for fixation cross at start of each trial (in s)
param.ISIrange_std_probe = [1 2];               % ISI between standard and probe (in s)
param.ISIrange_probe_resp = [1 2];              % ISI between probe and response (in s)
param.max_resp_time = 3;                        % Maximum response time (in s)
param.ITIrange = [1 2];                         % Range of Inter-trial interval (in s)
param.wait_end_time = 2;                        % Wait time at end of run (in s)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Specify other directories

param.degchange = 3;
dynbinauralDir = [expDir filesep 'Auditory_Stimuli' filesep 'Dynamic_Binaural' filesep 'Moving_Sounds'];
audiostimDir = [dynbinauralDir filesep num2str(param.sndfreq) '_Hz_' num2str(param.audiodist(1)) 'cm_' num2str(param.audiodur(1)) 'ms_SampF_' num2str(param.snd_sampfreq) filesep 'HRTFconv_DegChange' num2str(param.degchange)];
audiostimlistDir = [expDir filesep 'Stimuli_List'];

%% Save paths in one structure

paths.expDir = expDir;
paths.subjDir = subjDir;
paths.dynbinauralDir = dynbinauralDir;
paths.audiostimDir = audiostimDir;

%% Load stimuli list
if param.runno == 0
    param.stimlistno = param.seqno_array(1);
elseif param.runno ~= 0
    param.stimlistno = param.seqno_array(param.runno);
end

stimlist_fn = ['FTSA_StimListFull_L' num2str(param.stimlistno, '%02d') '.mat'];
stimlist_path = [audiostimlistDir filesep stimlist_fn];
load(stimlist_path)

param.stimlist_import = FTS_audstimlist_sp;
param.ntrials = size(param.stimlist_import,1);

%% Randomise trial sequence 

%%%%%%%%%%%%%%%%%% XQ edited till here %%%%%%%%%%%%%%%%%%%%

seedlist_path = [resultsDir filesep 'Seeds_Used.csv'];

% To randomise seeds in MATLAB so that the starting seed is different each
% time we start up MATLAB
% RandStream.setDefaultStream(RandStream('mt19937ar','seed',sum(100*clock))); % Use for older MATLABs
rng('shuffle') % Use for newer MATLABs

seed_seq = randperm(2^20);
param.seed_used = seed_seq(1);

% Check if file containing seed details exist, if not create file
if ~exist(seedlist_path, 'file')
    
    % Create dataset to store seed used
    seedinfo_all_ds = dataset({param.subj, 'Subject'}, ...
        {param.runno, 'RunNo'}, ...
        {param.exp_type, 'ExpType'}, ...
        {param.seed_used, 'SeedUsed'});
    
    % Write out seed information to .csv file
    export(seedinfo_all_ds, 'File', seedlist_path, 'Delimiter', ',')
    
elseif exist(seedlist_path, 'file')
    
    % If file exist, check if seed number has been used before
    seedinfo_all_ds = dataset('File', seedlist_path, 'Delimiter', ',');
    if sum(logical(seedinfo_all_ds.SeedUsed == param.seed_used)) > 0
        change_seed = 1;
        fprintf('\n Seed No %d Repeated, seed to be changed... \n', param.seed_used)
    else
        change_seed = 0;
    end
    
    % If seed has already been used, then need to change seed...
    while change_seed == 1
        rng('shuffle')
        seed_seq = randperm(2^20);
        param.seed_used = seed_seq(1);
        if sum(logical(seedinfo_all_ds.SeedUsed == param.seed_used)) > 0
            change_seed = 1;
            fprintf('\n Seed No %d Repeated, seed to be changed... \n', param.seed_used)
        else
            change_seed = 0;
        end
    end
    
    % Get new seed info into a dataset
    curr_seed_info = dataset({param.subj, 'Subject'}, ...
        {param.runno, 'RunNo'}, ...
        {param.exp_type, 'ExpType'}, ...
        {param.seed_used, 'SeedUsed'});
    seedinfo_all_ds = [seedinfo_all_ds; curr_seed_info]; % add to current seed information
    
    % Write out seed information to .csv file
    export(seedinfo_all_ds, 'File', seedlist_path, 'Delimiter', ',')
    
end

rng(param.seed_used)
rand_indx = randperm(param.ntrials);
stimlist_rand = param.stimlist_import(rand_indx,:);

%% Create arrays of trial-specific parameters

% Response options
param.resp_options = ...
    {'S', '=', 'L'; ...
    'L', 'S', '='; ...
    '=', 'L', 'S';};
nrespopt = size(param.resp_options,1);
nrep_resp = param.ntrials/nrespopt;
respopt_arr = repmat(param.resp_options, [nrep_resp,1]);
rand_indx_resp = randperm(param.ntrials);
respopt_rand = respopt_arr(rand_indx_resp,:);
param.respopt_rand = respopt_rand;

% Create ISIs and ITIs arrays
param.ISI1_cross_audio = (param.fix_cross_dur_range(2) - param.fix_cross_dur_range(1)).*rand(param.ntrials,1) + param.fix_cross_dur_range(1); % Between cross and appearance of STANDARD
param.ISI2_cross_audio = (param.fix_cross_dur_range(2) - param.fix_cross_dur_range(1)).*rand(param.ntrials,1) + param.fix_cross_dur_range(1); % Between cross and appearance of PROBE
param.ISI_1 = (param.ISIrange_std_probe(2) - param.ISIrange_std_probe(1)).*rand(param.ntrials,1) + param.ISIrange_std_probe(1); % Between Sample and Probe (in s)
param.ISI_2 = (param.ISIrange_probe_resp(2) - param.ISIrange_probe_resp(1)).*rand(param.ntrials,1) + param.ISIrange_probe_resp(1); % Between Probe and Response (in s)
param.ITI = (param.ITIrange(2) - param.ITIrange(1)).*rand(param.ntrials,1) + param.ITIrange(1);   % in s

% Update stimuli list
stimlist_rand = [stimlist_rand ...
    array2table(respopt_rand, 'VariableNames', {'RespIndex', 'RespMiddle', 'RespRing'}) ...
    array2table(param.ISI1_cross_audio, 'VariableNames', {'ISI_Cross_Std'}), ...
    array2table(param.ISI2_cross_audio, 'VariableNames', {'ISI_Cross_Probe'}), ...
    array2table(param.ISI_1, 'VariableNames', {'ISI_Std_Probe'}) ...
    array2table(param.ISI_2, 'VariableNames', {'ISI_Probe_Response'}) ...
    array2table(param.ITI, 'VariableNames', {'ITI'})];
param.stimlist_fullrun = stimlist_rand;

%% Open Screen

Screen('Preference', 'SkipSyncTests', 0) % 1 to skip, 0 to always do, use during script testing only

[window, wRect] = Screen('OpenWindow', param.screen_no, param.backgrd_col);

% check resolution entered to see if it's correct
if wRect(3) ~= param.x_res || wRect(4) ~= param.y_res
    fprintf('Wrong screen resolution specified! Real screen resolution is %d x %d. \n Changing x and y screen resolution specified... \n', wRect(3), wRect(4))
    param.x_res = wRect(3);
    param.y_res = wRect(4);
end

ppd_horz = calc_vis_angle(param.display_size_x, param.x_res, param.view_dist); % pixels per degree
ppd_vert = calc_vis_angle(param.display_size_y, param.y_res, param.view_dist); % pixels per degree
param.ppd = mean([ppd_horz; ppd_vert]); % average ppd

% Enable alpha blending with proper blend-function
Screen('BlendFunction', window, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
param.center_x = wRect(3)/2;
param.center_y = wRect(4)/2;
fps = Screen('FrameRate', window); % frames per second
ifi = Screen('GetFlipInterval', window);
if fps == 0
    fps = 1/ifi;
end

% Raise error if refresh frequency is not 60 Hz
% if ifi ~= 1/60
%     error('\nWARNING: REFRESH FREQUENCY IS NOT 60 HZ, STOPPING SCRIPT TO CHECK TIMING!!!\n')
% end

HideCursor;	% Hide the mouse cursor
Priority(MaxPriority(window));

% Do initial flip to get vbl...
vbl = Screen('Flip', window);

fprintf('\nPresentation screen parameters: FPS = %d, IFI = %d, VBL = %d \n', fps, ifi, vbl)

% Suppress characters into command window while listening is enabled 
% ListenChar(2)

% Store screen parameters
param.screenpresent_param = [fps ifi vbl];
param.actual_screenres = wRect;

%% Make fixation cross

[fix_cross, rect_cross] = CreateFixationCross(window, param);

%% Specify response device and get keycodes

% Response device
if ispc
    kbID = [];
elseif ismac
    kbID = GetKeyboardIndices('Apple Internal Keyboard / Trackpad');
end

param.kbID = kbID;

% Keycodes
KbName('UnifyKeyNames')

param.key_trig = KbName('space'); % For some reason, Psychtoolbox is unable to get input from 'Enter' with Apple keyboard. So space is use for practice
param.key_index = KbName(param.indxresp);
param.key_middle = KbName(param.middleresp);
param.key_ring = KbName(param.ringresp);
param.key_escape = KbName('ESCAPE'); % Escape key

% Setup internal keyboard queues for main task
keysInterest = zeros(1,256);
keysInterest([param.key_trig param.key_index param.key_middle param.key_ring param.key_escape]) = 1;
KbQueueCreate(kbID, keysInterest);

%% Initialise PsychPortAudio

InitializePsychSound(1);
PsychPortAudio('Verbosity',10)
pahandle = PsychPortAudio('Open', -1, [], param.sndmode, param.snd_sampfreq, param.snd_nchannels, 0);
PsychPortAudio('Volume', pahandle, param.snd_vol);

%% Set total trial number

if param.runno == 0
    nrun_trials = param.nprac_trials;
    runtype = 'Practice';
elseif param.runno ~= 0
    nrun_trials = param.ntrials;
    runtype = 'Test';
end

param.nrun_trials = nrun_trials;

%% Create arrays to minimise CPU workload

stim_playlist = cell(nrun_trials,2);
corr_key_arr = zeros([nrun_trials,1]);
corr_resp_arr = cell(nrun_trials,1);
corr_respfin_arr = cell(nrun_trials,1);
trialonset = zeros([nrun_trials,1]);
std_stimonset = zeros([nrun_trials,1]);
std_stimonsetstop = zeros([nrun_trials,1]);
std_stimend = zeros([nrun_trials,1]);
std_stimdur = zeros([nrun_trials,1]);
probe_stimonset = zeros([nrun_trials,1]);
probe_stimonsetstop = zeros([nrun_trials,1]);
probe_stimend = zeros([nrun_trials,1]);
probe_stimdur = zeros([nrun_trials,1]);
t_startresp = zeros([nrun_trials,1]);
t_endresp = zeros([nrun_trials,1]);
subjrespkc = zeros([nrun_trials,1]);
subjresp = cell(nrun_trials,1);
subjrespfin = cell(nrun_trials,1);
RT = zeros([nrun_trials,1]);

%% Save stimuli list and parameters before starting trials

expparam_fn = ['FTSA_PreExpInfo_S' num2str(param.subj, '%03d'), '_Run' num2str(param.runno, '%02d') '_' param.exp_runtype '.mat'];
expparam_path = [paths.subjresDir filesep expparam_fn];
save(expparam_path, 'param', 'paths', 'stimlist_rand')

%% Start of experiment

for i_trial = 1:nrun_trials
    
    % Standard Parameters
    stddur = stimlist_rand.Aud_Standard_Dur(i_trial);
    stddist = stimlist_rand.Aud_Standard_Dist(i_trial);
    stdstart = abs(stimlist_rand.Aud_Std_StartPt(i_trial));
    stdvel = round(stimlist_rand.Aud_Standard_Vel(i_trial));
    
    std_audio_fn = ['HRTFv9_degchange' num2str(param.degchange) '_pos' num2str(stdstart) 'cm_dur' num2str(stddur) 'ms_dist' num2str(stddist) 'cm_vel' num2str(round(stdvel)) 'ms-1.wav'];
    std_audio_fp = [audiostimDir filesep std_audio_fn];
    std_audio = audioread(std_audio_fp);
    stim_playlist{i_trial,1} = std_audio_fp;
    
    % Probe Parameters
    probedur = stimlist_rand.Aud_Probe_Dur(i_trial);
    probedist = stimlist_rand.Aud_Probe_Dist(i_trial);
    probestart = abs(stimlist_rand.Aud_Probe_StartPt(i_trial));
    probevel = round(stimlist_rand.Aud_Probe_Vel(i_trial));
    
    probe_audio_fn = ['HRTFv9_degchange' num2str(param.degchange) '_pos' num2str(probestart) 'cm_dur' num2str(probedur) 'ms_dist' num2str(probedist) 'cm_vel' num2str(round(probevel)) 'ms-1.wav'];
    probe_audio_fp = [audiostimDir filesep probe_audio_fn];
    probe_audio = audioread(probe_audio_fp);
    stim_playlist{i_trial,2} = probe_audio_fp;
    
    % Define trial information path
    trial_info_path = [paths.subjtrialDir filesep 'FTSA_TrialData_S' num2str(param.subj, '%03d') '_Run' num2str(param.runno, '%02d') '_Trial' num2str(i_trial, '%03d') '.mat'];
    
    % Preload buffers to minimise computer load now then audio files are
    % read
    stdbuffer = PsychPortAudio('CreateBuffer', [], std_audio');
    probebuffer = PsychPortAudio('CreateBuffer', [], probe_audio');
    
    % Get no of frames for duration of audio to synchronise fixation cross
    % with it later (Probably just need the probe?)
    nframes_std = fps*stddur/1000;
    nframes_probe = fps*probedur/1000;
    
    % General timing and response parameters
    ISI_crossstd = stimlist_rand.ISI_Cross_Std(i_trial);
    ISI_crossprobe = stimlist_rand.ISI_Cross_Probe(i_trial);
    ISI_std_probe = stimlist_rand.ISI_Std_Probe(i_trial);
    ISI_probe_resp = stimlist_rand.ISI_Probe_Response(i_trial);
    ITI = stimlist_rand.ITI(i_trial);
    
    indexresp = stimlist_rand.RespIndex{i_trial};
    middleresp = stimlist_rand.RespMiddle{i_trial};
    ringresp = stimlist_rand.RespRing{i_trial};
    resp_conf = {indexresp, middleresp, ringresp};
    
    % Print out trial related information
    if i_trial == 1
        fprintf('\n\n\nRUN NO: %d (%s) \nEXPERIMENT TYPE: %d (%s) \n\n\n', ...
            param.runno, runtype, param.exp_type, param.exp_runtype)
    end
    fprintf('\n\n#################### TRIAL NO %d ####################\n', i_trial)
    fprintf('\nStandard Start: %d cm \nStandard Duration: %d ms \nStandard Distance: %d cm \nStandard Velocity: %5.2f ms-1 \nStandard Audio File: %s\n', ... 
        stdstart, stddur, stddist, stdvel, std_audio_fn)
    fprintf('\nProbe Start: %d cm \nProbe Duration: %d ms \nProbe Distance: %d cm \nProbe Velocity: %5.2f ms-1 \nProbe Audio File: %s\n\n', ...
        probestart, probedur, probedist, probevel, probe_audio_fn)
    
    % Correct answer for current trial
    if param.exp_type == 1 % for temporal run, compare duration
        
        dur_diff = probedur - stddur;
        
        if dur_diff == 0
            
            corr_key_indx = find(ismember(resp_conf, '='));
            if corr_key_indx == 1
                corr_key = param.key_index;
                corr_respfin = 'Index';
            elseif corr_key_indx == 2
                corr_key = param.key_middle;
                corr_respfin = 'Middle';
            elseif corr_key_indx == 3
                corr_key = param.key_ring;
                corr_respfin = 'Ring';
            end
            corr_resp = 'Equals to';
            
        elseif dur_diff > 0
            
            corr_key_indx = find(ismember(resp_conf, 'L'));
            if corr_key_indx == 1
                corr_key = param.key_index;
                corr_respfin = 'Index';
            elseif corr_key_indx == 2
                corr_key = param.key_middle;
                corr_respfin = 'Middle';
            elseif corr_key_indx == 3
                corr_key = param.key_ring;
                corr_respfin = 'Ring';
            end
            corr_resp = 'Longer';
            
        elseif dur_diff < 0
            
            corr_key_indx = find(ismember(resp_conf, 'S'));
            if corr_key_indx == 1
                corr_key = param.key_index;
                corr_respfin = 'Index';
            elseif corr_key_indx == 2
                corr_key = param.key_middle;
                corr_respfin = 'Middle';
            elseif corr_key_indx == 3
                corr_key = param.key_ring;
                corr_respfin = 'Ring';
            end
            corr_resp = 'Shorter';
            
        end
        
    elseif param.exp_type == 2 % for spatial run, compare distance
        
        dist_diff = probedist - stddist;
        
        if dist_diff == 0
            
            corr_key_indx = find(ismember(resp_conf, '='));
            if corr_key_indx == 1
                corr_key = param.key_index;
                corr_respfin = 'Index';
            elseif corr_key_indx == 2
                corr_key = param.key_middle;
                corr_respfin = 'Middle';
            elseif corr_key_indx == 3
                corr_key = param.key_ring;
                corr_respfin = 'Ring';
            end
            corr_resp = 'Equals to';
            
        elseif dist_diff > 0
            
            corr_key_indx = find(ismember(resp_conf, 'L'));
            if corr_key_indx == 1
                corr_key = param.key_index;
                corr_respfin = 'Index';
            elseif corr_key_indx == 2
                corr_key = param.key_middle;
                corr_respfin = 'Middle';
            elseif corr_key_indx == 3
                corr_key = param.key_ring;
                corr_respfin = 'Ring';
            end
            corr_resp = 'Longer';
            
        elseif dist_diff < 0
            
            corr_key_indx = find(ismember(resp_conf, 'S'));
            if corr_key_indx == 1
                corr_key = param.key_index;
                corr_respfin = 'Index';
            elseif corr_key_indx == 2
                corr_key = param.key_middle;
                corr_respfin = 'Middle';
            elseif corr_key_indx == 3
                corr_key = param.key_ring;
                corr_respfin = 'Ring';
            end
            corr_resp = 'Shorter';
            
        end
    end
    
    fprintf('### TRIAL %d ###\nCorrect Response: %s\nCorrect Response Finger: %s\n', i_trial, corr_resp, corr_respfin)
    corr_key_arr(i_trial,1) = corr_key;
    corr_resp_arr{i_trial,1} = corr_resp;
    corr_respfin_arr{i_trial,1} = corr_respfin;
    
    % Send text parameters to screen mex function
    Screen('TextFont', window, 'Helvetica');
    Screen('TextSize', window, param.text_size);
    
    % Fill playbuffer with content of STANDARD audio:
    PsychPortAudio('FillBuffer', pahandle, stdbuffer);
    
    % Flush keyboard, start collection of queue
    KbQueueFlush(kbID);
    KbQueueStart(kbID);
    
    %%%%%%%%%%%%%%%%%%%%% TO RUN JUST FOR THE FIRST TRIAL %%%%%%%%%%%%%%%%%
    if i_trial == 1
        
        % Wait for input from scanner to trigger trial
        DrawFormattedText(window, sprintf('RUN %d START\n\n\n\nPress <SPACE> to continue.', param.runno), 'center', 'center', param.stim_col);
        [vbl] = Screen('Flip', window, vbl+0.5*ifi); %PTB should be ready to draw screen at half an IFI
        
        keyisdown = 0;
        while(~keyisdown)
            [keyisdown, firstpress] = KbQueueCheck(kbID);
            if keyisdown ~= 0
                if firstpress(param.key_trig)
                    keyisdown = 1;
                elseif firstpress(param.key_escape)
                    KbQueueRelease; % Release keyboard internal queue
                    PsychPortAudio('Close', pahandle); % Close the audio device
                    sca % close screen
                    Priority(0);
                    ShowCursor;
                    error('Escape key detected. Program exited!!!')
                else
                    keyisdown = 0;
                end
            end
        end
        
        % Clear the queue
        KbQueueFlush(kbID);
        
        % Some experiment instructions
        exp_info_text = sprintf(['In this experiment, you will be asked to make\neither DURATION or DISTANCE judgments of sounds.\n\n'...
            'For each trial, TWO moving sounds travelling from left to right\nwould be played to you consecutively.\n\n' ...
            'Always compare the SECOND sound to the FIRST sound\nwhen making your judgments.\n\n'...
            'Kindly maintain fixation on the presented fixation cross during the trial.\n\n\n' ...
            'Press <SPACE> to continue.'], upper(exp_runtype));
        DrawFormattedText(window, exp_info_text, 'center', 'center', param.stim_col, [], [], [], param.linespacing);
        [vbl] = Screen('Flip', window, vbl+0.5*ifi); %PTB should be ready to draw screen at half an IFI
        
        keyisdown = 0;
        while(~keyisdown)
            [keyisdown, firstpress] = KbQueueCheck(kbID);
            if keyisdown ~= 0
                if firstpress(param.key_trig)
                    keyisdown = 1;
                elseif firstpress(param.key_escape)
                    KbQueueRelease; % Release keyboard internal queue
                    PsychPortAudio('Close', pahandle); % Close the audio device
                    sca % close screen
                    Priority(0);
                    ShowCursor;
                    error('Escape key detected. Program exited!!!')
                else
                    keyisdown = 0;
                end
            end
        end
        
        % Clear the queue
        KbQueueFlush(kbID);
        
        % Run Information
        if param.exp_type == 1
            run_info_text = sprintf(['\n\nThe current run is a %s block.\n\nCompare the DURATION of the second sound stimulus to the first.\n\n'...
                'If second sound is SHORTER than first, choose S\n' ...
                'If second sound is EQUALS TO first, choose =\n'...
                'If second sound is LONGER than first, choose L'...
                '\n\nPlace the index, ring and middle finger of your %s hand\non "B", "N" and "M" respectively.' ...
                '\n\nResponse configuration would be shown at the end of the trial.\n\n\n' ...
                'Press <SPACE> to start the block.'], upper(exp_runtype), param.resphandtext);
        elseif param.exp_type == 2
            run_info_text = sprintf(['\n\n The current run is a %s block. \n\n Compare the DISTANCE TRAVELLED by the second sound stimulus to the first.\n\n'...
                'If second sound is SHORTER than first, choose S\n' ...
                'If second sound is EQUALS TO first, choose =\n'...
                'If second sound is LONGER than first, choose L'...
                '\n\nPlace the index, ring and middle finger of your %s hand\non "B", "N" and "M" respectively.' ...
                '\n\nResponse configuration would be shown at the end of the trial.\n\n\n' ...
                'Press <SPACE> to start the block.'], upper(exp_runtype), param.resphandtext);
        end
        
        DrawFormattedText(window, run_info_text, 'center', 'center', param.stim_col, [], [], [], param.linespacing);
        [vbl] = Screen('Flip', window, vbl+0.5*ifi);
        
        keyisdown = 0;
        while(~keyisdown)
            [keyisdown, firstpress] = KbQueueCheck(kbID);
            if keyisdown ~= 0
                if firstpress(param.key_trig)
                    pressed_keycode = find(firstpress > 0);
                    blk_start = firstpress(pressed_keycode);
                    keyisdown = 1;
                elseif firstpress(param.key_escape)
                    KbQueueRelease; % Release keyboard internal queue
                    PsychPortAudio('Close', pahandle); % Close the audio device
                    sca % close screen
                    Priority(0);
                    ShowCursor;
                    error('Escape key detected. Program exited!!!')
                else
                    keyisdown = 0;
                end
            end
        end
        
        % Clear the queue
        KbQueueFlush(kbID);        
        param.t0 = blk_start; % Start time for the entire block

    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Fixation Cross (STANDARD)
    Screen('FillRect', window, param.backgrd_col);
    Screen('DrawTexture', window, fix_cross, [], rect_cross);
    [vbl, crossonset] = Screen('Flip', window, vbl+0.2*ifi, 1);
    WaitSecs(ISI_crossstd);
    
    trialonset(i_trial,1) = crossonset;
    
    % STANDARD
    snd_start = PsychPortAudio('Start', pahandle, param.repetitions, param.snd_startcue, param.snd_waitfordevice);
    
    % Fixation cross will be presented during STANDARD playback
    for i_frame = 1:nframes_std
        Screen('FillRect', window, param.backgrd_col);
        Screen('DrawTexture', window, fix_cross, [], rect_cross);
        vbl = Screen('Flip', window, vbl+0.2*ifi, 1);
    end
    
    % Stop STANDARD playback
    [actualStartTime, endPositionSecs, xruns, estStopTime] = PsychPortAudio('Stop', pahandle, 1, 1);
    
    std_stimonset(i_trial,1) = snd_start;
    std_stimonsetstop(i_trial,1) = actualStartTime;
    std_stimend(i_trial,1) = estStopTime;
    std_stimdur(i_trial,1) = endPositionSecs;
    
    CheckForEscape(kbID, param.key_escape, pahandle)
    
    % Fill buffer with PROBE audio
    PsychPortAudio('FillBuffer', pahandle, probebuffer);
        
    % ISI between STANDARD and PROBE
    Screen('FillRect', window, param.backgrd_col);
    Screen('DrawTexture', window, fix_cross, [], rect_cross);
    vbl = Screen('Flip', window, vbl+0.2*ifi, 1);
    WaitSecs(ISI_std_probe);

%     % Fixation Cross (PROBE)
%     Screen('FillRect', window, param.backgrd_col);
%     Screen('DrawTexture', window, fix_cross, [], rect_cross);
%     vbl = Screen('Flip', window, vbl+0.2*ifi, 1);
%     WaitSecs(ISI_crossprobe);
    
    % PROBE
    snd_start = PsychPortAudio('Start', pahandle, param.repetitions, param.snd_startcue, param.snd_waitfordevice);
    
    % Fixation cross will be presented during STANDARD playback
    for i_frame = 1:nframes_probe
        Screen('FillRect', window, param.backgrd_col);
        Screen('DrawTexture', window, fix_cross, [], rect_cross);
        vbl = Screen('Flip', window, vbl+0.2*ifi, 1);
    end
    
    % Stop PROBE playback
    [actualStartTime, endPositionSecs, xruns, estStopTime] = PsychPortAudio('Stop', pahandle, 1, 1);
    
    probe_stimonset(i_trial,1) = snd_start;
    probe_stimonsetstop(i_trial,1) = actualStartTime;
    probe_stimend(i_trial,1) = estStopTime;
    probe_stimdur(i_trial,1) = endPositionSecs;
    
    CheckForEscape(kbID, param.key_escape, pahandle)
    
    % SOA between PROBE and RESPONSE
    Screen('FillRect', window, param.backgrd_col);
    Screen('DrawTexture', window, fix_cross, [], rect_cross);
    vbl = Screen('Flip', window, vbl+0.2*ifi, 1);
    WaitSecs(ISI_probe_resp);
    
    % Clear screen
    Screen('FillRect', window, param.backgrd_col);
    vbl = Screen('Flip', window, vbl+0.2*ifi, 1);
    WaitSecs(0.1);
    
    % Response
    Screen('TextSize', window, param.resp_conf_text_size);
    resp_text = sprintf('%s       %s       %s', resp_conf{1}, resp_conf{2}, resp_conf{3});
    DrawFormattedText(window, resp_text, 'center', 'center', param.stim_col);
    [vbl, resp_starttime] = Screen('Flip', window, vbl+0.5*ifi);
    
    % Flush Queue before capturing response
    KbQueueFlush(kbID);
    
    t_startresp(i_trial,1) = resp_starttime;
    
    % Capturing response
    t_endresp_temp = resp_starttime;
    keyisdown = 0;
    while ~keyisdown
        t_endresp_temp = GetSecs;
        [keyisdown, firstpress, firstrelease, lastpress, lastrelease] = KbQueueCheck(kbID);
        if t_endresp_temp - t_startresp > param.max_resp_time % End response section if more than maximum response time
            pressed_keycode = find(firstpress > 0);
            keyisdown = 1;
        elseif keyisdown
            if firstpress(param.key_index) || firstpress(param.key_middle) || firstpress(param.key_ring)
                pressed_keycode = find(firstpress > 0);
                curr_keypress_time = firstpress(pressed_keycode);
                keyisdown = 1;
            elseif firstpress(param.key_escape)
                KbQueueRelease; % Release keyboard internal queue
                PsychPortAudio('Close', pahandle); % Close the audio device
                sca % close screen
                Priority(0);
                ShowCursor;
                error('Escape key detected. Program exited!!!')
            else
                keyisdown = 0;
            end
        end
    end
    
    if isempty(pressed_keycode)
        pressed_keycode = 999;
        curr_keypress_time = NaN;
    end
    
    if pressed_keycode == param.key_index
        subjrespfin{i_trial,1} = 'Index';
        subjresp{i_trial,1} = resp_conf{1};
    elseif pressed_keycode == param.key_middle
        subjrespfin{i_trial,1} = 'Middle';
        subjresp{i_trial,1} = resp_conf{2};
    elseif pressed_keycode == param.key_ring
        subjrespfin{i_trial,1} = 'Ring';
        subjresp{i_trial,1} = resp_conf{3};
    elseif pressed_keycode == 999
        subjrespfin{i_trial,1} = 'Missed';
        subjresp{i_trial,1} = 'Missed';
    end
    
    t_endresp(i_trial,1) = curr_keypress_time;          
    subjrespkc(i_trial,1) = pressed_keycode;
    RT(i_trial,1) = curr_keypress_time - resp_starttime;
    
    % Print accuracy data in command window
    if pressed_keycode == 999
        fprintf('\n### TRIAL %d ###\nResponse Keycode: %d, MISSED response!!!', i_trial, pressed_keycode)
    elseif pressed_keycode == corr_key
        fprintf('\n### TRIAL %d ###\nResponse Keycode: %d, CORRECT response!!!', i_trial, pressed_keycode)
    elseif pressed_keycode ~= corr_key
        fprintf('\n### TRIAL %d ###\nResponse Keycode: %d, INCORRECT response!!!', i_trial, pressed_keycode)
    end
    
    if mod(i_trial,param.trialbreak) == 0
        % Clear the queue
        KbQueueFlush(kbID);
    else
        % Clear the queue
        KbQueueFlush(kbID);
        KbQueueStop(kbID);
    end
    
    % Quickly save trial information
    save(trial_info_path, 'corr_key_arr', 'corr_resp_arr', 'corr_respfin_arr', ...
        'std_stimonset', 'std_stimonsetstop', 'std_stimend', 'std_stimdur', ...
        'probe_stimonset', 'probe_stimonsetstop', 'probe_stimend', 'probe_stimdur', ...
        't_startresp', 't_endresp', 'subjresp', 'RT')
    
    % ITI
    Screen('TextSize', window, param.text_size);
    
    if mod(i_trial,param.trialbreak) == 0 && i_trial ~= nrun_trials % break within the block
        
        breakblk_text = sprintf('Break!\n\n\nPress <SPACE> to continue.');
        DrawFormattedText(window, breakblk_text, 'center', 'center', param.stim_col);
        [vbl] = Screen('Flip', window, vbl+0.5*ifi);
        
        keyisdown = 0;
        while(~keyisdown)
            [keyisdown, firstpress] = KbQueueCheck(kbID);
            if keyisdown ~= 0
                if firstpress(param.key_trig)
                    keyisdown = 1;
                elseif firstpress(param.key_escape)
                    KbQueueRelease; % Release keyboard internal queue
                    PsychPortAudio('Close', pahandle); % Close the audio device
                    sca % close screen
                    Priority(0);
                    ShowCursor;
                    error('Escape key detected. Program exited!!!')
                else
                    keyisdown = 0;
                end
            end
        end
        
        % Clear the queue
        KbQueueFlush(kbID);
        KbQueueStop(kbID);
        
    elseif i_trial == nrun_trials && param.runno == 4
        
        endofblk_text = sprintf('END OF EXPERIMENT.\n\nThe experimenter would attend to you shortly\nand provide further instructions.\n\nPlease wait patiently till then.');
        DrawFormattedText(window, endofblk_text, 'center', 'center', param.stim_col);
        [vbl] = Screen('Flip', window, vbl+0.5*ifi);
        WaitSecs(param.wait_end_time)
        
    elseif i_trial == nrun_trials && param.runno ~= 4
        
        endofblk_text = sprintf('END OF BLOCK.\n\nPlease take a rest.\n\nWe will continue with the next block shortly.');
        DrawFormattedText(window, endofblk_text, 'center', 'center', param.stim_col);
        [vbl] = Screen('Flip', window, vbl+0.5*ifi);
        WaitSecs(param.wait_end_time)
        
    elseif i_trial ~= nrun_trials
        
        Screen('FillRect', window, param.backgrd_col);
        [vbl] = Screen('Flip', window, vbl+0.5*ifi);
        WaitSecs(ITI)
        
    end
        
end

%% End routine

% ListenChar(1);

% Release keyboard internal queue
KbQueueStop(kbID);
KbQueueRelease;

% Close the audio device
PsychPortAudio('Close', pahandle);

sca % close screen
Priority(0);
ShowCursor;

%% Report accuracy

param.corr_indx = logical(corr_key_arr == subjrespkc);
param.accuracy = sum(param.corr_indx);

fprintf('\nACCURACY FOR RUN %d: %d/%d, %5.2f%% \n\n', param.runno, param.accuracy, nrun_trials, (param.accuracy/nrun_trials)*100)

%% Save responses

% Define filepath
results_matpath = [paths.subjresDir filesep 'FTSA_ExpData_S' num2str(param.subj, '%03d') '_Run' num2str(param.runno, '%02d') '_' param.exp_runtype '.mat'];
results_csvpath = [paths.subjresDir filesep 'FTSA_ExpData_S' num2str(param.subj, '%03d') '_Run' num2str(param.runno, '%02d') '_' param.exp_runtype '.csv'];

% Concatanate all necessary variables
blk_char_resp = horzcat(respopt_rand(1:nrun_trials,:), corr_resp_arr, corr_respfin_arr, subjresp, subjrespfin);
blk_subj_resp_array = [repmat(param.subj,[nrun_trials,1]) repmat(param.runno,[nrun_trials,1]) repmat(param.exp_type,[nrun_trials,1]) (1:1:nrun_trials)' ...
    trialonset std_stimonset std_stimonsetstop std_stimend std_stimdur ...
    probe_stimonset probe_stimonsetstop probe_stimend probe_stimdur ... 
    t_startresp t_endresp RT ...
    corr_key_arr subjrespkc];
blk_subj_resp_tbl = [array2table(repmat(param.subj,[nrun_trials,1]), 'VariableNames', {'Subject'}) ...
    cell2table(repmat({param.exp_runtype}, [nrun_trials,1]), 'VariableNames', {'Task'}) ...
    array2table(repmat(param.seqno,[nrun_trials,1]), 'VariableNames', {'SequenceOrder'}) ...
    array2table(repmat(param.stimlistno,[nrun_trials,1]), 'VariableNames', {'ListNo'}) ...
    array2table(repmat(param.runno,[nrun_trials,1]), 'VariableNames', {'RunNo'}) ...
    array2table((1:1:nrun_trials)', 'VariableNames', {'TrialNo'}) ...
    stimlist_rand(1:nrun_trials,8:end) ...
    array2table(blk_subj_resp_array(:,5:end-1), 'VariableNames', {'TrialOnset', 'StdOnset', 'StdOnsetByStop', 'StdOffset', 'StdDur', ...
    'ProbeOnset', 'ProbeOnsetByStop', 'ProbeOffset', 'ProbeDur', 'RespWinStart', 'RespTime', 'RT', 'CorrectKC'}) ...
    array2table(corr_resp_arr, 'VariableNames', {'CorrectResponse'}) ...
    array2table(corr_respfin_arr, 'VariableNames', {'CorrectFinger'}) ...
    array2table(subjrespkc, 'VariableNames', {'SubjRespKC'}) ...
    array2table(subjresp, 'VariableNames', {'SubjResponse'}) ...
    array2table(subjrespfin, 'VariableNames', {'SubjRespFinger'})];

% Write out results
save(results_matpath, 'param', 'paths', 'stim_playlist', 'blk_char_resp', 'blk_subj_resp_array', 'blk_subj_resp_tbl')
writetable(blk_subj_resp_tbl, results_csvpath, 'Delimiter', ',')



