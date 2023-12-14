%% (OPTIONAL) Run in server mode: Set Brainstorm nogui on completely headless mode (Baseline) 

% Set up the Brainstorm files
clear
addpath('/software/path/matlab/brainstorm3_v20220706');
BrainstormDbDir = '/private/path/cMMN/brainstorm_db';
% Start Brainstorm
if ~brainstorm('status')
    brainstorm server
end
bst_set('BrainstormDbDir',BrainstormDbDir)
% Select the correct protocol
ProtocolName = 'cMMN_EEG_SOA'; % Enter the name of your protocol
sProtocol.Comment = ProtocolName;
sProtocol.SUBJECTS = '/private/path/cMMN/brainstorm_db/cMMN_EEG_SOA/anat';
sProtocol.STUDIES = '/private/path/cMMN/brainstorm_db/cMMN_EEG_SOA/data';
db_edit_protocol('load',sProtocol);
% Get the protocol index
iProtocol = bst_get('Protocol', ProtocolName);
if isempty(iProtocol)
    error(['Unknown protocol: ' ProtocolName]);
end
% Select the current procotol
gui_brainstorm('SetCurrentProtocol', iProtocol);

%% Define variables cMMN_Baseline 

% Directories for server in local computer
running_in = 'local'; % 'server' 'local'

% Define paths
if strcmp(running_in,'server')
    root_dir = '/private/path/cMMN/User/cMMN_EEG_SOA'; 
    root_dir_bs = '/private/path/cMMN/brainstorm_db/cMMN_EEG_SOA'; 
    anat_path = '/private/path/cMMN/brainstorm_db/cMMN_EEG_SOA';
elseif strcmp(running_in,'local')
    root_dir = 'C:/private/path/User/cMMN_EEG_SOA'; 
    root_dir_bs = 'C:/private/path/brainstorm_db/cMMN_EEG_SOA'; 
    anat_path = 'C:/private/path/brainstorm_db/cMMN_EEG_SOA';
end
addpath([root_dir '/Scripts']);

% Define protocol name, in case we don't run it with server mode
ProtocolName = 'cMMN_EEG_SOA'; % 'cMMN_EEG' or 'cMMN_followup'
load([root_dir '/subject_array.mat'])

% Define participants
participant = {subject_array{:,1}};
% Whatever you want the GAVR folder to be named
gavr_name = 'GAVR_cMMN_SOA'; % 'GAVR_cMMN' 'GAVR_cMMN_avgref'; 
participant_group = {'C'};

% Experiment tasks
Exp_tasks = {'dualruleMMN_long_SOA'};
% Standard MMN Conditions
condition_standardMMN = {'1','2','3'};
condition_names_standardMMN = {'Standard','DeviantPitch','DeviantDuration'};
condition_mismatch_names_standardMMN = {'Standard','DeviantPitch','DeviantDuration','DeviantPitch-Standard','DeviantDuration-Standard'};
condition_short_labels_standardMMN = {'STD','PDev','DDev','Pitch_MMN','Dur_MMN'};
condition_mismatch_names_standardMMN_peaks = {'Standard_PD','Standard_DD','DeviantPitch','DeviantDuration','DeviantPitch-Standard','DeviantDuration-Standard'};
condition_short_labels_standardMMN_peaks = {'STD_PD','STD_DD','PDev','DDev','Pitch_MMN','Dur_MMN'};
% extratoneMMN Conditions
condition_extratoneMMN = {'1','2','3','4'};
condition_names_extratoneMMN = {'Standard1','Standard2','Standard3','DeviantExtra'};
condition_mismatch_names_extratoneMMN = {'Standard1','Standard2','Standard3','DeviantExtra','DeviantExtra-Standard1','DeviantExtra-Standard2','DeviantExtra-Standard3'};
condition_short_labels_extratoneMMN = {'STD1','STD2','STD3','DExtra','Extra_MMN1','Extra_MMN2','Extra_MMN3'};
% dualruleMMN Conditions
condition_dualruleMMN = {'1_12','1_18'};
condition_names_dualruleMMN = {'ABA','AAA'};
condition_mismatch_names_dualruleMMN = {'ABA','AAA','AAA-ABA'};
condition_short_labels_dualruleMMN = {'ABA','AAA','AAASABA'};
% Long SOA MMN Conditions
condition_long_SOA = {'1_12','1_18'};
condition_names_long_SOA = {'ABA','AAA'};
condition_mismatch_names_long_SOA = {'ABA','AAA','AAA-ABA'};
condition_short_labels_long_SOA = {'ABA','AAA','AAASABA'};
% All MMN-only conditions (for covariance trials)
condition_names_MMNonly_standardMMN = {'DeviantPitch-Standard','DeviantDuration-Standard'};
condition_names_MMNonly_extratoneMMN = {'DeviantExtra-Standard1','DeviantExtra-Standard2','DeviantExtra-Standard3'};
condition_names_MMNonly_dualruleMMN = {'DeviantDoubleRule-LeftPitchA','DeviantDoubleRule-RightPitchB'};
condition_names_MMNonly_long_SOA = {'DeviantDoubleRule_long_SOA-LeftPitchA_long_SOA','DeviantDoubleRule_long_SOA-RightPitchB_long_SOA'};
modality_data = {'EEG'};
epoch_wave_standardMMN = [-0.1, 0.3];
epoch_wave_extratoneMMN = [-0.1, 0.7];
epoch_wave_dualruleMMN = [-0.05, 0.950];
epoch_wave_long_SOA = [-0.05, 0.950];
epoch_baseline = [-0.05, 0]; % LLR baseline correction
time_noise_cov  = [-0.1, 0]; % time window from which to obtain noise covariance values
reject_EEG_absolute = [0, 50]; % absolute threshold
reject_MEG_GRAD_absolute = [0, 2500];
reject_MEG_MAG_absolute = [0, 2500];
LLR_highpass = 0; 
LLR_lowpass = 20; %
tranband_LLR = 0; % 
crit_sweeps = 30; % minimum number of surviving sweeps to discard EEG, MEG or BIMODAL
crit_percent = 50; % minimum percentage of surviving sweeps to discard EEG, MEG or BIMODAL
reref_option = 1; % 0 = NO 1 = YES Yes or no rereference
run_matlab_command = 1; % 0 = NO 1 = YES SELECT NO IF IMPORTING RAW DATA (No tsss or AMICA)
notch_filter = 0; % Apply a notch filter if indicated for specific subjects
notch_alpha = 10;
ref_EEG = 'M1, M2'; % Alternative: 'AVERAGE'
sensor_analysis = 1; % 1 = EEG only; 2 = MEG only; 3 = Combined EEG and MEG only; 4 = ALL
source_analysis = 2; % 1 = EEG only; 2 = MEG only; 3 = Combined EEG and MEG only. 4 = ALL
source_space = {'Surf'}; % Which sources to average (Volume,  Surface or both) 'Vol','Surf'
source_noise_option = 'reg'; % For 'NoiseMethod'; 'median'; 'diag'
source_noise_tag = 'Source_Regul'; % 'Source_Eigen' 'Source_diag'
source_orient_constrain = 0.2; % Could do orthogonal (0.0), or 0.2 (recommended in pipeline meeting)
source_inverse_measure = 'dspm2018'; % 'amplitude' or 'dspm2018'
delete_previous_file = 1; % 1 = delete previous head models and source kernels if reran these steps
compute_covariance = 1; % In case we want to repeat sources but not recompute noise covariance matrices
name_volume_source_grid = 'headmodel_grid.mat';
group_default_cortex = 'tess_cortex_pial_02.mat'; % cortical surface where all cortical sources of invididual subjects will be projected
group_default_volume = 'tess_cortex_mixed_02.mat'; % volume surface with brainstem and cerebellum

% If selected more than one, they will be averaged as a cluster
choice_channel_EEG = {'FCz'}; % Watch out because others may be missing for some subjects
choice_channel_MEG = {'MEG2411'};
Peaks_to_extract_standardMMN = {'MMN'};
Peaks_to_extract_extratoneMMN = {'FirstMMN','SecondMMN'};
Peaks_to_extract_dualruleMMN = {'FirstMMN','SecondMMN','ThirdMMN'};
% Time Windows for all conditions
for exptask = 1:length(Exp_tasks)
    % Define trigger numbers
    if strcmp(Exp_tasks{exptask},'standardMMN')
        condition_mismatch_names = condition_mismatch_names_standardMMN_peaks;  
        condition_short_labels = condition_short_labels_standardMMN_peaks;
        Peaks_to_extract = Peaks_to_extract_standardMMN;
        pitch_time_window = [120-50 120+50];
        duration_time_window = [160-50 160+50];
    elseif strcmp(Exp_tasks{exptask},'extratoneMMN')
        condition_mismatch_names = condition_mismatch_names_extratoneMMN;
        condition_short_labels = condition_short_labels_extratoneMMN;
        Peaks_to_extract = Peaks_to_extract_extratoneMMN;
        unique_time_window = {[145-50 145+50],[425-50 425+50]};
    elseif strcmp(Exp_tasks{exptask},'dualruleMMN')
        condition_mismatch_names = condition_mismatch_names_dualruleMMN;
        condition_short_labels = condition_short_labels_dualruleMMN;
        Peaks_to_extract = Peaks_to_extract_dualruleMMN;
        unique_time_window = {[180-50 180+50],[478-50 478+50],[776-50 776+50]};
    elseif contains(Exp_tasks{exptask},'long_SOA')
        condition_mismatch_names = condition_mismatch_names_long_SOA;
        condition_short_labels = condition_short_labels_long_SOA;
        Peaks_to_extract = Peaks_to_extract_dualruleMMN;
        unique_time_window = {[180-50 180+50],[478-50 478+50],[776-50 776+50]};
    else
        error(['No identifiable task']);
    end
    for mode = 1:length(modality_data)  
        for cmmn = 1:length(condition_mismatch_names)
        for pte = 1:length(Peaks_to_extract)
            if strcmp(Exp_tasks{exptask},'standardMMN') % For cases with different windows per condition
                if contains(condition_mismatch_names{cmmn},'DeviantPitch') % Both deviant and MMN
                    unique_time_window = pitch_time_window;
                elseif contains(condition_mismatch_names{cmmn},'DeviantDuration') % Both deviant and MMN
                    unique_time_window = duration_time_window;
                elseif strcmp(condition_mismatch_names{cmmn},'Standard_PD') % Standard with pitch deviant time window
                    unique_time_window = pitch_time_window;
                elseif strcmp(condition_mismatch_names{cmmn},'Standard_DD') % Standard with duration deviant time window
                    unique_time_window = duration_time_window;
                end
                eval(['time_window_scalp_' Peaks_to_extract{pte} '_' condition_short_labels{cmmn} '_' modality_data{mode} ' = unique_time_window;'])
                eval(['time_window_scalp_' Peaks_to_extract{pte} '_' condition_short_labels{cmmn} '_' modality_data{mode} ' = unique_time_window;'])
                eval(['time_window_source_left_' Peaks_to_extract{pte} '_' condition_short_labels{cmmn} ' = unique_time_window;'])
                eval(['time_window_source_right_' Peaks_to_extract{pte} '_' condition_short_labels{cmmn} ' = unique_time_window;'])
            else % For cases where all conditions within the task will have the same MMN time window
                eval(['time_window_scalp_' Peaks_to_extract{pte} '_' condition_short_labels{cmmn} '_' modality_data{mode} ' = unique_time_window{pte};'])
                eval(['time_window_scalp_' Peaks_to_extract{pte} '_' condition_short_labels{cmmn} '_' modality_data{mode} ' = unique_time_window{pte};'])
                eval(['time_window_source_left_' Peaks_to_extract{pte} '_' condition_short_labels{cmmn} ' = unique_time_window{pte};'])
                eval(['time_window_source_right_' Peaks_to_extract{pte} '_' condition_short_labels{cmmn} ' = unique_time_window{pte};'])
            end
        end   
        end
    end
end

participants_long_SOA = {'2786A','2791A','2792A','2793A','2794A','2801A','2804A','2809A','2814A','2815A','2821A','2822A','2823A'};

initialVars = who; % variables up until here, which won't be deleted afterwards
initialVars = who; % twice so that InitialVars itself is not deleted

%% Import raw data and run PSD (ALWAYS IN CNRL 99 OR 98)
 
tic
disp(' ');      
disp('-------------------------');  
disp('IMPORTING EEG/MEG DATA FOR cMMN');
disp(datetime)
disp('-------------------------');     
disp(' '); 

for p = 1:length(participant)
    % Check log info about the subject
    pos_subj = find(strcmp({subject_array{:,1}},participant{p})); %#ok<*CCAT1>
    if strcmp(subject_array{pos_subj,3},'needs_import')


    folders = dir(['/private/path/cMMN/analysis/' participant{p} '/']);
    infolder = find(endsWith({folders.name},'tsss_AMICA.fif') & (contains({folders.name},[participant{p} '_dualruleMMN_long_SOA'])));
        if isempty(infolder)
            error(['No cMMN AMICA files for ' participant{p}]);
        end

        for i = 1:length(infolder)
            line = infolder(i);
            file_name = ['/private/path/cMMN/analysis/' participant{p} '/' folders(line).name];  %#ok<*SAGROW>             
            
            disp(' ');      
            disp('-------------------------');  
            disp(['Importing data EEG/MEG cMMN data for ' participant{p}]);
            disp(datetime)
            disp(' '); 

            sFiles = [];
            % Process: Create link to raw file
            sFiles = bst_process('CallProcess', 'process_import_data_raw', sFiles, [], ...
                'subjectname',    participant{p}, ...
                'datafile',       {file_name, 'FIF'}, ...
                'channelreplace', 0, ...
                'channelalign',   0, ...
                'evtmode',        'value');   
            
%             disp(' ');      
%             disp('-------------------------');  
%             disp(['Calculating PSD for bad chans for ' participant{p}]);
%             disp(datetime)
%             disp(' '); 
%             
%             % Process: Power spectrum density (Welch) EEG ELECTRODES
%             sFiles_EEG = bst_process('CallProcess', 'process_psd', sFiles, [], ...
%                 'timewindow',  [], ...
%                 'win_length',  1, ...
%                 'win_overlap', 50, ...
%                 'sensortypes', 'EEG', ...
%                 'win_std',     0, ...
%                 'edit',        struct(...
%                      'Comment',         'Power', ...
%                      'TimeBands',       [], ...
%                      'Freqs',           [], ...
%                      'ClusterFuncTime', 'none', ...
%                      'Measure',         'power', ...
%                      'Output',          'all', ...
%                      'SaveKernel',      0));
% 
%             % Process: Power spectrum density (Welch) EEG ELECTRODES
%             sFiles_MEG_GRAD = bst_process('CallProcess', 'process_psd', sFiles, [], ...
%                 'timewindow',  [], ...
%                 'win_length',  1, ...
%                 'win_overlap', 50, ...
%                 'sensortypes', 'MEG GRAD', ...
%                 'win_std',     0, ...
%                 'edit',        struct(...
%                      'Comment',         'Power', ...
%                      'TimeBands',       [], ...
%                      'Freqs',           [], ...
%                      'ClusterFuncTime', 'none', ...
%                      'Measure',         'power', ...
%                      'Output',          'all', ...
%                      'SaveKernel',      0));
% 
% 
%             % Process: Power spectrum density (Welch) EEG ELECTRODES
%             sFiles_MEG_MAG = bst_process('CallProcess', 'process_psd', sFiles, [], ...
%                 'timewindow',  [], ...
%                 'win_length',  1, ...
%                 'win_overlap', 50, ...
%                 'sensortypes', 'MEG MAG', ...
%                 'win_std',     0, ...
%                 'edit',        struct(...
%                      'Comment',         'Power', ...
%                      'TimeBands',       [], ...
%                      'Freqs',           [], ...
%                      'ClusterFuncTime', 'none', ...
%                      'Measure',         'power', ...
%                      'Output',          'all', ...
%                      'SaveKernel',      0));
            
        end
        
        % If successful, update subject_array for this subject
        subject_array{pos_subj,3} = 'needs_events';
        save([root_dir '/subject_array.mat'],'subject_array')  
    end
end

clearvars('-except', initialVars{:});

disp 'DONE WITH IMPORTING DATA EEG/MEG FOR cMMN!!!'
disp(datetime)
toc

%% Remove ANY previous projectors and import events (ALWAYS IN CNRL 99 OR 98) 

tic
disp(' ');      
disp('-------------------------');  
disp('REMOVING PROJECTORS AND IMPORTING EVENTS (cMMN)');  
disp(datetime)
disp('-------------------------');     
disp(' '); 

for p = 1:length(participant)
    
    % Check log info about the subject
    pos_subj = find(strcmp({subject_array{:,1}},participant{p}));
    if strcmp(subject_array{pos_subj,3},'needs_events')
   
    for exptask = 1:length(Exp_tasks) 
        
        % Define trigger numbers
        if strcmp(Exp_tasks{exptask},'standardMMN')
            trigger_codes = '1, 2, 3, boundary';
            condition = condition_standardMMN;
        elseif strcmp(Exp_tasks{exptask},'extratoneMMN')
            trigger_codes = '1, 2, 3, 4, boundary';
            condition = condition_extratoneMMN;
        elseif strcmp(Exp_tasks{exptask},'dualruleMMN')
            trigger_codes = '1, 2, 8, boundary';
            condition = condition_dualruleMMN;
        elseif contains(Exp_tasks{exptask},'long_SOA') && any(ismember(participants_long_SOA,participant{p}))
            trigger_codes = '1, 2, 8, boundary';
            condition = condition_long_SOA;
        elseif contains(Exp_tasks{exptask},'long_SOA') && ~any(ismember(participants_long_SOA,participant{p}))
            continue;
        else
            error(['No identifiable task']);
        end
        
        
        folders = dir([root_dir_bs '/data/' participant{p} '/']);
        infolder = find(endsWith({folders.name},'tsss_AMICA') & (...
                    contains({folders.name},[participant{p} '_' Exp_tasks{exptask} '_run']) | ...
                    contains({folders.name},[participant{p} '_' Exp_tasks{exptask} '_b'])));

        if isempty(infolder)
            error(['No cMMN AMICA folders for ' participant{p}]);
        end

        % Remove active projectors
        for i = 1:length(infolder) % because there may be b1 and b2 or several runs
            line = infolder(i);
            sFiles = [root_dir_bs '/data/' participant{p} '/' folders(line).name '/channel_vectorview306_acc1.mat'];  %#ok<*SAGROW>             
            eval(['load ' sFiles])
            variableInfo = who('-file',sFiles);

            disp(' ');      
            disp('-------------------------');  
            disp(['Removing active projectors before importing events: ' participant{p} ' ' Exp_tasks{exptask}]);
            disp(datetime)
            disp(' '); 

            if ~isempty(Projector)
                num_proj = size(Projector,2);
                for npj = 1:num_proj
                    Projector(npj).Status = 0;
                end
            end
            save(sFiles,variableInfo{:});            
        end

        % Now we are going to import events, so reload the subject
        disp(' ');      
        disp('-------------------------');
        disp(['loading participant ' participant{p}]);
        disp(datetime)
        disp(' '); 

        prot_subs = bst_get('ProtocolSubjects');
        current_sub = find(strcmp({prot_subs.Subject.Name}, participant{p}));
        db_reload_conditions(current_sub);

        % Import events now
        folders = dir([root_dir_bs '/data/' participant{p} '/']);
        infolder = find(endsWith({folders.name},'tsss_AMICA') & (...
                    contains({folders.name},[participant{p} '_' Exp_tasks{exptask} '_run']) | ...
                    contains({folders.name},[participant{p} '_' Exp_tasks{exptask} '_b'])));
        if isempty(infolder)
            continue % Go to next subject
        end
        sFiles = {};
        for i= 1:length(infolder)
            line = infolder(i);
            sFiles{i} = [participant{p} '/' folders(line).name '/data_0raw_' folders(line).name(5:end) '.mat'];
        end            
        if isempty(sFiles)
            error(['No AMICA filenames for ' participant{p}])
        end

        % Only if this participant has EEG channels
        if ~strcmp(subject_array{pos_subj,4},'no_EEG_chans')
            
            % Process: CNRL rename EEG Channels
            if strcmp(subject_array{pos_subj,17},'new')
                sFiles = bst_process('CallProcess', 'process_cnrl_rename_eeg', sFiles, [], ...
                'action', 1);  % triux_chanlist.mat
            else
                sFiles = bst_process('CallProcess', 'process_cnrl_rename_eeg_old', sFiles, [], ...
                'action', 2);  % easycap
            end

            % Ensure FCz is healthy
            for i = 1:length(sFiles)
                % Identify position of FCz in this dataset
                filechan = sFiles(i).ChannelFile;
                load([root_dir_bs '/data/' filechan],'Channel');
                pos_fcz = find(strcmp({Channel.Name},'FCz'));
                % Load channel flag file
                fileflag = sFiles(i).FileName;
                load([root_dir_bs '/data/' fileflag]);
                variableInfo = who('-file',[root_dir_bs '/data/' fileflag]);
                ChannelFlag(pos_fcz) = 1;
                F.channelflag(pos_fcz) = 1;
                save([root_dir_bs '/data/' fileflag],variableInfo{:});
            end
            
            % For people with the new caps, M1 MUST be a good channel
            if strcmp(subject_array{pos_subj,17},'new') 
                for i = 1:length(sFiles)
                    % Identify position of M1 in this dataset
                    filechan = sFiles(i).ChannelFile;
                    load([root_dir_bs '/data/' filechan],'Channel');
                    pos_m1 = find(strcmp({Channel.Name},'M1'));
                    % Load channel flag file
                    fileflag = sFiles(i).FileName;
                    load([root_dir_bs '/data/' fileflag]);
                    variableInfo = who('-file',[root_dir_bs '/data/' fileflag]);
                    ChannelFlag(pos_m1) = 1;
                    F.channelflag(pos_m1) = 1;
                    save([root_dir_bs '/data/' fileflag],variableInfo{:});
                    % Do the same with M2 too just in case
                    pos_m2 = find(strcmp({Channel.Name},'M2'));
                    % Load channel flag file
                    fileflag = sFiles(i).FileName;
                    load([root_dir_bs '/data/' fileflag]);
                    variableInfo = who('-file',[root_dir_bs '/data/' fileflag]);
                    ChannelFlag(pos_m2) = 1;
                    F.channelflag(pos_m2) = 1;
                    save([root_dir_bs '/data/' fileflag],variableInfo{:});
                end
                % Reload subject just in case
                disp(' ');      
                disp('-------------------------');
                disp(['loading participant ' participant{p}]);
                disp(datetime)
                disp(' '); 

                prot_subs = bst_get('ProtocolSubjects');
                current_sub = find(strcmp({prot_subs.Subject.Name}, participant{p}));
                db_reload_conditions(current_sub);
            end
            
            % Process: Set channels type
            sFiles = bst_process('CallProcess', 'process_channel_settype', sFiles, [], ...
            'sensortypes', 'HEOG,VEOG', ...
            'newtype',     'EOG');

            % Process: Set channels type
            sFiles = bst_process('CallProcess', 'process_channel_settype', sFiles, [], ...
            'sensortypes', 'ECG', ...
            'newtype',     'ECG');

            % Process: Set channels type: addition on December 2020
            sFiles = bst_process('CallProcess', 'process_channel_settype', sFiles, [], ...
                'sensortypes', 'M1,M2', ...
                'newtype',     'EEG');
            
        
            if reref_option == 1                    
                % Process: Re-reference EEG
                sFiles = bst_process('CallProcess', 'process_eegref', sFiles, [], ...
                'eegref',      ref_EEG, ...
                'sensortypes', 'EEG');  
            end


        end

        if run_matlab_command == 1
        % Process: Run Matlab command
        sFiles = bst_process('CallProcess', 'process_matlab_eval', sFiles, [], ...
            'matlab',      ['% Available variables: Data, TimeVector' 10 '' 10 'TimeVector = TimeVector - TimeVector(1);' 10 ''], ...
            'sensortypes', '');
        end
        
        % Process: Import events file
        sFiles = bst_process('CallProcess', 'process_CNRL_evt_import', sFiles, [], ...
            'evtname', trigger_codes);

        % Process: Convert to simple event
        sFiles = bst_process('CallProcess', 'process_evt_simple', sFiles, [], ...
            'eventname', 'boundary', ...
            'method',    1);  % Keep the start of the events

        % Process: Rename event
        sFiles = bst_process('CallProcess', 'process_evt_rename', sFiles, [], ...
            'src',  'boundary', ...
            'dest', 'boundaryStart');

        % Process: Import events file
        sFiles = bst_process('CallProcess', 'process_CNRL_evt_import', sFiles, [], ...
            'evtname', 'boundary');

        % Process: Convert to simple event
        sFiles = bst_process('CallProcess', 'process_evt_simple', sFiles, [], ...
            'eventname', 'boundary', ...
            'method',    3);  % Keep the end of the events

        % Process: Rename event
        sFiles = bst_process('CallProcess', 'process_evt_rename', sFiles, [], ...
            'src',  'boundary', ...
            'dest', 'boundaryEnd');

        for c = 1:length(condition)
        % Process: Remove simultaneous End
        sFiles = bst_process('CallProcess', 'process_evt_remove_simult', sFiles, [], ...
            'remove', condition{c}, ...
            'target', 'boundaryEnd', ...
            'dt',     0.05, ...
            'rename', 0);

        % Process: Remove simultaneous Start
        sFiles = bst_process('CallProcess', 'process_evt_remove_simult', sFiles, [], ...
            'remove', condition{c}, ...
            'target', 'boundaryStart', ...
            'dt',     0.3, ...
            'rename', 0);
        end

        % Process: scale MEG values (only for files needing it)
        if strcmp(subject_array{pos_subj,4},'scale')
            sFiles = bst_process('CallProcess', 'process_scale', sFiles, [], ...
                'factor',      3, ...
                'sensortypes', 'MEG');
        end
    end
    
    % If successful, update subject_array for this subject
    subject_array{pos_subj,3} = 'needs_filter';
    save([root_dir '/subject_array.mat'],'subject_array')
    
    end
end

clearvars('-except', initialVars{:});
disp 'DONE WITH REMOVING PROJECTORS AND IMPORTING EVENTS (cMMN)!!!'
disp(datetime)
toc

%% Filter

tic
disp(' ');      
disp('-------------------------');  
disp('FILTERING (cMMN)');  
disp(datetime)
disp('-------------------------');     
disp(' '); 

for p = 1:length(participant)
    % Check log info about the subject
    pos_subj = find(strcmp({subject_array{:,1}},participant{p}));
    if strcmp(subject_array{pos_subj,3},'needs_filter')
    
    % Reload subject first
    disp(' ');      
    disp('-------------------------');
    disp(['loading participant ' participant{p}]);
    disp(datetime)
    disp(' '); 

    prot_subs = bst_get('ProtocolSubjects');
    current_sub = find(strcmp({prot_subs.Subject.Name}, participant{p}));
    db_reload_conditions(current_sub);
        
    % Find files    
    folders = dir([root_dir_bs '/data/' participant{p} '/']);
    if strcmp(subject_array{pos_subj,4},'scale')
        infolder = find(contains({folders.name},'raw') & (contains({folders.name},[participant{p} '_standardMMN']) | contains({folders.name},[participant{p} '_extratoneMMN']) | contains({folders.name},[participant{p} '_dualruleMMN'])) & endsWith({folders.name},'matlab_scale'));
    else
        infolder = find(contains({folders.name},'raw') & (contains({folders.name},[participant{p} '_standardMMN']) | contains({folders.name},[participant{p} '_extratoneMMN']) | contains({folders.name},[participant{p} '_dualruleMMN'])) & endsWith({folders.name},'matlab'));
    end
    if isempty(infolder)
        error(['No MMN matlab or matlab_scale folders for ' participant{p}]);
    end  
        
    sFiles = {};
    for i= 1:length(infolder)
        line = infolder(i);
        sFiles{i} = [participant{p} '/' folders(line).name '/data_0raw_' folders(line).name(5:end) '.mat'];
    end

    disp(' ');      
    disp('-------------------------');  
    disp(['Filtering ' participant{p}]);
    disp(datetime)
    disp(' '); 
    
    % If we allow notch filter
    if notch_filter == 1
        % If this participant requires alpha notch
        alpha_indicator = find(strcmp(subject_array(pos_subj,:),'notch_alpha'));
        if ~isempty(alpha_indicator)
            % Process: Notch filter: 60Hz
            sFiles = bst_process('CallProcess', 'process_notch', sFiles, [], ...
                'sensortypes', '', ...
                'freqlist',    notch_alpha, ...
                'cutoffW',     1, ...
                'useold',      0, ...
                'read_all',    0);
        end
        % Add as many others as needed after this
    end
    
    try
    % Process: LLR filter (low pass)
    sFiles = bst_process('CallProcess', 'process_bandpass', sFiles, [], ...
        'sensortypes', '', ...
        'highpass',    LLR_highpass, ... 
        'lowpass',     LLR_lowpass, ... 
        'tranband',    tranband_LLR, ...
        'attenuation', 'strict', ...  % 60dB
        'ver',         '2019', ...  % 2019
        'mirror',      0, ...
        'read_all',    0);
    
    % If successful, update subject_array for this subject
    subject_array{pos_subj,3} = 'needs_stim_resp';
    save([root_dir '/subject_array.mat'],'subject_array')    
    catch
        error(['No filtering performed for ' participant{p}]);
    end
    end  
end

clearvars('-except', initialVars{:});
disp 'DONE WITH FILTERING (cMMN)!!!'
disp(datetime)
toc

%% Combine stim/response

tic
disp(' ');      
disp('-------------------------');  
disp('COMBINING STIM/RESPONSE (cMMN)');  
disp(datetime)
disp('-------------------------');     
disp(' '); 

for p = 1:length(participant)
    % Check log info about the subject
    pos_subj = find(strcmp({subject_array{:,1}},participant{p}));
    if strcmp(subject_array{pos_subj,3},'needs_stim_resp')
    
    % Reload subject first
    disp(' ');      
    disp('-------------------------');
    disp(['loading participant ' participant{p}]);
    disp(datetime)
    disp(' '); 

    prot_subs = bst_get('ProtocolSubjects');
    current_sub = find(strcmp({prot_subs.Subject.Name}, participant{p}));
    db_reload_conditions(current_sub);
       
    for exptask = 1:length(Exp_tasks)
            
        % Find files    
        folders = dir([root_dir_bs '/data/' participant{p} '/']);
        infolder = find(contains({folders.name},'raw') & endsWith({folders.name},'low') & (...
                    contains({folders.name},[participant{p} '_' Exp_tasks{exptask} '_run']) | ...
                    contains({folders.name},[participant{p} '_' Exp_tasks{exptask} '_b'])));
        
        if isempty(infolder)
            error(['No folder names ending in low for ' participant{p} ' ' Exp_tasks{exptask}]);
        end  

        sFiles = {};
        for i= 1:length(infolder)
            line = infolder(i);
            sFiles{i} = [participant{p} '/' folders(line).name '/data_0raw_' folders(line).name(5:end) '.mat'];
        end

        disp(' ');      
        disp('-------------------------');
        disp(['Combining events to epoch data: participant ' participant{p} ' ' Exp_tasks{exptask}]);
        disp(datetime)
        disp(' ');  
        
        % Process: Combine stim/response
        sFiles = bst_process('CallProcess', 'process_evt_combine', sFiles, [], ...
            'combine', ['1_12, ignore, 1, 2' 10 '1_18, ignore, 1, 8'], ...
            'dt',      1);

   
    
    end
    
    % If successful, update subject_array for this subject
    subject_array{pos_subj,3} = 'needs_epoch';
    save([root_dir '/subject_array.mat'],'subject_array') 
    
    end  
end

clearvars('-except', initialVars{:});
disp 'COMBINING STIM/RESPONSE (cMMN)!!!'
disp(datetime)
toc

%% Epoch and amplitude threshold

tic
disp(' ');      
disp('-------------------------');  
disp('EPOCHING (cMMN)');  
disp(datetime)
disp('-------------------------');     
disp(' '); 

for p = 1:length(participant)
    % Check log info about the subject
    pos_subj = find(strcmp({subject_array{:,1}},participant{p}));
    if strcmp(subject_array{pos_subj,3},'needs_epoch')
    
    % Reload subject first
    disp(' ');      
    disp('-------------------------');
    disp(['loading participant ' participant{p}]);
    disp(datetime)
    disp(' '); 

    prot_subs = bst_get('ProtocolSubjects');
    current_sub = find(strcmp({prot_subs.Subject.Name}, participant{p}));
    db_reload_conditions(current_sub);
       
    for exptask = 1:length(Exp_tasks)
    
        % Define trigger numbers
        if strcmp(Exp_tasks{exptask},'standardMMN')
            condition = condition_standardMMN;
            condition_names = condition_names_standardMMN;
            epoch_wave = epoch_wave_standardMMN;
        elseif strcmp(Exp_tasks{exptask},'extratoneMMN')
            condition = condition_extratoneMMN;
            condition_names = condition_names_extratoneMMN;
            epoch_wave = epoch_wave_extratoneMMN;
        elseif strcmp(Exp_tasks{exptask},'dualruleMMN')
            condition = condition_dualruleMMN;
            condition_names = condition_names_dualruleMMN;
            epoch_wave = epoch_wave_dualruleMMN;
        elseif contains(Exp_tasks{exptask},'long_SOA') && any(ismember(participants_long_SOA,participant{p}))
            condition = condition_long_SOA;
            condition_names = condition_names_long_SOA;
            epoch_wave = epoch_wave_long_SOA;
        elseif contains(Exp_tasks{exptask},'long_SOA') && ~any(ismember(participants_long_SOA,participant{p}))
            % If we are in long SOA iteration but this subject does not
            % have it, move onto the next task
            continue;
        else
            error(['No identifiable task']);
        end
        
        % Find files    
        folders = dir([root_dir_bs '/data/' participant{p} '/']);
        infolder = find(contains({folders.name},'raw') & endsWith({folders.name},'low') & (...
                    contains({folders.name},[participant{p} '_' Exp_tasks{exptask} '_run']) | ...
                    contains({folders.name},[participant{p} '_' Exp_tasks{exptask} '_b'])));
        
        if isempty(infolder)
            error(['No folder names ending in low for ' participant{p} ' ' Exp_tasks{exptask}]);
        end  

        sFiles = {};
        for i= 1:length(infolder)
            line = infolder(i);
            sFiles{i} = [participant{p} '/' folders(line).name '/data_0raw_' folders(line).name(5:end) '.mat'];
        end

        disp(' ');      
        disp('-------------------------');
        disp(['Renaming events to epoch data: participant ' participant{p} ' ' Exp_tasks{exptask}]);
        disp(datetime)
        disp(' ');  

        % Rename events (if they are already renamed, nothing will happen)
        for c =1:length(condition)
        sFiles = bst_process('CallProcess', 'process_evt_rename', sFiles, [], ...
            'src',  condition{c}, ...
            'dest', condition_names{c});
        end    

        disp(' ');      
        disp('-------------------------');  
        disp(['Making epochs for ' participant{p} ' ' Exp_tasks{exptask}]);
        disp(datetime)
        disp(' ');  

        % Create string for next step
        if size(condition_names,2) == 2
            string_cond_names = [condition_names{1} ', ' condition_names{2}];
        elseif size(condition_names,2) == 3
            string_cond_names = [condition_names{1} ', ' condition_names{2} ', ' condition_names{3}];
        elseif size(condition_names,2) == 4
            string_cond_names = [condition_names{1} ', ' condition_names{2} ', ' condition_names{3} ', ' condition_names{4}];
        else
            error(['size of condition names does not make sense']);   
        end
        
        % Process: epoch data for normal epochs
        sFiles = bst_process('CallProcess', 'process_import_data_event', sFiles, [], ...
            'subjectname',  participant{p}, ...
            'condition',    '', ...
        ...%    'datafile',     RawFiles, ...
            'eventname',    string_cond_names, ...
            'timewindow',   [], ...
            'epochtime',    epoch_wave, ...
            'createcond',   1, ...
            'ignoreshort',  1, ...
            'channelalign', 0, ...
            'usectfcomp',   0, ...
            'usessp',       1, ...
            'freq',         [], ...
            'baseline',     []); 

        disp(' ');      
        disp('-------------------------');  
        disp(['Baseline correcting epochs for ' participant{p} ' ' Exp_tasks{exptask}]);
        disp(datetime)
        disp(' ');  

        % Process: DC offset correction: [-50ms,1ms]
        sFiles = bst_process('CallProcess', 'process_baseline_norm', sFiles, [], ...
            'baseline',    epoch_baseline, ...
            'sensortypes', '', ...
            'method',      'bl', ...  % DC offset correction:    x_std = x - &mu;
            'overwrite',   1);

        % EEG and BIMODAL cleaning and averages will only take place if 
        % that subject doesn't have a 'no_EEG_chans' label
        if ~strcmp(subject_array{pos_subj,4},'no_EEG_chans')

        disp(' ');      
        disp('-------------------------');  
        disp(['Cleaning epochs for ' participant{p} '(EEG) ' Exp_tasks{exptask}]);
        disp(datetime)
        disp(' '); 

        % Process: Detect bad trials: Absolute threshold EEG
        sFiles_EEG = bst_process('CallProcess', 'process_CNRL_detectbad', sFiles, [], ...
            'timewindow', [], ...
            'meggrad',    [0, 0], ...
            'megmag',     [0, 0], ...
            'eeg',        reject_EEG_absolute, ...
            'ieeg',       [0, 0], ...
            'eog',        [0, 0], ...
            'ecg',        [0, 0], ...
            'rejectmode', 2);  % Reject the entire trial

        % SENSOR AVERAGE EEG  
        disp(' ');      
        disp('-------------------------');  
        disp(['Averaging epochs for ' participant{p} '(EEG) ' Exp_tasks{exptask}]);
        disp(datetime)
        disp(' '); 

        sFiles_EEG = bst_process('CallProcess', 'process_average', sFiles_EEG, [], ...
            'avgtype',         5, ...  % By trial group (folder average)
            'avg_func',        1, ...  % Arithmetic average:  mean(x)
            'weighted',        0, ...
            'keepevents', 0);

        % Process: Add tag
        sFiles_EEG = bst_process('CallProcess', 'process_add_tag', sFiles_EEG, [], ...
            'tag',           'EEG_average', ...
            'output',        2);  % Add to file name (1 to add a tag)

        % Process: Set name
        sFiles_EEG = bst_process('CallProcess', 'process_set_comment', sFiles_EEG, [], ...
            'tag',           'EEG_average', ...
            'isindex',       1);

        % Reset BadTrials variable and reload folder (Std, DevP, DevDur) 
        for c = 1:length(condition_names)
            trials_file = [root_dir_bs '/data/' participant{p} '/' condition_names{c} '/brainstormstudy.mat'];
            if ~exist(trials_file,'file')
                % If there is no file/folder, continue
                warning(['No ' condition_names{c} ' folder for ' participant{p}]);
                continue;
            end  
            eval(['load ' trials_file]); % load them
            variableInfo = who('-file',trials_file);
            % Reset them and save back to original file
            BadTrials = cell([], 1); % exact structure it has when empty
            save(trials_file,variableInfo{:});
            % Reload folder
            [~, iStudies] = bst_get('StudyWithCondition', [participant{p} filesep condition_names{c}]);
            db_reload_studies(iStudies);
        end

        end
    
    end
    
    % If successful, update subject_array for this subject
    subject_array{pos_subj,3} = 'needs_mismatch';
    save([root_dir '/subject_array.mat'],'subject_array') 
    
    end  
end

clearvars('-except', initialVars{:});
disp 'DONE WITH EPOCHING (cMMN)!!!'
disp(datetime)
toc

%% Obtain Mismatch 

tic
disp(' ');      
disp('-------------------------');  
disp('OBTAINING MISMATCH (DEV - STD: cMMN)');  
disp(datetime)
disp('-------------------------');     
disp(' '); 

% Obtain MMN
for p = 1:length(participant)
    % Check log info about the subject
    pos_subj = find(strcmp({subject_array{:,1}},participant{p}));
    if strcmp(subject_array{pos_subj,3},'needs_mismatch')
    
    for exptask = 1:length(Exp_tasks)
       
        % Define trigger numbers
        if strcmp(Exp_tasks{exptask},'standardMMN')
            % condition_mismatch_names_standardMMN = {'Standard','DeviantPitch','DeviantDuration','DeviantPitch-Standard','DeviantDuration-Standard'};
            Standard_conditions = condition_mismatch_names_standardMMN(1);
            Deviant_conditions = condition_mismatch_names_standardMMN(2:3);
            condition_names = condition_names_standardMMN;
        elseif strcmp(Exp_tasks{exptask},'extratoneMMN')
            % condition_mismatch_names_extratoneMMN = {'Standard1','Standard2','Standard3','DeviantExtra','DeviantExtra-Standard1','DeviantExtra-Standard2','DeviantExtra-Standard3'};
            Standard_conditions = condition_mismatch_names_extratoneMMN(1:3);
            Deviant_conditions = condition_mismatch_names_extratoneMMN(4);
            condition_names = condition_names_extratoneMMN;
        elseif strcmp(Exp_tasks{exptask},'dualruleMMN')
            % condition_mismatch_names_dualruleMMN = {'LeftPitchA','RightPitchB','DeviantDoubleRule','DeviantDoubleRule-LeftPitchA','DeviantDoubleRule-RightPitchB'};
            Standard_conditions = condition_mismatch_names_dualruleMMN(1);
            Deviant_conditions = condition_mismatch_names_dualruleMMN(2);
            condition_names = condition_names_dualruleMMN;
        elseif contains(Exp_tasks{exptask},'long_SOA') && any(ismember(participants_long_SOA,participant{p}))
            % condition_mismatch_names_long_SOA = {'LeftPitchA_long_SOA','RightPitchB_long_SOA','DeviantDoubleRule_long_SOA','DeviantDoubleRule_long_SOA-LeftPitchA_long_SOA','DeviantDoubleRule_long_SOA-RightPitchB_long_SOA'};
            Standard_conditions = condition_mismatch_names_long_SOA(1:2);
            Deviant_conditions = condition_mismatch_names_long_SOA(3);
            condition_names = condition_names_long_SOA;
        elseif contains(Exp_tasks{exptask},'long_SOA') && ~any(ismember(participants_long_SOA,participant{p}))
            % If we are in long SOA iteration but this subject does not
            % have it, move onto the next task
            continue;
        else
            error(['No identifiable task']);
        end
        
        % Reload subject epoch folders first
        disp(' ');      
        disp('-------------------------');
        disp(['loading epoch folders for participant ' participant{p}]);
        disp(datetime)
        disp(' '); 
        for c = 1:length(condition_names)
            [~, iStudies] = bst_get('StudyWithCondition', [participant{p} filesep condition_names{c}]);
            db_reload_studies(iStudies);
        end

        % Do once for each modality (EEG, MEG, BIMODAL)
        for mode = 1:length(modality_data)
            
        % For each Deviant - Standard combination
        for stdc = 1:length(Standard_conditions)
            % Define standard sFile
            files = dir([root_dir_bs '/data/' participant{p} '/' Standard_conditions{stdc}]);
            if isempty(files)
                error(['No ' Standard_conditions{stdc} ' files for ' participant{p}]);
            end
            infolder = find(endsWith({files.name},[modality_data{mode} '_average.mat'])); 
            if isempty(infolder)
               % In case, for instance, there is no EEG data for this subject
               warning(['No ' Standard_conditions{stdc} ' ' modality_data{mode} ' average for ' participant{p}]);
               continue;
            end  
            if length(infolder) > 1
                error(['More than one ' Standard_conditions{stdc} ' ' modality_data{mode} ' average for ' participant{p}]);
            end
            filename = [participant{p} '/' Standard_conditions{stdc} '/' files(infolder).name];
            STD_sFiles = {}; STD_sFiles{1} = filename; 
            % Now, for each deviant for this standard file
            for devc = 1:length(Deviant_conditions)
                files = dir([root_dir_bs '/data/' participant{p} '/' Deviant_conditions{devc}]);
                if isempty(files)
                    error(['No ' Deviant_conditions{devc} ' files for ' participant{p}]);
                end
                infolder = find(endsWith({files.name},[modality_data{mode} '_average.mat'])); 
                if isempty(infolder)
                    % In case, for instance, there is no EEG data for this subject
                    warning(['No ' Deviant_conditions{devc} ' ' modality_data{mode} ' average for ' participant{p}]);
                    continue;
                end  
                if length(infolder) > 1
                    error(['More than one ' Deviant_conditions{devc} ' ' modality_data{mode} ' average for ' participant{p}]);
                end
                filename = [participant{p} '/' Deviant_conditions{devc} '/' files(infolder).name];
                DEV_sFiles = {}; DEV_sFiles{1} = filename;
                
                disp(' ');      
                disp('-------------------------');  
                disp(['Calculating ' Deviant_conditions{devc} ' - ' Standard_conditions{stdc} ' MMN for ' participant{p} '(' modality_data{mode} ')']);
                disp(datetime)
                disp(' '); 

                % Process: Difference: A-B
                sFiles_mismatch = bst_process('CallProcess', 'process_diff_ab', DEV_sFiles, STD_sFiles);

                % Give outcome file a new name
                sFiles_mismatch = bst_process('CallProcess', 'process_add_tag', sFiles_mismatch, [], ...
                'tag',           ['_' modality_data{mode} '_' Deviant_conditions{devc} '-' Standard_conditions{stdc}  '_MMN'], ...
                'output',        2);  % Add to file name (1 to add a tag)

                % Process: Set name
                sFiles_mismatch = bst_process('CallProcess', 'process_set_comment', sFiles_mismatch, [], ...
                'tag',           [modality_data{mode} '_' Deviant_conditions{devc} '-' Standard_conditions{stdc}  '_MMN'], ...
                'isindex',       1);
                
            end
        end
        end
    end
    
    % If successful, update subject_array for this subject
    subject_array{pos_subj,3} = 'needs_forward';
    save([root_dir '/subject_array.mat'],'subject_array')    
    end  
end

clearvars('-except', initialVars{:});
disp 'DONE WITH OBTAINING MISMATCH (DEV - STD: cMMN)!!!'
disp(datetime)
toc

%% GAVR scalp level (Brainstorm) 

% May want to reload the entire protocol first to fix errors

tic
disp(' ');      
disp('-------------------------');  
disp('GAVR SCALP LEVEL (cMMN)');  
disp(datetime)
disp('-------------------------');     
disp(' '); 

% First, reload all folders for each subject
for p = 1:length(participant)    

    % This step is important as new epoch folders for MMN will be created
    % and won't be read by brainstorm unless we reload the whole subject
    % functional data
    disp(' ');      
    disp('-------------------------');
    disp(['loading participant ' participant{p}]);
    disp(datetime)
    disp(' ');
    
    prot_subs = bst_get('ProtocolSubjects');
    current_sub = find(strcmp({prot_subs.Subject.Name}, participant{p}));
    db_reload_conditions(current_sub);
end

% Delete duplicates/keep the most recent in epoch folders
for p = 1:length(participant)
    % Check log info about the subject
    pos_subj = find(strcmp({subject_array{:,1}},participant{p}));
    
    % If subject does not have any of these markers, don't include it in
    % the sweep count (because it may be a bad subject, etc)
    if ~strcmp(subject_array{pos_subj,3},'needs_forward') && ~strcmp(subject_array{pos_subj,3},'needs_sources') && ~strcmp(subject_array{pos_subj,3},'sources_finished')
        continue;
    end
    for exptask = 1:length(Exp_tasks)
        
        % Define trigger numbers
        if strcmp(Exp_tasks{exptask},'standardMMN')
            condition_mismatch_names = condition_mismatch_names_standardMMN;
        elseif strcmp(Exp_tasks{exptask},'extratoneMMN')
            condition_mismatch_names = condition_mismatch_names_extratoneMMN;
        elseif strcmp(Exp_tasks{exptask},'dualruleMMN')
            condition_mismatch_names = condition_mismatch_names_dualruleMMN;
        elseif contains(Exp_tasks{exptask},'long_SOA') && any(ismember(participants_long_SOA,participant{p}))
            condition_mismatch_names = condition_mismatch_names_long_SOA;
        elseif contains(Exp_tasks{exptask},'long_SOA') && ~any(ismember(participants_long_SOA,participant{p}))
            continue;
        else
            error('No identifiable task');
        end
        
        % Only for the mismatch folders
        for c = 1:length(condition_mismatch_names)
        % Find folder
        files = dir([root_dir_bs '/data/' participant{p} '/' condition_mismatch_names{c}]);
        if isempty(files)
            error(['No ' condition_mismatch_names{c} ' files for ' participant{p}]);
        end
        % Do once for each modality (EEG, MEG, BIMODAL)
        for mode = 1:length(modality_data)
        % Retrieve average/s files
        infolder = find(contains({files.name},[modality_data{mode}]) & (...
        endsWith({files.name}, '_MMN.mat') | endsWith({files.name}, '_average.mat')));

        if isempty(infolder)
           % In case, for instance, there is no EEG data for this subject
           warning(['No ' condition_mismatch_names{c} ' ' modality_data{mode} ' average for ' participant{p}]);
           continue; 
        elseif length(infolder) == 1 % there are no duplicates to delete
            continue;
        elseif length(infolder) > 1 % there are duplicates to delete
            for i = 1:length(infolder)-1 % all except last one in list (older)
                line = infolder(i);
                delete([root_dir_bs '/data/' participant{p} '/' condition_mismatch_names{c} '/' files(line).name]); % always has this format
            end
        end
        end
        end 
    end
end

% GAVR for each modality and condition (separately for FE/C)
for mode = 1:length(modality_data)

    for exptask = 1:length(Exp_tasks)

    % Define trigger numbers
    if strcmp(Exp_tasks{exptask},'standardMMN')
        condition_mismatch_names = condition_mismatch_names_standardMMN;
    elseif strcmp(Exp_tasks{exptask},'extratoneMMN')
        condition_mismatch_names = condition_mismatch_names_extratoneMMN;
    elseif strcmp(Exp_tasks{exptask},'dualruleMMN')
        condition_mismatch_names = condition_mismatch_names_dualruleMMN;
    elseif contains(Exp_tasks{exptask},'long_SOA') && any(ismember(participants_long_SOA,participant{p}))
        condition_mismatch_names = condition_mismatch_names_long_SOA;
    elseif contains(Exp_tasks{exptask},'long_SOA') && ~any(ismember(participants_long_SOA,participant{p}))
        continue;
    else
        error(['No identifiable task']);
    end

    for c = 1:length(condition_mismatch_names)
        for pg = 1:length(participant_group)
        sFiles = {};
        for p = 1:length(participant)
        % Identify patients/controls
        pos_subj = find(strcmp({subject_array{:,1}},participant{p}));

        % If subject does not have any of these markers, don't include (because it may be a bad subject, etc)
        if ~strcmp(subject_array{pos_subj,3},'needs_forward') && ~strcmp(subject_array{pos_subj,3},'needs_sources') && ~strcmp(subject_array{pos_subj,3},'sources_finished')
            continue;
        end

        if ~strcmp(subject_array{pos_subj,2},participant_group{pg})
            continue; % so only include participants that correspond to the group
        end

        % Check if we should include it
        subject_row = {subject_array{pos_subj,:}}; 
        subject_row = subject_row(~cellfun('isempty', subject_row')); % to avoid empty cells in contains
        bad_signal = find(contains(subject_row,['bad_' modality_data{mode}]));
        if strcmp(modality_data{mode},'BIMODAL') % In case of bimodal, check if any EEG or MEG is bad
            bad_signal_EEG = find(contains(subject_row,'bad_EEG'));
            bad_signal_MEG = find(contains(subject_row,'bad_MEG'));
            % If either one of these is bad, mark bimodal as bad
            if ~isempty(bad_signal_EEG) || ~isempty(bad_signal_MEG)
                bad_signal = 1;
            end
        end
        if ~isempty(bad_signal)
            warning([participant{p} ' not included in ' modality_data{mode} ' scalp average']);
            continue; % to next subject
        end     

        files = dir([root_dir_bs '/data/' participant{p} '/' condition_mismatch_names{c}]);
        if isempty(files)
            % If both EEG and MEG were bad there will be no mismatch folder at all
            warning(['No ' condition_mismatch_names{c} ' files for ' participant{p}]);
            continue;
        end
        
        infolder = find(contains({files.name},[modality_data{mode}]) & (...
        endsWith({files.name}, '_MMN.mat') | endsWith({files.name}, '_average.mat')));
        
        if isempty(infolder)
           % In case, for instance, there is no EEG data for this subject
           warning(['No ' condition_mismatch_names{c} ' ' modality_data{mode} ' average for ' participant{p}]);
           continue;
        end  
        if length(infolder) > 1
            error(['More than one ' condition_mismatch_names{c} ' ' modality_data{mode} ' average for ' participant{p}]);
        end
        sFiles{p} = [participant{p} '/' condition_mismatch_names{c} '/' files(infolder).name];
        end

        sFiles = sFiles(~cellfun('isempty', sFiles')); % to avoid empty cells
        if isempty(sFiles)
            error(['No files to perform GAVR for ' condition_mismatch_names{c} ' ' modality_data{mode}]);
        end

        gavr_n = num2str(length(sFiles));

        disp(' ');
        disp('-------------------------');  
        disp(['GAVR scalp data for ' condition_mismatch_names{c} ' ' modality_data{mode} ' ' participant_group{pg}]);
        disp(datetime)
        disp(' ');

        % If stated, find and delete any previous GAVR SENSOR data
        if delete_previous_file == 1
            % check if there is already GAVR source in Group analysis folder
            folders_delete = dir([root_dir_bs '/data/Group_analysis/' condition_mismatch_names{c}]);
            infolder_delete = find(contains({folders_delete.name},'GAVR_SENSOR_')...
            & endsWith({folders_delete.name}, [modality_data{mode} '_' condition_mismatch_names{c} '_' participant_group{pg} '_n' gavr_n '.mat']));
            if ~isempty(infolder_delete) % file exists, therefore delete it
               delete([root_dir_bs '/data/Group_analysis/' condition_mismatch_names{c} '/' folders_delete(infolder_delete).name]);
            end
        end

        % Average using default function (because we are using vertices, not channels)
        sFiles = bst_process('CallProcess', 'process_average', sFiles, [], ...
            'avgtype',         1, ...  % Everything
            'avg_func',        1, ...  % Arithmetic average:  mean(x)
            'weighted',        0, ...
            'scalenormalized', 0);

        % Process: Add tag
        sFiles = bst_process('CallProcess', 'process_add_tag', sFiles, [], ...
            'tag',           ['GAVR_SENSOR_' modality_data{mode} '_' condition_mismatch_names{c} '_' participant_group{pg} '_n' gavr_n], ...
            'output',        2);  % Add to file name (1 to add a tag)

        % Process: Set name
        sFiles = bst_process('CallProcess', 'process_set_comment', sFiles, [], ...
            'tag',           ['GAVR_SENSOR_' modality_data{mode} '_' condition_mismatch_names{c} '_' participant_group{pg} '_n' gavr_n], ...
            'isindex',       1);
        end
    end   
    end
end

clearvars('-except', initialVars{:});
disp 'DONE WITH GAVR SCALP LEVEL (cMMN)!!!'
disp(datetime)
toc

%% Exract scalp data out of brainstorm and GAVR outside

tic
disp(' ');      
disp('-------------------------');  
disp('EXTRACTING SCALP VALUES OUT OF BRAINSTORM (cMMN)');  
disp(datetime)
disp('-------------------------');     
disp(' '); 

% In this case, modality data will always be EEG and MEG only, since it's scalp
modality_data = {'EEG'};

% Select which channels/clusters are going to be used from matrix
load([root_dir '/Scripts/Standard_EEG_channel_order.mat']);
% EEG channels
if size(choice_channel_EEG,2) > 1 % cluster was selected
    channel_pos_EEG = [];
    for i = 1:size(choice_channel_EEG,2) % Find all the indices
        channel_pos_EEG(i) = find(strcmp(std_chan_order, choice_channel_EEG)); 
    end
else % single channel
   channel_pos_EEG = find(strcmp(std_chan_order, choice_channel_EEG)); 
end
load([root_dir '/Scripts/Standard_MEG_channel_order.mat']);
% MEG sensors
if size(choice_channel_MEG,2) > 1 % cluster was selected
    channel_pos_MEG = [];
    for i = 1:size(choice_channel_MEG,2) % Find all the indices
        channel_pos_MEG(i) = find(strcmp(std_chan_order, choice_channel_MEG)); 
    end
else % single sensor
   channel_pos_MEG = find(strcmp(std_chan_order, choice_channel_MEG)); 
end

% Extract individual responses of each participant and condition
for p = 1:length(participant)
    for mode = 1:length(modality_data)
        pos_subj = find(strcmp({subject_array{:,1}},participant{p}));

        % If subject does not have any of these markers, don't include it in
        % the sweep count (because it may be a bad subject, etc)
        if ~strcmp(subject_array{pos_subj,3},'needs_forward') && ~strcmp(subject_array{pos_subj,3},'needs_sources') && ~strcmp(subject_array{pos_subj,3},'sources_finished')
            continue;
        end
        
        % Check if we should include it
        subject_row = {subject_array{pos_subj,:}}; 
        subject_row = subject_row(~cellfun('isempty', subject_row')); % to avoid empty cells in contains
        bad_signal = find(contains(subject_row,['bad_' modality_data{mode}]));
        if ~isempty(bad_signal)
            warning([participant{p} ' scalp data not extracted for ' modality_data{mode}]);
            continue; % to next subject
        end    

        if ~exist([root_dir '/Scalp/' participant{p}], 'dir')
            mkdir([root_dir '/Scalp/'], [participant{p}]);
        end
        
        for exptask = 1:length(Exp_tasks)
        
        % Define trigger numbers
        if strcmp(Exp_tasks{exptask},'standardMMN')
            condition_mismatch_names = condition_mismatch_names_standardMMN;
            condition_short_labels = condition_short_labels_standardMMN;
        elseif strcmp(Exp_tasks{exptask},'extratoneMMN')
            condition_mismatch_names = condition_mismatch_names_extratoneMMN;
            condition_short_labels = condition_short_labels_extratoneMMN;
        elseif strcmp(Exp_tasks{exptask},'dualruleMMN')
            condition_mismatch_names = condition_mismatch_names_dualruleMMN;
            condition_short_labels = condition_short_labels_dualruleMMN;
        elseif contains(Exp_tasks{exptask},'long_SOA') && any(ismember(participants_long_SOA,participant{p}))
            condition_mismatch_names = condition_mismatch_names_long_SOA;
            condition_short_labels = condition_short_labels_long_SOA;
        elseif contains(Exp_tasks{exptask},'long_SOA') && ~any(ismember(participants_long_SOA,participant{p}))
            continue;
        else
            error(['No identifiable task']);
        end
        
            for c = 1:length(condition_mismatch_names)  

            files = dir([root_dir_bs '/data/' participant{p} '/' condition_mismatch_names{c}]);
            if isempty(files)
                % If both EEG and MEG were bad there will be no mismatch folder at all
                warning(['No ' condition_mismatch_names{c} ' files for ' participant{p}]);
                continue;
            end
            
            % Retrieve average/s files
            infolder = find(contains({files.name},[modality_data{mode}]) & (...
            endsWith({files.name}, '_MMN.mat') | endsWith({files.name}, '_average.mat')));
            
            if isempty(infolder)
               % In case, for instance, there is no EEG data for this subject
               warning(['No ' condition_mismatch_names{c} ' ' modality_data{mode} ' average for ' participant{p}]);
               continue;
            end  
            if length(infolder) > 1
                error(['More than one ' condition_mismatch_names{c} ' ' modality_data{mode} ' average for ' participant{p}]);
            end
            file_name = [root_dir_bs '/data/' participant{p} '/' condition_mismatch_names{c} '/' files(infolder).name];
            channel_name = [root_dir_bs '/data/' participant{p} '/' condition_mismatch_names{c} '/channel_vectorview306_acc1.mat'];
            
            % This way so that we can preserve the structure for the plotting script
            eval(['load ' file_name])
            % Load also channel file to check indices (since some subjects only have MEG channels
            eval(['load ' channel_name])
            chan_index = find(contains({Channel.Type},modality_data{mode}));
            chan_list = {Channel(chan_index).Name};
            F = F(chan_index,:);
            % Reorder channel signal to match standard (since some subjects
            % were recorded with 60 scalp electrodes and some others with
            % 63, and in different order
            load([root_dir '/Scripts/Standard_' modality_data{mode} '_channel_order.mat']); % std_chan_order
            ordered_F = [];
            for sco = 1:length(std_chan_order)
                pos_in_chan_list = find(strcmp(chan_list,std_chan_order{sco}));
                ordered_F(sco,1:length(F(pos_in_chan_list,:))) = F(pos_in_chan_list,:);
            end
            F = ordered_F;
            save([root_dir '/Scalp/' participant{p} '/' modality_data{mode} '_' condition_mismatch_names{c} '.mat'], 'F');     
            end
        end  
    end
end

% GAVR outside brainstorm (and obtain Matrix)
for mode = 1:length(modality_data)
    
    for exptask = 1:length(Exp_tasks)
        
        % Define trigger numbers
        if strcmp(Exp_tasks{exptask},'standardMMN')
            condition_mismatch_names = condition_mismatch_names_standardMMN;
            condition_short_labels = condition_short_labels_standardMMN;
        elseif strcmp(Exp_tasks{exptask},'extratoneMMN')
            condition_mismatch_names = condition_mismatch_names_extratoneMMN;
            condition_short_labels = condition_short_labels_extratoneMMN;
        elseif strcmp(Exp_tasks{exptask},'dualruleMMN')
            condition_mismatch_names = condition_mismatch_names_dualruleMMN;
            condition_short_labels = condition_short_labels_dualruleMMN;
        elseif contains(Exp_tasks{exptask},'long_SOA')
            condition_mismatch_names = condition_mismatch_names_long_SOA;
            condition_short_labels = condition_short_labels_long_SOA;
        else
            error(['No identifiable task']);
        end
        
        for c = 1:length(condition_mismatch_names)
            for pg = 1:length(participant_group)   
            file_matrix = {};
            for p = 1:length(participant)
                pos_subj = find(strcmp({subject_array{:,1}},participant{p}));

                % If subject does not have any of these markers, don't include it in
                % the sweep count (because it may be a bad subject, etc)
                if ~strcmp(subject_array{pos_subj,3},'needs_forward') && ~strcmp(subject_array{pos_subj,3},'needs_sources') && ~strcmp(subject_array{pos_subj,3},'sources_finished')
                    continue;
                end

                % Identify patients/controls
                if ~strcmp(subject_array{pos_subj,2},participant_group{pg})
                    continue; % so only include participants that correspond to the group
                end

                % Again, check if we should include it
                subject_row = {subject_array{pos_subj,:}}; 
                subject_row = subject_row(~cellfun('isempty', subject_row')); % to avoid empty cells in contains
                bad_signal = find(contains(subject_row,['bad_' modality_data{mode}]));
                if ~isempty(bad_signal)
                    warning([participant{p} ' scalp data not averaged for ' modality_data{mode}]);
                    continue; % to next subject
                end  
                
                if contains(Exp_tasks{exptask},'long_SOA') && ~any(ismember(participants_long_SOA,participant{p}))
                    continue;
                end

                load([root_dir '/Scalp/' participant{p} '/' modality_data{mode} '_' condition_mismatch_names{c} '.mat']);
                if strcmp(modality_data{mode},'EEG')
                    F = F*1e6; % Transform from Volts to microVolts
                elseif strcmp(modality_data{mode},'MEG') % MEG
                    F = F*1e15; % Transform from Tesla to Femtotesla
                end
                if size(F,1) > 306 && strcmp(modality_data{mode},'MEG')
                    disp(['WARNING!! MEG DATA HAS MORE THAN 306 CHANNELS in ' participant{p}]);
                    F(307:end,:)=[];
                elseif size(F,1) > 62 && strcmp(modality_data{mode},'EEG')
                    disp(['WARNING!! EEG DATA HAS MORE THAN 62 CHANNELS in ' participant{p}]);                   
                    F(63:end,:)=[];
                end
                file_matrix{p} = F;
            end

            % Delete empty subjects
            file_matrix = file_matrix(~cellfun('isempty', file_matrix')); % to avoid empty cells in contains
            % Generate files
            Matrix = cat(3,file_matrix{:});
            gavr = mean(Matrix,3);
            STD = squeeze(std(Matrix,0,3));
            num_valid_subjects = size(Matrix,3);
            STD_ERR = STD/sqrt(num_valid_subjects);

            % Save files
            if ~exist([root_dir '/Scalp/' gavr_name], 'dir')
                mkdir([root_dir '/Scalp/' gavr_name], 'gavr');
                mkdir([root_dir '/Scalp/' gavr_name], 'std_dev');
                mkdir([root_dir '/Scalp/' gavr_name], 'std_err');
            end
            F = gavr;
            eval(['save(''' root_dir '/Scalp/' gavr_name '/gavr/' modality_data{mode} '_' condition_mismatch_names{c} '_' participant_group{pg} ''', ''F'');'])
            F = STD;
            eval(['save(''' root_dir '/Scalp/' gavr_name '/std_dev/' modality_data{mode} '_' condition_mismatch_names{c} '_' participant_group{pg} ''', ''F'');'])
            F = STD_ERR;
            eval(['save(''' root_dir '/Scalp/' gavr_name  '/std_err/' modality_data{mode} '_' condition_mismatch_names{c} '_' participant_group{pg} ''', ''F'');'])
            % Save the matrix for peak average extraction too
            eval(['save(''' root_dir '/Scalp/' gavr_name '/Matrix_' modality_data{mode} '_' condition_mismatch_names{c} '_' participant_group{pg} ''', ''Matrix'');'])  
           end
        end    
    end
end

% Obtain peak matrices
for mode = 1:length(modality_data)
    eval(['chan_num = channel_pos_' modality_data{mode} ';'])
for pg = 1:length(participant_group)   
for exptask = 1:length(Exp_tasks)
        
    % Define trigger numbers
    if strcmp(Exp_tasks{exptask},'standardMMN')
        condition_mismatch_names = condition_mismatch_names_standardMMN_peaks;
        condition_short_labels = condition_short_labels_standardMMN_peaks;
        % Define time samples, which we are gonna need
        plot_baseline = epoch_wave_standardMMN(1)*1000;
        plot_post = epoch_wave_standardMMN(2)*1000;
        s_r = 1000;
        time_samples=linspace(plot_baseline,plot_post,((((plot_baseline*(-1)) + plot_post)/1000)*s_r) +1);
        Peaks_to_extract = Peaks_to_extract_standardMMN;
    elseif strcmp(Exp_tasks{exptask},'extratoneMMN')
        condition_mismatch_names = condition_mismatch_names_extratoneMMN;
        condition_short_labels = condition_short_labels_extratoneMMN;
        % Define time samples, which we are gonna need
        plot_baseline = epoch_wave_extratoneMMN(1)*1000;
        plot_post = epoch_wave_extratoneMMN(2)*1000;
        s_r = 1000;
        time_samples=linspace(plot_baseline,plot_post,((((plot_baseline*(-1)) + plot_post)/1000)*s_r) +1);
        Peaks_to_extract = Peaks_to_extract_extratoneMMN;
    elseif strcmp(Exp_tasks{exptask},'dualruleMMN')
        condition_mismatch_names = condition_mismatch_names_dualruleMMN;
        condition_short_labels = condition_short_labels_dualruleMMN;
        % Define time samples, which we are gonna need
        plot_baseline = epoch_wave_dualruleMMN(1)*1000;
        plot_post = epoch_wave_dualruleMMN(2)*1000;
        s_r = 1000;
        time_samples=linspace(plot_baseline,plot_post,((((plot_baseline*(-1)) + plot_post)/1000)*s_r) +1);
        Peaks_to_extract = Peaks_to_extract_dualruleMMN;
    elseif contains(Exp_tasks{exptask},'long_SOA')
        condition_mismatch_names = condition_mismatch_names_long_SOA;
        condition_short_labels = condition_short_labels_long_SOA;
        % Define time samples, which we are gonna need
        plot_baseline = epoch_wave_long_SOA(1)*1000;
        plot_post = epoch_wave_long_SOA(2)*1000;
        s_r = 1000;
        time_samples=linspace(plot_baseline,plot_post,((((plot_baseline*(-1)) + plot_post)/1000)*s_r) +1);
        Peaks_to_extract = Peaks_to_extract_dualruleMMN;
    else
        error(['No identifiable task']);
    end

    % Obtain amplitude values from Matrix, selecting channel (conditions)
    for c = 1:length(condition_mismatch_names)
         
        % Since standard is divided in two, do this
        if strcmp(condition_mismatch_names{c},'Standard_PD')
            load([root_dir '/Scalp/' gavr_name '/Matrix_' modality_data{mode} '_Standard_' participant_group{pg} '.mat'])
        elseif strcmp(condition_mismatch_names{c},'Standard_DD')
            load([root_dir '/Scalp/' gavr_name '/Matrix_' modality_data{mode} '_Standard_' participant_group{pg} '.mat'])
        elseif strcmp(condition_mismatch_names{c},'DeviantPitch')
            load([root_dir '/Scalp/' gavr_name '/Matrix_' modality_data{mode} '_DeviantPitch_' participant_group{pg} '.mat'])
        elseif strcmp(condition_mismatch_names{c},'DeviantDuration')
            load([root_dir '/Scalp/' gavr_name '/Matrix_' modality_data{mode} '_DeviantDuration_' participant_group{pg} '.mat'])
        elseif strcmp(condition_mismatch_names{c},'DeviantPitch-Standard')
            load([root_dir '/Scalp/' gavr_name '/Matrix_' modality_data{mode} '_DeviantPitch-Standard_' participant_group{pg} '.mat'])
        elseif strcmp(condition_mismatch_names{c},'DeviantDuration-Standard')
            load([root_dir '/Scalp/' gavr_name '/Matrix_' modality_data{mode} '_DeviantDuration-Standard_' participant_group{pg} '.mat'])
        else
            % Load matrix of data
            load([root_dir '/Scalp/' gavr_name '/Matrix_' modality_data{mode} '_' condition_mismatch_names{c} '_' participant_group{pg} '.mat'])
        end
        
        % Ensure correct number of channels
        if strcmp(modality_data{mode}, 'EEG')
            if size(Matrix,1) ~= 62
                error('watch out, matrix does not have 62 channels')
            end
        elseif strcmp(modality_data{mode}, 'MEG')
            if size(Matrix,1) ~= 306
                error('watch out, matrix does not have 306 channels')
            end
        end

        % Retrieve values from channels
        Amplitudes = Matrix(chan_num,:,:);
        % if it's a cluster of channels, average first
        if size(Amplitudes,1) ~= 1
           Amplitudes = squeeze(mean(Amplitudes,1));
        else
           Amplitudes = squeeze(Amplitudes);
        end

        for pe = 1:length(Peaks_to_extract)   
            
            eval(['time_window = time_window_scalp_' Peaks_to_extract{pe} '_' condition_short_labels{c} '_' modality_data{mode} ';'])

            [~,closestIndex] = min(abs(time_samples-time_window(1)));
            init_time = closestIndex;
            [~,closestIndex] = min(abs(time_samples-time_window(2)));
            end_time = closestIndex;
            F = mean(Amplitudes(init_time:end_time,:),1);

            if ~exist([root_dir '/Scalp/' gavr_name '/Peaks'],'dir')
             mkdir([root_dir '/Scalp/' gavr_name '/'], 'Peaks');
            end

            % size(choice_channel_EEG,2) > 1
            if length(chan_num) >1 % Group of electrodes/sensors
                save([root_dir '/Scalp/' gavr_name '/Peaks/' modality_data{mode} '_' condition_short_labels{c} '_' Peaks_to_extract{pe} '_' participant_group{pg} '_cluster.mat'],'F');
            elseif length(chan_num) == 1
               if strcmp(modality_data{mode},'EEG')
                   save([root_dir '/Scalp/' gavr_name '/Peaks/' modality_data{mode} '_' condition_short_labels{c} '_' Peaks_to_extract{pe} '_' participant_group{pg} '_' choice_channel_EEG{:} '.mat'],'F');
               else
                   save([root_dir '/Scalp/' gavr_name '/Peaks/' modality_data{mode} '_' condition_short_labels{c} '_' Peaks_to_extract{pe} '_' participant_group{pg} '_' choice_channel_MEG{:} '.mat'],'F');
               end
            end

        end
    end
end
end
end

% SPSS Matrix
for mode = 1:length(modality_data)
eval(['chan_num = choice_channel_' modality_data{mode} ';'])
temporal = []; iter = 1;
header_cMMN = {}; pos_head = 1;
for pg = 1:length(participant_group)
for exptask = 1:length(Exp_tasks)
        
    % Define trigger numbers
    if strcmp(Exp_tasks{exptask},'standardMMN')
        condition_mismatch_names = condition_mismatch_names_standardMMN_peaks;
        condition_short_labels = condition_short_labels_standardMMN_peaks;
        % Define time samples, which we are gonna need
        plot_baseline = epoch_wave_standardMMN(1)*1000;
        plot_post = epoch_wave_standardMMN(2)*1000;
        s_r = 1000;
        time_samples=linspace(plot_baseline,plot_post,((((plot_baseline*(-1)) + plot_post)/1000)*s_r) +1);
        Peaks_to_extract = Peaks_to_extract_standardMMN;
    elseif strcmp(Exp_tasks{exptask},'extratoneMMN')
        condition_mismatch_names = condition_mismatch_names_extratoneMMN;
        condition_short_labels = condition_short_labels_extratoneMMN;
        % Define time samples, which we are gonna need
        plot_baseline = epoch_wave_extratoneMMN(1)*1000;
        plot_post = epoch_wave_extratoneMMN(2)*1000;
        s_r = 1000;
        time_samples=linspace(plot_baseline,plot_post,((((plot_baseline*(-1)) + plot_post)/1000)*s_r) +1);
        Peaks_to_extract = Peaks_to_extract_extratoneMMN;
    elseif strcmp(Exp_tasks{exptask},'dualruleMMN')
        condition_mismatch_names = condition_mismatch_names_dualruleMMN;
        condition_short_labels = condition_short_labels_dualruleMMN;
        % Define time samples, which we are gonna need
        plot_baseline = epoch_wave_dualruleMMN(1)*1000;
        plot_post = epoch_wave_dualruleMMN(2)*1000;
        s_r = 1000;
        time_samples=linspace(plot_baseline,plot_post,((((plot_baseline*(-1)) + plot_post)/1000)*s_r) +1);
        Peaks_to_extract = Peaks_to_extract_dualruleMMN;
    elseif contains(Exp_tasks{exptask},'long_SOA')
        condition_mismatch_names = condition_mismatch_names_long_SOA;
        condition_short_labels = condition_short_labels_long_SOA;
        % Define time samples, which we are gonna need
        plot_baseline = epoch_wave_long_SOA(1)*1000;
        plot_post = epoch_wave_long_SOA(2)*1000;
        s_r = 1000;
        time_samples=linspace(plot_baseline,plot_post,((((plot_baseline*(-1)) + plot_post)/1000)*s_r) +1);
        Peaks_to_extract = Peaks_to_extract_dualruleMMN;
    else
        error(['No identifiable task']);
    end
    
    for c = 1:length(condition_mismatch_names)
    for pe = 1:length(Peaks_to_extract)      
        if size(chan_num,2) > 1 % Cluster
            load ([root_dir '/Scalp/' gavr_name '/Peaks/' modality_data{mode} '_' condition_short_labels{c} '_' Peaks_to_extract{pe} '_' participant_group{pg} '_cluster.mat']);
        else % Single channel
            load ([root_dir '/Scalp/' gavr_name '/Peaks/' modality_data{mode} '_' condition_short_labels{c} '_' Peaks_to_extract{pe} '_' participant_group{pg} '_' chan_num{:} '.mat']);
        end
        temporal(1:length(F),iter) = F;
        iter = iter + 1;
        header_cMMN{pos_head} = [participant_group{pg} '_' condition_short_labels{c} '_' Peaks_to_extract{pe}]; 
        pos_head = pos_head + 1;  
    end
    end
end
end

if ~exist([root_dir '/Scalp/' gavr_name '/Statistics'],'dir')
 mkdir([root_dir '/Scalp/' gavr_name '/'], 'Statistics');
end
 

SPSS_matrix = [header_cMMN; num2cell(temporal)];
% Add subject list
% SPSS_matrix = [cell(size(SPSS_matrix,1),1) SPSS_matrix];
% SPSS_matrix{1,1} = 'subj';
% for lat = 1:length(lat_header_cMMN)
%     SPSS_matrix{1+lat,1} = lat_header_cMMN{lat};
% end

if size(chan_num,2) > 1 % Cluster
    save([root_dir '/Scalp/' gavr_name '/Statistics/SPSS_Matrix_' modality_data{mode} '_cluster.mat'],'SPSS_matrix');
    xlswrite([root_dir '/Scalp/' gavr_name '/Statistics/SPSS_Matrix_' modality_data{mode} '_cluster.xlsx'],SPSS_matrix);
else % Single channel
    save([root_dir '/Scalp/' gavr_name '/Statistics/SPSS_Matrix_' modality_data{mode} '_' chan_num{:} '.mat'],'SPSS_matrix');
    xlswrite([root_dir '/Scalp/' gavr_name '/Statistics/SPSS_Matrix_' modality_data{mode} '_' chan_num{:} '.xlsx'],SPSS_matrix);
end
end

% Reset modality data to its original self
modality_data = {'EEG'};

clearvars('-except', initialVars{:});
disp 'DONE WITH EXTRACTING SCALP VALUES OUT OF BRAINSTORM (cMMN)!!!'
disp(datetime)
toc

%% Define Plot-specific variables

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% SPECIFIC PARAMETERS OF THIS SECTION %%%%%%%%%%%%%%%%%%%%%%
set(0,'defaultfigurecolor',[1 1 1]); % I want white backgrounds
abreviated_cond_titles = {'STD','PDEV','DDEV','PitchMMN','DurMMN'};
hemisphere = {'L','R'};
% ROIs = {'MBelt','A1','LBelt','PBelt','STSdp','STSvp'}; % The ones you want to plot...
% ROIs = {'MBelt','A1','LBelt','PBelt','IFG'};
ROIs = {'A1'}; % The ones you want to plot... 'A1','AUDCORTEX' up to six
% 'MBelt','A1','LBelt','PBelt','STSdp','STSvp'
%... based on scouts extracted (Has to be the same than 'Extract source scouts' section)
% 'A1','AUDCORTEX','pOFC','IFG'
scout_positions = {'52_L','A1_L','A4_L','A5_L','LBelt_L', ...
    'MBelt_L','OP4_L','PBelt_L','RI_L','STSdp_L','STSvp_L', ...
    '52_R','A1_R','A4_R','A5_R','LBelt_R','MBelt_R', ...
    'OP4_R','PBelt_R','RI_R','STSdp_R','STSvp_R',...
    '45_L','IFSa_L','IFSp_L','45_R','IFSa_R','IFSp_R' ...
    'OFC_L','pOFC_L','OFC_R','pOFC_R','IFG_L','IFG_R','OFC_L_merged','OFC_R_merged'...
    'AUDCORTEX_L','AUDCORTEX_R','AUDCORTEX3ROI_L','AUDCORTEX3ROI_R'};

normalization_option_tag = 'NORM'; % Or '3D'
scout_function_string = 'mean'; % 'mean','PCA','Max','Std'
dev_GAVR = 2; % 1 = Standard Deviation 2 = Standard Error
s_r = 1000;
overlay_average = 'YES'; % 'NO' or 'YES' for overlay of GAVR over butterfly plot
large_figures = 'NO'; % 'NO' or 'YES' for full screen or regular figure dimensions

% Colors
color_group = {[255 0 0],[0 0 0]}; % FE C
color_condition = {[0 0 0],[255 0 0],[0 128 255]};
color_hem = {[0 0 255],[255 0 0]}; % L R
color_roi = {[255 0 255],[255 0 0],[150 0 0],...
    [255 128 0],[0 153 76],[0 0 153]}; % Depends on how many you put above
transparency = 0.1;
shaded_areas = 0; % 0 = NO shaded areas; 1 = YES
font_size_plots = 24; 
line_width_plots = 2; 
font_size_legend = 24;

% Define shaded areas based on built time-windows
for exptask = 1:length(Exp_tasks)
    if strcmp(Exp_tasks{exptask},'standardMMN')
        condition_short_labels = condition_short_labels_standardMMN_peaks;
        Peaks_to_extract = Peaks_to_extract_standardMMN;
    elseif strcmp(Exp_tasks{exptask},'extratoneMMN')
        condition_short_labels = condition_short_labels_extratoneMMN;
        Peaks_to_extract = Peaks_to_extract_extratoneMMN;
    elseif strcmp(Exp_tasks{exptask},'dualruleMMN')
        condition_short_labels = condition_short_labels_dualruleMMN;
        Peaks_to_extract = Peaks_to_extract_dualruleMMN;
    elseif strcmp(Exp_tasks{exptask},'dualruleMMN_long_SOA')
        condition_short_labels = condition_short_labels_dualruleMMN;   
        Peaks_to_extract = Peaks_to_extract_dualruleMMN;
    else
        error(['No identifiable task']);
    end
    for cmmn = 1:length(condition_short_labels)
        for mode = 1:length(modality_data)
        for pe = 1:length(Peaks_to_extract)
        eval(['shaded_areas_scalp_' condition_short_labels{cmmn} '_' modality_data{mode} '{pe} = [time_window_scalp_' Peaks_to_extract{pe} '_' condition_short_labels{cmmn} '_' ...
            modality_data{mode} '(1) time_window_scalp_' Peaks_to_extract{pe} '_' condition_short_labels{cmmn} '_' modality_data{mode} '(2) time_window_scalp_' Peaks_to_extract{pe} '_' condition_short_labels{cmmn} ...
            '_' modality_data{mode} '(2) time_window_scalp_' Peaks_to_extract{pe} '_' condition_short_labels{cmmn} '_' modality_data{mode} '(1)];']);
        end
        eval(['shaded_areas_source_left_' condition_short_labels{cmmn} '{pe} = [time_window_source_left_' Peaks_to_extract{pe} '_' condition_short_labels{cmmn} ...
            '(1) time_window_source_left_' Peaks_to_extract{pe} '_' condition_short_labels{cmmn} '(2) time_window_source_left_' Peaks_to_extract{pe} '_' condition_short_labels{cmmn} ...
            '(2) time_window_source_left_' Peaks_to_extract{pe} '_' condition_short_labels{cmmn}  '(1)];']);
        eval(['shaded_areas_source_right_' condition_short_labels{cmmn} '{pe} = [time_window_source_right_' Peaks_to_extract{pe} '_' condition_short_labels{cmmn} ...
            '(1) time_window_source_right_' Peaks_to_extract{pe} '_' condition_short_labels{cmmn} '(2) time_window_source_right_' Peaks_to_extract{pe} '_' condition_short_labels{cmmn} ...
            '(2) time_window_source_right_' Peaks_to_extract{pe} '_' condition_short_labels{cmmn}  '(1)];']);
        end
    end
end

save_figures_option = 0; % 0 = NO 1 = YES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Plot GAVR scalp Dual Rule cMMN

% Change modality_data list depending on the choice of sources
modality_data = {'EEG'};

% Select which channels/clusters are going to be used from matrix
load([root_dir '/Scripts/Standard_EEG_channel_order.mat']);
% EEG channels
if size(choice_channel_EEG,2) > 1 % cluster was selected
    channel_pos_EEG = [];
    for i = 1:size(choice_channel_EEG,2) % Find all the indices
        channel_pos_EEG(i) = find(strcmp(std_chan_order, choice_channel_EEG)); 
    end
else % single channel
   channel_pos_EEG = find(strcmp(std_chan_order, choice_channel_EEG)); 
end
load([root_dir '/Scripts/Standard_MEG_channel_order.mat']);
% MEG sensors
if size(choice_channel_MEG,2) > 1 % cluster was selected
    channel_pos_MEG = [];
    for i = 1:size(choice_channel_MEG,2) % Find all the indices
        channel_pos_MEG(i) = find(strcmp(std_chan_order, choice_channel_MEG)); 
    end
else % single sensor
   channel_pos_MEG = find(strcmp(std_chan_order, choice_channel_MEG)); 
end

condition_mismatch_names = {'ABA','AAA','AAA-ABA'}; 
condition_short_labels = {'ABA','AAA','AAASABA'};
abreviated_cond_titles = condition_short_labels;
for mode = 1:length(modality_data)
    % New plot 
    if strcmp(large_figures,'YES')
        figure('units','normalized','outerposition',[0 0 1 1]);
    elseif strcmp(large_figures,'NO')
        figure;
    end
    h(1) = subplot(1,1,1);
    % Define channel number
    eval(['chan_num = channel_pos_' modality_data{mode} ';'])
    current_max = []; current_min = []; pos_ylim = 1; % To adjust ylim later
    
for c = 1:length(condition_mismatch_names)
    if strcmp(condition_mismatch_names{c},'ABA')
    plot_baseline = -50;
    plot_post = 950;
    time_samples=linspace(plot_baseline,plot_post,((((plot_baseline*(-1)) + plot_post)/1000)*s_r) +1);
    elseif strcmp(condition_mismatch_names{c},'AAA')    
    plot_baseline = -50;
    plot_post = 950;
    time_samples=linspace(plot_baseline,plot_post,((((plot_baseline*(-1)) + plot_post)/1000)*s_r) +1);
    elseif strcmp(condition_mismatch_names{c},'AAA-ABA')
    plot_baseline = -50;
    plot_post = 950;
    time_samples=linspace(plot_baseline,plot_post,((((plot_baseline*(-1)) + plot_post)/1000)*s_r) +1);
    end
    
    plot_pos = 1;
    
for pg = 1:length(participant_group)

    % Load data matrix
    load([root_dir '/Scalp/' gavr_name '/gavr/' modality_data{mode} '_' condition_mismatch_names{c} '_' participant_group{pg} '.mat']);
    % Ensure number of channels is correct
    if size(F,1) ~= 306 && strcmp(modality_data{mode},'MEG')
        error(['MEG DATA MATRIX HAS MORE THAN 306 CHANNELS in ' condition_mismatch_names{c} ' ' participant_group{pg}]);
    elseif size(F,1) ~= 62 && strcmp(modality_data{mode},'EEG')
        error(['EEG DATA MATRIX HAS MORE THAN 62 CHANNELS in ' condition_mismatch_names{c} ' ' participant_group{pg}]);
    end

    % Load gavr variable
    average = F(chan_num,:);
    % If cluster was selected, average amplitudes before plotting
    if size(average,1) ~= 1
       average = squeeze(mean(average,1));
    else
       average = squeeze(average);
    end
    % Load stdev or stderr
    if dev_GAVR == 1 % standard deviation
        load([root_dir '/Scalp/' gavr_name '/std_dev/' modality_data{mode} '_' condition_mismatch_names{c} '_' participant_group{pg} '.mat']);
        dev = F(chan_num,:);
    elseif dev_GAVR == 2 % standard error
        load([root_dir '/Scalp/' gavr_name '/std_err/' modality_data{mode} '_' condition_mismatch_names{c} '_' participant_group{pg} '.mat']);
        dev = F(chan_num,:);
    end  
    % Set data ready for plot
    curve1 = average + dev;
    curve2 = average - dev;
    time_samples_2 = [time_samples, fliplr(time_samples)];
    inBetween = [curve1, fliplr(curve2)];
    % Also, grab values for ylim later
    current_max(pos_ylim) = max(curve1);
    current_min(pos_ylim) = min(curve2);
    pos_ylim = pos_ylim + 1; % iterations of participants and conditions
    
    % Now plot
    hold (h(plot_pos),'on')
    fill(h(plot_pos),time_samples_2, inBetween, (color_condition{c}/256), 'FaceAlpha', transparency, 'LineStyle', 'none','HandleVisibility','off');
    plot(h(plot_pos),time_samples, average, 'color', (color_condition{c}/256), 'LineWidth', 1.5);
    
    % Add title (unique of every subplot)
    hold (h(plot_pos),'on')
%     current_title = [abreviated_cond_titles{c}];
%     current_title = strrep(current_title,'_',' ');
%     title (h(plot_pos),current_title);
    
end

end 
% Set ylim
% We set ylim here separataley for EEG and MEG
y_lim_min = min(current_min);
y_lim_max = max(current_max);
for i = 1
    hold (h(i),'on')
    ylim(h(i),[y_lim_min,y_lim_max]);
    line(h(i),[0 0], [y_lim_min y_lim_max],'Color','black')
    line(h(i),[330 330], [y_lim_min y_lim_max],'Color','black')
    line(h(i),[660 660], [y_lim_min y_lim_max],'Color','black')
end
    
    % xlim of DualRuleMMN
    for i = 1
        hold (h(i),'on')
        xlim(h(i),[-50,950]);
        line(h(i),[-50 950],[0 0], 'Color','black')
        set(h(i),'FontName','Arial','FontSize',font_size_plots,'LineWidth', line_width_plots);
        plot(h(i),xlim,[0 0], '-k')
    end
    
if strcmp(modality_data{mode},'EEG')
    hold (h(1),'on')
    ylabel(h(1),'Amplitude (µV)'); 
    xlabel(h(1),'Time (ms)');
    leg = legend(h(1),condition_mismatch_names); % If we wanted it in only one
    set(leg,'FontSize',font_size_legend);
    shaded_areas_scalp_DR_MMN = shaded_areas_scalp_AAASABA_EEG;
elseif strcmp(modality_data{mode},'MEG')
    hold (h(1),'on')
    ylabel(h(1),'Magnetic field (fT)'); 
    xlabel(h(1),'Time (ms)');
    legend(h(1),condition_mismatch_names); % If we wanted it in only one
    set(leg,'FontSize',font_size_legend);
    shaded_areas_scalp_DR_MMN = shaded_areas_scalp_AAASABA_MEG;
end
current_title = ['Grand average sensor level ' modality_data{mode}];
current_title = strrep(current_title,'_',' ');
suptitle(current_title)
    
% Patch MMN areas
if shaded_areas == 1

for i = 1 
    hold (h(i),'on')
    for shad = 1:length(shaded_areas_scalp_DR_MMN)
        gray = [0 0 0];
        patch(h(i),shaded_areas_scalp_DR_MMN{shad},[min(ylim)*[1 1] max(ylim)*[1 1]],gray,'FaceAlpha', transparency,'EdgeAlpha',0.1,'HandleVisibility','off')
    end
end
end

end   

% Reset modality_data variable to its original self
modality_data = {'EEG'};
