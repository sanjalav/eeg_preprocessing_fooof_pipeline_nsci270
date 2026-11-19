% Originally written by Tobias Schwippel, Francesca Pupillo, Mengsen Zhang
% Amended by Athena Stein, Akshada Subramaniyan and Sanjala Vasudev
% Date: 11.18.2026
% Preprocessing Script for SSDE resting-state EEG preprocessing pipeline 
% NSCI 270 Final Project



%% Step 1: Import & Verify Raw Dataset 
% Clean workspace
close all; clc; clear;

% Toolboxes
ROOT_DESKTOP = '/Users/sanjalavasudev/Desktop/TOOLBOXES/';
addpath(fullfile(ROOT_DESKTOP, 'eeglab2026.0.0'));
ROOT_ONEDRIVE = '/Users/sanjalavasudev/Library/CloudStorage/OneDrive-UniversityofNorthCarolinaatChapelHill/Shortcuts/PandA-tACS - toolboxes';
addpath(fullfile(ROOT_ONEDRIVE, 'dsputils'));
addpath(fullfile(ROOT_ONEDRIVE, 'plotutils'));

% Start EEGLAB
eeglab;

% Define paths and base name
filePath = '/Users/sanjalavasudev/Desktop/SSDE_P101_RSEEG_EO_AS_20190308_113428.mff';
[~, base_name, ~] = fileparts(filePath);
output_folder = fullfile('/Users/sanjalavasudev/Desktop', 'EEG_Processed_Outputs');
if ~exist(output_folder, 'dir')
    mkdir(output_folder);
end

% Check if EEG is already in the workspace, otherwise import raw .mff file
if ~exist('EEG', 'var') || isempty(EEG.data)
    fprintf('ℹ️ No EEG structure found in workspace. Importing raw .mff file...\n');
    EEG = pop_mffimport({filePath}, {'code'}, 0, 0);
    [ALLEEG, EEG, CURRENTSET] = pop_newset([], EEG, 0, 'gui', 'off');
else
    fprintf('✅ Using dataset currently loaded in the EEGLAB workspace.\n');
end

% Save the raw baseline dataset checkpoint
[ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, length(ALLEEG), 'savenew', fullfile(output_folder, [base_name '_raw.set']), 'gui', 'off');
fprintf('✅ Successfully saved raw baseline dataset: %s\n', [base_name '_raw.set']);



      %% Step 1.5: Compute, Inspect, and Save Raw PSD (Clean Y-Limits & Inset)
        fprintf('\n------------------------------------------------------------\n');
        fprintf(' [Step 1.5] Computing and Plotting Raw PSD\n');
        fprintf('------------------------------------------------------------\n');
        
        % Define paths and base name
        filePath = '/Users/sanjalavasudev/Desktop/SSDE_P101_RSEEG_EO_AS_20190308_113428.mff';
        [~, base_name, ~] = fileparts(filePath);
        output_folder = fullfile('/Users/sanjalavasudev/Desktop', 'EEG_Processed_Outputs');
        
        if ~exist('EEG', 'var') || isempty(EEG.data)
            raw_set_path = fullfile(output_folder, [base_name '_raw.set']);
            if exist(raw_set_path, 'file')
                EEG = pop_loadset('filename', [base_name '_raw.set'], 'filepath', output_folder);
                fprintf('✅ Loaded raw dataset checkpoint from disk.\n');
            else
                error('❌ No raw dataset found. Please save your raw checkpoint first!');
            end
        end
        
        % Parameters for Welch's PSD
        Fs = EEG.srate;
        winSize = 2 * Fs; 
        win = hanning(winSize); 
        lag = round(Fs / 8); 
        overlap = winSize - lag;
        nfft = 2^nextpow2(winSize);
        
        % Compute Welch's PSD
        [pxx, foi] = pwelch(detrend(EEG.data'), win, overlap, nfft, Fs);
        pxx_db = 10 * log10(pxx + eps);
        
        % Open new figure window
        figure('Name', ['Raw PSD - ' base_name], 'Color', 'w', 'Position', [100, 100, 900, 600]);
        
        % Plot with thin lines matching your inspection scripts
        hLines = plot(foi, pxx_db, 'LineWidth', 0.8);
        for i = 1:length(hLines)
            hLines(i).Tag = EEG.chanlocs(i).labels;
            hLines(i).ButtonDownFcn = @(x,~) disp(x.Tag);
        end
        
        % Set axis bounds (cutting Y at -50 to 60 as requested)
        xlim([1 100]); 
        ylim([-50 60]);
        ylabel('Log Power Spectral Density 10*log_{10}(\mu V^2/Hz)'); 
        xlabel('Frequency [Hz]');
        title(sprintf('Raw PSD - %s', base_name), 'Interpreter', 'none', 'FontSize', 11, 'FontWeight', 'bold');
        grid on;
        
        % Position the channel layout inset in the top right
        axes('Position', [0.75, 0.65, 0.18, 0.22]);
        topoplot([], EEG.chanlocs, 'style', 'blank', 'electrodes', 'on', 'chaninfo', EEG.chaninfo);
        title('Channel Layout', 'FontSize', 9);
        
        % Save figures
        psdFigPath = fullfile(output_folder, [base_name '_raw_PSD.png']);
        saveas(gcf, psdFigPath, 'png');
        savefig(fullfile(output_folder, [base_name '_raw_PSD.fig']));
        fprintf('✅ Saved clean raw PSD plot to: %s\n', psdFigPath);
        drawnow;


%% Step 2. Center-Cut the Data
fprintf('\n------------------------------------------------------------\n');
% Set the resting state duration explicitly first
rsDuration = 240; 
fprintf(' [Step 2] Center-Cutting Data to %d Seconds\n', rsDuration);
fprintf('------------------------------------------------------------\n');

% Convert EEGLAB times from milliseconds to seconds
sampleRate = EEG.srate;
EEG_time_sec = EEG.times / 1000; 
numTimePts = length(EEG_time_sec);

% Get the exact center time point of the recording in seconds
centerTime = EEG_time_sec(round(numTimePts / 2));

% Cut around the center - rsDuration/2 in both directions
cutTiming = [(centerTime - (rsDuration / 2)) ...
    (centerTime + (rsDuration / 2) - (1 / sampleRate))];

% Safety check
cutTiming(1) = max(cutTiming(1), EEG_time_sec(1));
cutTiming(2) = min(cutTiming(2), EEG_time_sec(end));

% Select the specified time window
EEG = pop_select(EEG, 'time', cutTiming);

% Update EEGLAB dataset history and redraw
[ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, length(ALLEEG), 'savenew', fullfile(output_folder, [base_name '_centercut.set']), 'gui', 'off');
eeglab redraw;
fprintf('✅ Successfully center-cut dataset to %.2f seconds.\n', rsDuration);



        %% Step 2.5: Compute and Save Center-Cut PSD (Custom Style)
        fprintf('\n------------------------------------------------------------\n');
        fprintf(' [Step 2.5] Computing and Plotting Center-Cut PSD\n');
        fprintf('------------------------------------------------------------\n');
        
        Fs = EEG.srate;
        winSize = 2 * Fs; 
        win = hanning(winSize); 
        overlap = winSize - round(Fs / 8);
        nfft = 2^nextpow2(winSize);
        
        [pxx, foi] = pwelch(detrend(EEG.data'), win, overlap, nfft, Fs);
        pxx_db = 10 * log10(pxx + eps);
        
        figure('Name', ['Center-Cut PSD - ' base_name], 'Color', 'w', 'Position', [100, 100, 900, 600]);
        hLines = plot(foi, pxx_db, 'LineWidth', 0.8);
        for i = 1:length(hLines)
            hLines(i).Tag = EEG.chanlocs(i).labels;
            hLines(i).ButtonDownFcn = @(x,~) disp(x.Tag);
        end
        
        xlim([1 100]); 
        ylim([-50 60]);
        ylabel('Log Power Spectral Density 10*log_{10}(\mu V^2/Hz)'); 
        xlabel('Frequency [Hz]');
        title(sprintf('Center-Cut PSD (240s) - %s', base_name), 'Interpreter', 'none', 'FontSize', 11, 'FontWeight', 'bold');
        grid on;
        
        axes('Position', [0.75, 0.65, 0.18, 0.22]);
        topoplot([], EEG.chanlocs, 'style', 'blank', 'electrodes', 'on', 'chaninfo', EEG.chaninfo);
        title('Channel Layout', 'FontSize', 9);
        
        % Save with the matching Step 2 naming convention
        saveas(gcf, fullfile(output_folder, [base_name '_Step2_Centercut_PSD.png']), 'png');
        savefig(fullfile(output_folder, [base_name '_Step2_Centercut_PSD.fig']));
        fprintf('✅ Saved Step 2 centercut PSD plot.\n');





%% Step 3. Remove Auxiliary Channels
fprintf('\n------------------------------------------------------------\n');
fprintf(' [Step 3] Removing Auxiliary Channels\n');
fprintf('------------------------------------------------------------\n');

% Select only the primary scalp channels (keeping 1 up to 128)
EEG = pop_select(EEG, 'nochannel', 129:EEG.nbchan);

% Update EEGLAB dataset history and save/redraw
[ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, length(ALLEEG), 'savenew', fullfile(output_folder, [base_name '_noaux.set']), 'gui', 'off');
eeglab redraw;

fprintf('✅ Removed auxiliary channels. Remaining channel count: %d\n', EEG.nbchan);



        %% Step 3.5: Clean 128-Channel Layout & Post-Auxiliary Removal PSD 
        fprintf('\n------------------------------------------------------------\n');
        fprintf(' [Step 3.5] Plotting 128-Channel Layout and Post-Aux PSD\n');
        fprintf('------------------------------------------------------------\n');

        % 1. Plot and save clean 128-channel dot layout
        figure('Name', ['128-Channel Layout - ' base_name], 'Color', 'w', 'Position', [100, 100, 600, 600]);
        topoplot([], EEG.chanlocs, 'style', 'blank', 'electrodes', 'on', 'chaninfo', EEG.chaninfo);
        title(sprintf('128-Channel Scalp Layout: %s', base_name), 'Interpreter', 'none', 'FontSize', 11, 'FontWeight', 'bold');

        layoutFigPath = fullfile(output_folder, [base_name '_Step3_128_channel_layout.png']);
        saveas(gcf, layoutFigPath, 'png');
        savefig(fullfile(output_folder, [base_name '_Step3_128_channel_layout.fig']));
        fprintf('✅ Saved clean 128-channel layout plot.\n');
        drawnow;

        % 2. Compute, plot, and save post-auxiliary removal PSD 
        Fs = EEG.srate;
        winSize = 2 * Fs; 
        win = hanning(winSize); 
        overlap = winSize - round(Fs / 8);
        nfft = 2^nextpow2(winSize);

        [pxx, foi] = pwelch(detrend(EEG.data'), win, overlap, nfft, Fs);
        pxx_db = 10 * log10(pxx + eps);

        figure('Name', ['Post-Aux PSD - ' base_name], 'Color', 'w', 'Position', [100, 100, 900, 600]);
        hLines = plot(foi, pxx_db, 'LineWidth', 0.8);
        for i = 1:length(hLines)
            hLines(i).Tag = EEG.chanlocs(i).labels;
            hLines(i).ButtonDownFcn = @(x,~) disp(x.Tag);
        end

        xlim([1 100]); 
        ylim([-50 60]);
        ylabel('Log Power Spectral Density 10*log_{10}(\mu V^2/Hz)'); 
        xlabel('Frequency [Hz]');
        title(sprintf('Post-Auxiliary Removal PSD (128 Ch): %s', base_name), 'Interpreter', 'none', 'FontSize', 11, 'FontWeight', 'bold');
        grid on;

        axes('Position', [0.75, 0.65, 0.18, 0.22]);
        topoplot([], EEG.chanlocs, 'style', 'blank', 'electrodes', 'on', 'chaninfo', EEG.chaninfo);
        title('Channel Layout', 'FontSize', 9);

        saveas(gcf, fullfile(output_folder, [base_name '_Step3_noaux_PSD.png']), 'png');
        savefig(fullfile(output_folder, [base_name '_Step3_noaux_PSD.fig']));
        fprintf('✅ Saved post-auxiliary removal PSD plot.\n');
        drawnow;


%% Step 4. Filter and Downsample
fprintf('\n------------------------------------------------------------\n');
fprintf(' [Step 4] Filtering and Downsampling Data\n');
fprintf('------------------------------------------------------------\n');

eegFilename = [base_name '_filtered.set'];
eegFile = fullfile(output_folder, eegFilename);

if exist(eegFile, 'file') ~= 2
    highPassFilter = 1;
    lowPassFilter = 80;
    notchFilter = [58 62];
    downSampleRate = 200;
    
    filterOrder = round(EEG.srate * 3.3);
    
    fprintf('Applying high-pass (%d Hz) and low-pass (%d Hz) filters...\n', highPassFilter, lowPassFilter);
    EEG = pop_eegfiltnew(EEG, highPassFilter, lowPassFilter, filterOrder, 0, [], 0);
    
    fprintf('Applying notch filter for power line noise ([%d %d] Hz)...\n', notchFilter(1), notchFilter(2));
    EEG = pop_eegfiltnew(EEG, notchFilter(1), notchFilter(2), filterOrder, 1, [], 0);
    
    fprintf('Resampling data to %d Hz...\n', downSampleRate);
    EEG = pop_resample(EEG, downSampleRate);
    
    [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, length(ALLEEG), 'savenew', eegFile, 'gui', 'off');
    eeglab redraw;
    fprintf('✅ Saved filtered dataset: %s\n', eegFilename);
else
    fprintf('ℹ️ Step 4 filtered file already exists. Loading from disk...\n');
    EEG = pop_loadset(eegFile);
end


       %% Step 4.5: Sanity Check & PSD Plot (Post-Filter)
        fprintf('\n------------------------------------------------------------\n');
        fprintf(' [Step 4.5] Computing and Plotting Post-Filter PSD\n');
        fprintf('------------------------------------------------------------\n');
        
        Flim = [1 60]; taper = 1;
        foi = Flim(1):0.1:Flim(2); Fs = EEG.srate;
        winSize = 2 * Fs; 
        win = hanning(winSize); 
        
        lag = round(Fs / 8); overlap = winSize - lag;
        
        fprintf('Computing Power Spectral Density (PSD)...\n');
        tic;
        [pxx, ~] = pwelch(detrend(EEG.data'), win, overlap, foi, Fs);
        toc;
        
        figure('Name', ['Post-Filter PSD - ' base_name], 'Color', 'w');
        
        % Convert to Decibels (Log Scale) for clean visualization
        pxx_db = 10 * log10(pxx + eps);
        
        freqRange = foi >= 8 & foi <= 12;
        maxPower = max(max(pxx_db(freqRange, :)));
        minPower = min(min(pxx_db(freqRange, :)));
        
        % Plot channels using standard MATLAB colors
        hLines = plot(foi, pxx_db);
        for i = 1:length(hLines)
            hLines(i).Tag = EEG.chanlocs(i).labels;
            hLines(i).ButtonDownFcn = @(x,~) disp(x.Tag);
        end
        
        xlim(Flim); ylim([minPower - 10, maxPower + 20]);
        ylabel('Log Power Spectral Density 10*log_{10}(\mu V^2/Hz)'); 
        xlabel('Frequency [Hz]');
        title(sprintf('After Filter and Downsample - %s', base_name), 'Interpreter', 'none');
        
        % Clean layout inset
        axes('Position', [0.72, 0.62, 0.20, 0.25]);
        topoplot([], EEG.chanlocs, 'style', 'blank', 'electrodes', 'on', 'chaninfo', EEG.chaninfo);
        title('Channel Layout', 'FontSize', 9);
        
        step4_png = fullfile(output_folder, [base_name '_Step4_Filter_PSD.png']);
        step4_fig = fullfile(output_folder, [base_name '_Step4_Filter_PSD.fig']);
        saveas(gcf, step4_png, 'png');
        savefig(step4_fig);
        fprintf('✅ Saved Step 4.5 PSD figure files to: %s\n', output_folder);
        
        drawnow;


        % 1. Set working directory and load the last saved filtered dataset
        filePath = '/Users/sanjalavasudev/Desktop/EEG_Processed_Outputs';
        EEG = pop_loadset('filename', 'SSDE_P101_RSEEG_EO_AS_20190308_113428_filtered.set', 'filepath', filePath);

        % 2. Run bad channel rejection and interpolation code here
        % (e.g., eegplot for visual inspection, flagging bad channels, and running pop_interp)
        % Example placeholder:
        % EEG = eeg_interp(EEG, bad_chans, 'spherical');

        % 3. Save the post-interpolation checkpoint
        [ALLEEG, EEG, CURRENTSET] = pop_saveset( EEG, ...
            'filename', 'SSDE_P101_RSEEG_EO_AS_20190308_113428_interp_round1.set', ...
            'filepath', filePath );


%% Step 5. Two Rounds of Bad Channel Rejection and Interpolation
fprintf('\n------------------------------------------------------------\n');
fprintf(' [Step 5] Two Rounds of Bad Channel Rejection & Interpolation\n');
fprintf('------------------------------------------------------------\n');

prev_eegFile = fullfile(output_folder, [base_name '_filtered.set']);
eegFilename  = [base_name '_filtered_interp.set'];
eegFile      = fullfile(output_folder, eegFilename);

if exist(eegFile, 'file') ~= 2
    EEG = pop_loadset(prev_eegFile);
    
    % --- Round 1: First Bad Channel Rejection ---
    fprintf('Launching eegplot for Round 1 bad channel inspection...\n');
    pop_eegplot(EEG, 1, 1, 1);
    inputStr = input('Round 1 - Which electrodes to remove? (e.g. 3, 4, 5, 20 or "none"): ', 's');
    
    if strcmpi(inputStr, 'none') || isempty(inputStr)
        elecToRej = [];
        fprintf('No channels selected for removal in Round 1.\n');
    else
        elecToRej = str2num(inputStr); 
        fprintf('Interpolating %d bad channels in Round 1...\n', length(elecToRej));
    end
    EEG = pop_interp(EEG, elecToRej, 'spherical');
    
    % --- Sanity Check & PSD Plot: Round 1 ---
    fprintf('Computing PSD after Round 1 interpolation...\n');
    [pxx, ~] = pwelch(detrend(EEG.data'), win, overlap, foi, Fs);
    pxx_db = 10 * log10(pxx + eps);
    
    freqRange = find(foi >= Flim(1) & foi <= Flim(2));
    
    figure('Name', ['Interpolation Round 1 - ' base_name], 'Color', 'w');
    maxPower = max(max(pxx_db(freqRange, :)));
    minPower = min(min(pxx_db(freqRange, :)));
    
    hLines = plot(foi, pxx_db);
    for i = 1:length(hLines)
        hLines(i).Tag = EEG.chanlocs(i).labels;
        hLines(i).ButtonDownFcn = @(x,~) disp(x.Tag);
    end
    
    xlim(Flim); ylim([minPower - 10, maxPower + 20]);
    ylabel('Log Power Spectral Density 10*log_{10}(\mu V^2/Hz)'); 
    xlabel('Frequency [Hz]');
    title(sprintf('After Interpolation Round 1 - %s', base_name), 'Interpreter', 'none');
    
    axes('Position', [0.72, 0.62, 0.20, 0.25]);
    topoplot([], EEG.chanlocs, 'style', 'blank', 'electrodes', 'on', 'chaninfo', EEG.chaninfo);
    title('Channel Layout', 'FontSize', 9);
    
    step5_r1_png = fullfile(output_folder, [base_name '_Step5_Interp_Round1_PSD.png']);
    step5_r1_fig = fullfile(output_folder, [base_name '_Step5_Interp_Round1_PSD.fig']);
    saveas(gcf, step5_r1_png, 'png');
    savefig(step5_r1_fig);
    fprintf('✅ Saved Step 5 Round 1 PSD figure files.\n');
    drawnow;
    
    pop_saveset(EEG, ...
        'filename', [base_name '_interp_round1.set'], ...
        'filepath', output_folder);
    fprintf('✅ Saved Step 5 Round 1 dataset checkpoint (.set/.fdt).\n');


  
    % --- Round 2: Second Bad Channel Rejection ---
    fprintf('Launching eegplot for Round 2 bad channel inspection...\n');
    pop_eegplot(EEG, 1, 1, 1);
    uiwait(gcf);

    % Opens a popup dialog box instead of locking the command window
    answer = inputdlg('Which electrodes to remove? (e.g. 3, 4, 5, 20 or "none"):', ...
        'Round 2 Bad Channels', [1 50], {'none'});

    if isempty(answer)
        inputStr2 = 'none'; % Default if they hit cancel
    else
        inputStr2 = answer{1};
    end

    if strcmpi(inputStr2, 'none') || isempty(inputStr2)
        elecToRej2 = [];
        fprintf('No channels selected for removal in Round 2.\n');
    else
        elecToRej2 = str2num(inputStr2);
        fprintf('Interpolating %d bad channels in Round 2...\n', length(elecToRej2));
    end
    EEG = pop_interp(EEG, elecToRej2, 'spherical');

    
    % --- Sanity Check & PSD Plot: Round 2 ---
    fprintf('Computing PSD after Round 2 interpolation...\n');
    [pxx, ~] = pwelch(detrend(EEG.data'), win, overlap, foi, Fs);
    pxx_db = 10 * log10(pxx + eps);

    freqRange = find(foi >= Flim(1) & foi <= Flim(2));

    figure('Name', ['Interpolation Round 2 - ' base_name], 'Color', 'w');
    maxPower = max(max(pxx_db(freqRange, :)));
    minPower = min(min(pxx_db(freqRange, :)));

    hLines = plot(foi, pxx_db);
    for i = 1:length(hLines)
        hLines(i).Tag = EEG.chanlocs(i).labels;
        hLines(i).ButtonDownFcn = @(x,~) disp(x.Tag);
    end

    xlim(Flim); ylim([minPower - 10, maxPower + 20]);
    ylabel('Log Power Spectral Density 10*log_{10}(\mu V^2/Hz)'); 
    xlabel('Frequency [Hz]');
    title(sprintf('After Interpolation Round 2 - %s', base_name), 'Interpreter', 'none');

    axes('Position', [0.72, 0.62, 0.20, 0.25]);
    topoplot([], EEG.chanlocs, 'style', 'blank', 'electrodes', 'on', 'chaninfo', EEG.chaninfo);
    title('Channel Layout', 'FontSize', 9);

    step5_r2_png = fullfile(output_folder, [base_name '_Step5_Interp_Round2_PSD.png']);
    step5_r2_fig = fullfile(output_folder, [base_name '_Step5_Interp_Round2_PSD.fig']);
    saveas(gcf, step5_r2_png, 'png');
    savefig(step5_r2_fig);
    fprintf('✅ Saved Step 5 Round 2 PSD figure files.\n');
    
    % --- SAVE ROUND 2 CHECKPOINT ---
    pop_saveset(EEG, ...
        'filename', [base_name '_interp_round2.set'], ...
        'filepath', output_folder);
    eeglab redraw;
    fprintf('✅ Completed Round 2 bad channel rejection and saved checkpoint: %s_interp_round2.set\n', base_name);
    drawnow;


    %% Step 6: Global Average Re-referencing
    fprintf('\n------------------------------------------------------------\n');
    fprintf(' [Step 6] Global Average Re-referencing\n');
    fprintf('------------------------------------------------------------\n');

    prev_eegFile = fullfile(output_folder, [base_name '_interp_round2.set']);
    rerefFilename = [base_name '_reref.set'];
    rerefFile     = fullfile(output_folder, rerefFilename);

    if exist(rerefFile, 'file') ~= 2
        EEG = pop_loadset(prev_eegFile);

        % Re-reference all active channels to the common average reference
        EEG = pop_reref(EEG, []);

        % --- Sanity Check & PSD Plot: Re-referencing ---
        fprintf('Computing PSD after Common Average Re-referencing...\n');
        [pxx, ~] = pwelch(detrend(EEG.data'), win, overlap, foi, Fs);
        pxx_db = 10 * log10(pxx + eps);

        freqRange = find(foi >= Flim(1) & foi <= Flim(2));

        figure('Name', ['Re-referenced PSD - ' base_name], 'Color', 'w');
        maxPower = max(max(pxx_db(freqRange, :)));
        minPower = min(min(pxx_db(freqRange, :)));

        hLines = plot(foi, pxx_db);
        for i = 1:length(hLines)
            hLines(i).Tag = EEG.chanlocs(i).labels;
            hLines(i).ButtonDownFcn = @(x,~) disp(x.Tag);
        end

        xlim(Flim); ylim([minPower - 10, maxPower + 20]);
        ylabel('Log Power Spectral Density 10*log_{10}(\mu V^2/Hz)'); 
        xlabel('Frequency [Hz]');
        title(sprintf('After Re-referencing - %s', base_name), 'Interpreter', 'none');

        axes('Position', [0.72, 0.62, 0.20, 0.25]);
        topoplot([], EEG.chanlocs, 'style', 'blank', 'electrodes', 'on', 'chaninfo', EEG.chaninfo);
        title('Channel Layout', 'FontSize', 9);

        step6_png = fullfile(output_folder, [base_name '_Step6_Reref_PSD.png']);
        step6_fig = fullfile(output_folder, [base_name '_Step6_Reref_PSD.fig']);
        saveas(gcf, step6_png, 'png');
        savefig(step6_fig);
        fprintf('✅ Saved Step 6 Re-referencing PSD figure files.\n');

        % --- SAVE Re-referenced CHECKPOINT ---
        pop_saveset(EEG, ...
            'filename', rerefFilename, ...
            'filepath', output_folder);

        [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET, 'gui', 'off');
        eeglab redraw;
        fprintf('✅ Completed re-referencing and saved: %s\n', rerefFilename);
        drawnow;
    else
        fprintf('ℹ️ Step 6 re-reference file already exists. Loading from disk...\n');
        EEG = pop_loadset(rerefFile);
    end



    %% Step 7A: Run ICA & Save Checkpoint
    fprintf('\n------------------------------------------------------------\n');
    fprintf(' [Step 7A] Running Independent Component Analysis (ICA)\n');
    fprintf('------------------------------------------------------------\n');

    % Ensure your setup/path block at the top of your script has been run 
    % so 'output_folder' and 'base_name' are in your workspace!
    prev_eegFile = fullfile(output_folder, [base_name '_reref.set']);
    eegFilename  = [base_name '_ica_step7A.set']; 
    eegFile      = fullfile(output_folder, eegFilename);

    if exist(eegFile, 'file') ~= 2
        EEG = pop_loadset(prev_eegFile);
        [ALLEEG, EEG, ~] = eeg_store(ALLEEG, EEG);

        threshold = 0.04;
        ChannelValues = svd(EEG.data);
        inputRank = sum(ChannelValues(:) > threshold);
        fprintf('Calculated data rank via SVD (threshold = %.2f): %d\n', threshold, inputRank);

        fprintf('Running extended infomax ICA (this may take a moment)...\n');
        tic;
        EEG = pop_runica(EEG, 'icatype', 'runica', 'extended', 1, 'pca', inputRank);
        toc;

        % Save ICA checkpoint using pop_saveset
        pop_saveset(EEG, 'filename', eegFilename, 'filepath', output_folder);
        [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, length(ALLEEG), 'gui', 'off');
        eeglab redraw;
        fprintf('✅ ICA decomposition complete and saved: %s\n', eegFilename);

        % --- Static Component Topographies for Project Report ---
        fprintf('Generating and saving clean static ICA component topographies...\n');
        fig = figure('Name', ['ICA Component Topos - ' base_name], 'Color', 'w', 'Position', [100, 100, 1000, 800]);
        num_comps_to_plot = min(20, size(EEG.icaweights, 1));
        pop_topoplot(EEG, 0, 1:num_comps_to_plot, 'ICA Component Topographies', [4 5], 0, 'electrodes', 'on');

        drawnow; % Forces clean rendering

        step7_topo_png = fullfile(output_folder, [base_name '_Step7A_ICA_Topos.png']);
        step7_topo_fig = fullfile(output_folder, [base_name '_Step7A_ICA_Topos.fig']);

        % Use exportgraphics instead of saveas to avoid capturing window UI lines
        exportgraphics(fig, step7_topo_png, 'Resolution', 300);
        savefig(fig, step7_topo_fig);
        fprintf('✅ Saved clean static ICA component topography figure files.\n');
        close(fig); 

    else
        fprintf('ℹ️ Step 7A ICA file already exists. Loading from disk...\n');
        EEG = pop_loadset(eegFile);
    end



    
    %% Step 7B: Automatic ICA Component Classification & Save Checkpoint
    fprintf('\n------------------------------------------------------------\n');
    fprintf(' [Step 7B] Automatic ICA Component Classification & Removal\n');
    fprintf('------------------------------------------------------------\n');

    prev_eegFile = fullfile(output_folder, [base_name '_ica_step7A.set']);
    eegFilename  = [base_name '_ica_step7B.set']; 
    eegFile      = fullfile(output_folder, eegFilename);

    if exist(eegFile, 'file') ~= 2
        EEG = pop_loadset(prev_eegFile);
        [ALLEEG, EEG, ~] = eeg_store(ALLEEG, EEG);

        % Run ICLabel classification
        fprintf('Running ICLabel to classify independent components...\n');
        EEG = pop_iclabel(EEG, 'default');

        % Apply lab-specific probability thresholds for artifact rejection
        fprintf('Applying lab artifact probability thresholds via pop_icflag...\n');
        EEG = pop_icflag(EEG, [0 0; ...     % Brain
            0.5 1; ...   % Muscle
            0.5 1; ...   % Eye
            0.5 1; ...   % Heart
            0.8 1; ...   % Line Noise
            0.5 1; ...   % Channel Noise
            0 0]);       % Other

        rejected_ICs_automatic = find(EEG.reject.gcompreject > 0);
        count_rejected_ICs_automatic = length(rejected_ICs_automatic);
        fprintf('🔍 Automatically flagged %d artifact components for removal.\n', count_rejected_ICs_automatic);

        % Remove flagged ICs
        EEG = pop_subcomp(EEG, rejected_ICs_automatic, 0);

        % Save final cleaned dataset checkpoint using pop_saveset
        pop_saveset(EEG, 'filename', eegFilename, 'filepath', output_folder);
        [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, length(ALLEEG), 'gui', 'off');
        eeglab redraw;
        fprintf('✅ Step 7B cleaned dataset complete and saved: %s\n', eegFilename);

        fprintf('\n------------------------------------------------------------\n');
        fprintf(' [Step 7B] Generating Validation Plots & Saving Checkpoint\n');
        fprintf('------------------------------------------------------------\n');

        eegFilename  = [base_name '_ica_step7B.set']; 
        eegFile      = fullfile(output_folder, eegFilename);

        % Load your already-saved cleaned dataset if not already in workspace
        if ~exist('EEG', 'var') || isempty(EEG.data)
            EEG = pop_loadset(eegFile);
        end

        %Standalone PSD & Multi-Panel Topography Validation
        fprintf('\n------------------------------------------------------------\n');
        fprintf(' [Step 7B] Generating Standalone PSD and Combined Topographies\n');
        fprintf('------------------------------------------------------------\n');

        eegFilename  = [base_name '_ica_step7B.set']; 
        eegFile      = fullfile(output_folder, eegFilename);

        if ~exist('EEG', 'var') || isempty(EEG.data)
            EEG = pop_loadset(eegFile);
        end

        % Parameters
        Flim = [1 45]; 
        foi = Flim(1):0.1:Flim(2); 
        Fs = EEG.srate;
        winSize = 2 * Fs; 
        win = hanning(winSize); 
        lag = round(Fs / 8); 
        overlap = winSize - lag;

        % Compute PSD
        fprintf('Computing Power Spectral Density (PSD) post-ICA...\n');
        [pxx, ~] = pwelch(detrend(EEG.data'), win, overlap, foi, Fs);
        pxx_db = 10 * log10(pxx + eps);

        % Band power calculations
        alpha_idx = (foi >= 8) & (foi <= 12);
        mean_alpha = mean(pxx_db(alpha_idx, :), 1)';

        beta_idx = (foi >= 13) & (foi <= 30);
        mean_beta = mean(pxx_db(beta_idx, :), 1)';

        broad_idx = (foi >= 1) & (foi <= 45);
        mean_broad = mean(pxx_db(broad_idx, :), 1)';


        % --- FIGURE 1: Standalone Cleaned PSD with Inset ---
        fig_psd = figure('Name', ['Post-ICA Cleaned PSD - ' base_name], 'Color', 'w', 'Position', [100, 100, 900, 600]);
        freqRange = find(foi >= 8 & foi <= 12);
        maxPower = max(max(pxx_db(freqRange, :)));
        minPower = min(min(pxx_db(freqRange, :)));

        hLines = plot(foi, pxx_db);
        for i = 1:length(hLines)
            hLines(i).Tag = EEG.chanlocs(i).labels;
            hLines(i).ButtonDownFcn = @(x,~) disp(x.Tag);
        end
        xlim(Flim); 
        ylim([minPower - 10, maxPower + 20]);
        ylabel('Log Power Spectral Density 10*log_{10}(\mu V^2/Hz)'); 
        xlabel('Frequency [Hz]');
        title(sprintf('Post-ICA Cleaned PSD - %s', base_name), 'Interpreter', 'none', 'FontSize', 11, 'FontWeight', 'bold');
        grid on;

        % Channel layout inset
        axes('Position', [0.72, 0.62, 0.20, 0.25]);
        topoplot([], EEG.chanlocs, 'style', 'blank', 'electrodes', 'on', 'chaninfo', EEG.chaninfo);
        title('Channel Layout', 'FontSize', 9);

        % Save standalone PSD
        psd_png = fullfile(output_folder, [base_name '_Step7B_Cleaned_PSD.png']);
        psd_fig = fullfile(output_folder, [base_name '_Step7B_Cleaned_PSD.fig']);
        saveas(fig_psd, psd_png, 'png');
        savefig(fig_psd, psd_fig);


        % --- FIGURE 2: Combined 3-Panel Topography Figure (Alpha, Beta, Broad-Band) ---
        fig_topos = figure('Name', ['Multi-Band Topographies - ' base_name], 'Color', 'w', 'Position', [150, 150, 1300, 450]);

        % 1. Alpha Topography (8-12 Hz)
        subplot(1, 3, 1);
        topoplot(mean_alpha, EEG.chanlocs, 'maplimits', 'maxmin', 'electrodes', 'on', 'colormap', jet, 'headrad', 'rim');
        colorbar;
        title(sprintf('Alpha Band (8-12 Hz)\n%s', base_name), 'Interpreter', 'none', 'FontSize', 10, 'FontWeight', 'bold');

        % 2. Beta Topography (13-30 Hz)
        subplot(1, 3, 2);
        topoplot(mean_beta, EEG.chanlocs, 'maplimits', 'maxmin', 'electrodes', 'on', 'colormap', jet, 'headrad', 'rim');
        colorbar;
        title(sprintf('Beta Band (13-30 Hz)\n%s', base_name), 'Interpreter', 'none', 'FontSize', 10, 'FontWeight', 'bold');

        % 3. Broad-Band Topography (1-45 Hz)
        subplot(1, 3, 3);
        topoplot(mean_broad, EEG.chanlocs, 'maplimits', 'maxmin', 'electrodes', 'on', 'colormap', jet, 'headrad', 'rim');
        colorbar;
        title(sprintf('Broad-Band (1-45 Hz)\n%s', base_name), 'Interpreter', 'none', 'FontSize', 10, 'FontWeight', 'bold');

        % Save combined topographies
        topos_png = fullfile(output_folder, [base_name '_Step7B_Combined_Topos.png']);
        topos_fig = fullfile(output_folder, [base_name '_Step7B_Combined_Topos.fig']);
        saveas(fig_topos, topos_png, 'png');
        savefig(fig_topos, topos_fig);

        fprintf('✅ Saved standalone PSD and combined 3-panel topography files to: %s\n', output_folder);
        drawnow;

        fprintf('\n============================================================\n');
        fprintf(' 🎉 Fully automated SSDE preprocessing pipeline finished successfully!\n');
        fprintf('============================================================\n');




%% ================================================================
%  Master Figure 1 Generator: Dataset Loading & 4-Panel Plot 
% ================================================================
clear; clc; close all;

% Point EEGLAB to your local desktop version and initialize
ROOT_DESKTOP = '/Users/sanjalavasudev/Desktop/TOOLBOXES/';
addpath(fullfile(ROOT_DESKTOP, 'eeglab2026.0.0'));
[ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab; % Initialize EEGLAB

output_folder = '/Users/sanjalavasudev/Desktop/EEG_Processed_Outputs';
base_name = 'SSDE_P101_RSEEG_EO_AS_20190308_113428';

fprintf('Initializing master 4-panel Figure 1 window...\n');
master_fig = figure('Color', 'w', 'Units', 'inches', 'Position', [1 1 18 4.5]);
clf(master_fig);

Flim = [1 60]; % Frequency range to display
Ylim_unified = [-60 60]; % Uniform Y-axis scale for clean comparison

% --- PANEL A: Raw Baseline ---
fprintf('Processing Panel A: Raw Baseline...\n');
sub1 = subplot(1, 4, 1, 'Parent', master_fig);
EEG_raw = pop_loadset('filename', [base_name '_raw.set'], 'filepath', output_folder);
[pxx, foi] = pwelch(detrend(EEG_raw.data'), 2*EEG_raw.srate, EEG_raw.srate, 2^nextpow2(2*EEG_raw.srate), EEG_raw.srate);
pxx_db = 10 * log10(pxx + eps);
plot(sub1, foi, pxx_db, 'LineWidth', 0.6);
title(sub1, 'A. Raw Baseline', 'FontSize', 10, 'FontWeight', 'bold');
xlabel(sub1, 'Frequency [Hz]', 'FontSize', 9);
ylabel(sub1, '10*log_{10}(\mu V^2/Hz)', 'FontSize', 9);
xlim(sub1, Flim); ylim(sub1, Ylim_unified); grid(sub1, 'on'); box(sub1, 'off');

% --- PANEL B: Filter & Notch ---
fprintf('Processing Panel B: Filter & Notch...\n');
sub2 = subplot(1, 4, 2, 'Parent', master_fig);
EEG_filt = pop_loadset('filename', [base_name '_filtered.set'], 'filepath', output_folder);
[pxx, foi] = pwelch(detrend(EEG_filt.data'), 2*EEG_filt.srate, EEG_filt.srate, 2^nextpow2(2*EEG_filt.srate), EEG_filt.srate);
pxx_db = 10 * log10(pxx + eps);
plot(sub2, foi, pxx_db, 'LineWidth', 0.6);
title(sub2, 'B. Filter & Notch', 'FontSize', 10, 'FontWeight', 'bold');
xlabel(sub2, 'Frequency [Hz]', 'FontSize', 9);
ylabel(sub2, '10*log_{10}(\mu V^2/Hz)', 'FontSize', 9);
xlim(sub2, Flim); ylim(sub2, Ylim_unified); grid(sub2, 'on'); box(sub2, 'off');

% --- PANEL C: Interpolation Round 2 ---
fprintf('Processing Panel C: Interpolation...\n');
sub3 = subplot(1, 4, 3, 'Parent', master_fig);
EEG_interp = pop_loadset('filename', [base_name '_interp_round2.set'], 'filepath', output_folder);
[pxx, foi] = pwelch(detrend(EEG_interp.data'), 2*EEG_interp.srate, EEG_interp.srate, 2^nextpow2(2*EEG_interp.srate), EEG_interp.srate);
pxx_db = 10 * log10(pxx + eps);
plot(sub3, foi, pxx_db, 'LineWidth', 0.6);
title(sub3, 'C. Interpolation', 'FontSize', 10, 'FontWeight', 'bold');
xlabel(sub3, 'Frequency [Hz]', 'FontSize', 9);
ylabel(sub3, '10*log_{10}(\mu V^2/Hz)', 'FontSize', 9);
xlim(sub3, Flim); ylim(sub3, Ylim_unified); grid(sub3, 'on'); box(sub3, 'off');

% --- PANEL D: ICA Cleaned ---
fprintf('Processing Panel D: ICA Cleaned...\n');
sub4 = subplot(1, 4, 4, 'Parent', master_fig);
EEG_clean = pop_loadset('filename', [base_name '_ica_step7B.set'], 'filepath', output_folder);
[pxx, foi] = pwelch(detrend(EEG_clean.data'), 2*EEG_clean.srate, EEG_clean.srate, 2^nextpow2(2*EEG_clean.srate), EEG_clean.srate);
pxx_db = 10 * log10(pxx + eps);
plot(sub4, foi, pxx_db, 'LineWidth', 0.6);
title(sub4, 'D. ICA Cleaned', 'FontSize', 10, 'FontWeight', 'bold');
xlabel(sub4, 'Frequency [Hz]', 'FontSize', 9);
ylabel(sub4, '10*log_{10}(\mu V^2/Hz)', 'FontSize', 9);
xlim(sub4, Flim); ylim(sub4, Ylim_unified); grid(sub4, 'on'); box(sub4, 'off');

% Unified super title
sgtitle('Figure 1: Step-by-Step EEG Preprocessing & Signal Stabilization Pipeline', 'FontSize', 12, 'FontWeight', 'bold');

% Save master figure directly to your desktop folder
master_png = fullfile(output_folder, 'Figure1_4Panel_Preprocessing.png');
master_fig_file = fullfile(output_folder, 'Figure1_4Panel_Preprocessing.fig');

exportgraphics(master_fig, master_png, 'Resolution', 300);
savefig(master_fig, master_fig_file);

fprintf('\n✅ Success! Polished master 4-panel Figure 1 generated and saved to:\n  %s\n\n', master_png);