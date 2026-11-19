% Author: Athena Stein edited by Sanjala Vasudev
% Date: 11.18.2026
% FOOOF Script for NSCI 270 Project
% SSDE resting-state data




%% Step 1: Import & Setup FOOOF Environment & Output Directory
fprintf('\n------------------------------------------------------------\n');
fprintf(' [Step 1] Initializing Toolboxes & Loading Cleaned Dataset\n');
fprintf('------------------------------------------------------------\n');

% Clean workspace
close all; clc; clear;

% Point EEGLAB to local desktop version and add Toolboxes
ROOT_DESKTOP = '/Users/sanjalavasudev/Desktop/TOOLBOXES/';
addpath(fullfile(ROOT_DESKTOP, 'eeglab2026.0.0'));
ROOT_ONEDRIVE = '/Users/sanjalavasudev/Library/CloudStorage/OneDrive-UniversityofNorthCarolinaatChapelHill/Shortcuts/PandA-tACS - toolboxes';
addpath(fullfile(ROOT_ONEDRIVE, 'dsputils'));
addpath(fullfile(ROOT_ONEDRIVE, 'plotutils'));
addpath(fullfile(ROOT_ONEDRIVE, 'fooof_mat'));
addpath(fullfile(ROOT_ONEDRIVE, 'Riddler_Toolbox'));
addpath(fullfile(ROOT_ONEDRIVE, 'TSST'));

% Start EEGLAB
eeglab;

% Define input paths 
base_name = 'SSDE_P101_RSEEG_EO_AS_20190308_113428';
input_folder = fullfile('/Users/sanjalavasudev/Desktop', 'EEG_Processed_Outputs');

% Define output folder for FOOOF results
fooof_output_folder = fullfile('/Users/sanjalavasudev/Desktop', 'FOOOF nsci 270');
if ~exist(fooof_output_folder, 'dir')
    mkdir(fooof_output_folder);
    fprintf('Created new FOOOF outputs directory: %s\n', fooof_output_folder);
end

% Load final preprocessed clean dataset 
set_filename = [base_name '_ica_step7B.set'];
if ~exist('EEG', 'var') || isempty(EEG.data)
    fprintf('Loading fully cleaned dataset from preprocessing folder...\n');
    EEG = pop_loadset('filename', set_filename, 'filepath', input_folder);
    [ALLEEG, EEG, CURRENTSET] = pop_newset([], EEG, 0, 'gui', 'off');
else
    fprintf('Using dataset currently loaded in the EEGLAB workspace.\n');
end
fprintf('Successfully loaded cleaned dataset: %s (Channels: %d, Rate: %dHz)\n', set_filename, EEG.nbchan, EEG.srate);



%% Step 2: Configure ROI Definitions, Spectral Parameters & FOOOF Config
fprintf('\n------------------------------------------------------------\n');
fprintf(' [Step 2] Configuring Regional ROIs, Pipeline Modes & FOOOF Settings\n');
fprintf('------------------------------------------------------------\n');

% Choose pipeline mode & analysis config
PIPELINE_MODE = 'WELCH';     % 'WELCH' or 'TUKEY'
DO_FFT_CHECK  = true;        % estimate slope from FFT amplitude
subject = base_name;

% FOOOF settings (fixed)
Flim                    = [2 40];
f_range                 = [2 40];
settings                = struct();
settings.max_n_peaks    = 10;
settings.peak_threshold = 1;
settings.peak_width_limits = [0.5 12];
settings.aperiodic_mode = 'fixed';

% Random gallery / persistence config
SAVE_ALL_FITS       = true;   % persist FOOOF result + residuals + spectra
N_RANDOM_PER_ROI    = 3;      % how many random examples per ROI per state to show
RNG_SEED            = 42;     % reproducible sampling

% Extra PSD persistence (non-FOOOFed)
FULL_SAVE_RANGE     = [1 59]; % save raw PSD in this range for downstream band-power

% Define Rest States
restStates = {'EO'};          % State cell array for Step 3A loop

% ROI definitions (index-based)
LPFC_broad = [18 19 20 22 23 24 26 27 28 29 33 34 35];
LPFC       = [24 23 27 19 20 28];
RPFC       = [124 123 117 118 3 4];
CZ_roi     = [7 106 31 80 55];
Fz_roi     = [11 16 12 5 18 10];
P3_roi     = [52 51 58 59];
P4_roi     = [92 96 97 91];

% NEW ROI definitions (label-based)
central_left_labels  = {'E35','E40','E41','E36'};       % central left
central_right_labels = {'E109','E103','E104','E110'};  % central right
roi_names = {'LPFC_broad','LPFC','RPFC','CZ','Fz','P3','P4', ...
    'central_left','central_right'};

% PSD Parameters & Setup per Mode
Fs = EEG.srate;
switch upper(PIPELINE_MODE)
    case 'WELCH'
        win      = hamming(2048);
        noverlap = 1024;
        nfft     = 16384;
        stamp    = sprintf('Welch:\nFs=%d\nWin=2048\nNOverlap=1024\nNFFT=16384', Fs);
        compute_psd = @(x) pwelch(detrend(x), win, noverlap, nfft, Fs);
    case 'TUKEY'
        winSize  = max(256, 2 * Fs);
        taper    = 1;                    % Tukey alpha=1 (Hann)
        win      = tukeywin(winSize, taper);
        lag      = round(Fs/8);
        foi      = Flim(1):0.1:Flim(2);
        stamp    = sprintf('Tukey:\nFs=%d\nWin=%d\nLag=%d\nfoi=2:0.1:40', Fs, winSize, lag);
        compute_psd = @(x) pwelch(detrend(x), win, winSize - lag, foi, Fs);
    otherwise
        error('Unknown PIPELINE_MODE: %s', PIPELINE_MODE);
end

% Full-range raw PSD for non-FOOOFed band analyses (1-59 Hz)
compute_psd_full = @(x) pwelch(detrend(x), hamming(2048), 1024, 16384, Fs);
fprintf('Configured pipeline mode (%s), %d regional ROIs, and spectral parameters.\n', PIPELINE_MODE, numel(roi_names));




%% Step 3A: Run Regional ROI Analysis, FOOOF Fitting & Save Outputs
fprintf('\n------------------------------------------------------------\n');
fprintf(' [Step 3A] Running Regional ROI Analysis & Saving ROI Outputs\n');
fprintf('------------------------------------------------------------\n');

% Pre-flight workspace sanity checks
if ~exist('EEG', 'var') || isempty(EEG)
    error('ERROR: EEG structure is missing from the workspace. Run Step 1 first.');
end
if ~exist('restStates', 'var') || isempty(restStates)
    restStates = {'EO'}; 
    fprintf('restStates was not defined; defaulting to {''EO''}.\n');
end
if ~exist('subject', 'var')
    subject = base_name;
end

MASTER_ROI = table();
gallery_entries = table(); 
out_root = fullfile(fooof_output_folder, PIPELINE_MODE);
if exist(out_root, 'dir') ~= 7
    mkdir(out_root);
    fprintf('Created directory: %s\n', out_root);
end

for r = 1:numel(restStates)
    restState = restStates{r};   
    max_chan = EEG.nbchan;
    
    % ROI bounds struct lookup
    roi_defs = struct();
    roi_defs.LPFC_broad   = intersect(LPFC_broad(:), (1:max_chan).');
    roi_defs.LPFC         = intersect(LPFC(:),       (1:max_chan).');
    roi_defs.RPFC         = intersect(RPFC(:),       (1:max_chan).');
    roi_defs.CZ           = intersect(CZ_roi(:),     (1:max_chan).');
    roi_defs.Fz           = intersect(Fz_roi(:),     (1:max_chan).');
    roi_defs.P3           = intersect(P3_roi(:),     (1:max_chan).');
    roi_defs.P4           = intersect(P4_roi(:),     (1:max_chan).');
    
    allLabs = {EEG.chanlocs.labels};
    idx_cL = [];
    for L = 1:numel(central_left_labels)
        hit = find(strcmpi(allLabs, central_left_labels{L}));
        if ~isempty(hit), idx_cL(end+1) = hit; end
    end
    roi_defs.central_left = idx_cL(:);
    
    idx_cR = [];
    for L = 1:numel(central_right_labels)
        hit = find(strcmpi(allLabs, central_right_labels{L}));
        if ~isempty(hit), idx_cR(end+1) = hit; end
    end
    roi_defs.central_right = idx_cR(:);
    
    out_sub_state = fullfile(out_root, restState);
    out_roi  = fullfile(out_sub_state, 'ROI');
    out_ele  = fullfile(out_sub_state, 'Electrodes');
    if exist(out_roi,'dir') ~= 7, mkdir(out_roi); end
    if exist(out_ele,'dir') ~= 7, mkdir(out_ele); end
    
    roi_metrics = struct([]);     
    roi_compare = [];             
    
    for i = 1:numel(roi_names)
        nm   = roi_names{i};
        raw_idxs = roi_defs.(nm); 
        idxs = double(raw_idxs(:))'; 
        
        fprintf('Processing ROI: %s (Channels: %s)\n', nm, mat2str(idxs));
        
        if isempty(idxs)
            fprintf('   Skipped %s (no matching channels found)\n', nm);
            continue; 
        end
        
        % 1. Extraction
        sub_data = EEG.data(idxs, :);
        sig_roi  = mean(double(sub_data), 1);
        
        % 2. PSD Computation with error isolation
        try
            [pxx_roi, fvec] = compute_psd(sig_roi);
        catch psd_err
            fprintf('CRASH in compute_psd for %s: %s (Line %d)\n', nm, psd_err.message, psd_err.stack(1).line);
            keyboard;
        end
        
        fvec    = double(fvec(:));
        pxx_roi = double(pxx_roi(:));
        fmask = (fvec >= f_range(1)) & (fvec <= f_range(2));
        f_in  = fvec(fmask);
        p_in  = pxx_roi(fmask);
        if ~issorted(f_in)
            [f_in, ord] = sort(f_in,'ascend'); p_in = p_in(ord);
        end
        
        % 3. FOOOF fitting via standalone Python call
        freqs  = f_in(:)';
        powers = p_in(:)';
        
        % Save input spectra for Python script execution
        save(fullfile(fooof_output_folder, 'temp_fooof_input.mat'), 'freqs', 'powers', 'f_range');
        
        % Call python script via terminal command
        python_executable = '/opt/anaconda3/bin/python3';
        python_script     = fullfile('/Users/sanjalavasudev/Desktop', 'run_fooof_batch.py');
        
        [status, cmdout] = system(sprintf('"%s" "%s"', python_executable, python_script));
        
        if status ~= 0
            fprintf('CRASH in Python script execution: %s\n', cmdout);
            keyboard;
        end
        
        % Load fitted results back into MATLAB
        res_py = load(fullfile(fooof_output_folder, 'temp_fooof_results.mat'));
        
        off     = res_py.offset;
        slp     = res_py.slope;
        expnt   = res_py.exponent;
        resid   = res_py.residuals(:);
        R2      = res_py.r_squared;
        RMSE    = sqrt(mean(resid.^2));
        MAD     = mean(abs(resid));
        RSTD    = std(resid);
        
        % Reconstruct result structure for plotting
        res = struct();
        res.freqs            = freqs(:);
        res.power_spectrum   = powers(:);
        res.ap_fit           = res_py.ap_fit(:);
        res.fooofed_spectrum = res_py.fooofed_spectrum(:);
        res.r_squared        = R2;
        
        % Optional FFT check
        slope_fft = NaN; r2_fft = NaN; d_fft_welch = NaN;
        if DO_FFT_CHECK
            N = numel(sig_roi);
            X = fft(sig_roi);
            half_N     = floor(N/2)+1;
            freqs_half = (0:half_N-1).' * (Fs / N);
            amp_half   = abs(X(1:half_N)).' / N;
            pow_half   = amp_half(:).^2;
            fmaskH = (freqs_half >= f_range(1)) & (freqs_half <= f_range(2));
            fH_in  = freqs_half(fmaskH);
            pH_in  = pow_half(fmaskH);
            if ~issorted(fH_in)
                [fH_in, ordH] = sort(fH_in,'ascend'); pH_in = pH_in(ordH);
            end
            freqsH = fH_in(:)';
            powersH = pH_in(:)';
            
            % Save temporary input for FFT check
            freqs = freqsH; powers = powersH;
            save(fullfile(fooof_output_folder, 'temp_fooof_input.mat'), 'freqs', 'powers', 'f_range');
            system(sprintf('"%s" "%s"', python_executable, python_script));
            resH_py = load(fullfile(fooof_output_folder, 'temp_fooof_results.mat'));
            
            slope_fft  = resH_py.slope;
            r2_fft     = resH_py.r_squared;
            d_fft_welch= slope_fft - slp;
            roi_compare= [roi_compare; {subject, restState, nm, slp, R2, slope_fft, r2_fft, d_fft_welch}];
        end
        
        % ROI summary row & saving
        m = struct('subject',subject, 'state',restState, 'roi',string(nm), 'channel',string(''), ...
                   'offset',off, 'slope',slp, 'exponent',expnt, 'r_squared',R2, ...
                   'rmse',RMSE, 'resid_mad',MAD, 'resid_std',RSTD, ...
                   'slope_fft',slope_fft, 'r2_fft',r2_fft, ...
                   'delta_slope_fft_minus_psd',d_fft_welch, ...
                   'freqs_min',min(res.freqs), 'freqs_max',max(res.freqs));
        roi_metrics = [roi_metrics; m];
        
        fig1 = figure('Name', sprintf('FOOOF Fit – %s', nm), 'Color','w','Visible','off');
        plot(res.freqs, res.power_spectrum, 'k', 'DisplayName', 'Original PSD'); hold on;
        plot(res.freqs, res.ap_fit, 'b--', 'DisplayName', 'Aperiodic Fit');
        plot(res.freqs, res.fooofed_spectrum, 'r', 'DisplayName', 'Full Model');
        xlim(Flim); xlabel('Frequency (Hz)'); ylabel('Power [log10 \muV^2/Hz]');
        title(sprintf('%s | %s | R^2=%.2f', subject, nm, R2), 'Interpreter','none');
        legend('show','Location','southwest');
        annotation('textbox',[0.68 0.72 0.3 0.2], 'String', stamp, ...
            'EdgeColor','none','FontName','Consolas','BackgroundColor',[1 1 1 0.6]);
        saveas(fig1, fullfile(out_roi, sprintf('%s_%s_%s_fit.png', subject, restState, nm)));
        savefig(fig1, fullfile(out_roi, sprintf('%s_%s_%s_fit.fig', subject, restState, nm)));
        close(fig1);
        
        [pxx_fullRoi, f_fullRoi] = compute_psd_full(sig_roi);
        mask_fullRoi  = f_fullRoi >= FULL_SAVE_RANGE(1) & f_fullRoi <= FULL_SAVE_RANGE(2);
        f_fullRoi     = f_fullRoi(mask_fullRoi);
        psd_fullRoi   = pxx_fullRoi(mask_fullRoi);
        
        if SAVE_ALL_FITS
            spectra_bundle = struct( ...
                'freqs',            res.freqs(:), ...
                'psd',              res.power_spectrum(:), ...
                'aperiodic_fit',    res.ap_fit(:), ...
                'fooofed_spectrum', res.fooofed_spectrum(:), ...
                'residuals',        resid(:), ...
                'freq_range',       f_range, ...
                'pipeline',         PIPELINE_MODE, ...
                'roi',              nm, ...
                'state',            restState, ...
                'subject',          subject, ...
                'settings',         settings, ...
                'stamp',            stamp, ...
                'saved_full_range', FULL_SAVE_RANGE, ...
                'freqs_full',       f_fullRoi(:), ...
                'psd_full',         psd_fullRoi(:) ...
            );
            save(fullfile(out_roi, sprintf('%s_%s_%s_spectra_bundle.mat', subject, restState, nm)), 'spectra_bundle');
            res_struct = res; residuals = resid;
            save(fullfile(out_roi, sprintf('%s_%s_%s_fit_residuals.mat', subject, restState, nm)), 'res_struct', 'residuals');
            newrow = table(string(subject), string(restState), string(nm), ...
                string(fullfile(out_roi, sprintf('%s_%s_%s_fit_residuals.mat', subject, restState, nm))), ...
                'VariableNames', {'subject','state','roi','resfile'});
            gallery_entries = [gallery_entries; newrow];
        end
    end
    
    fprintf('Step 3A Complete: Successfully finished all regional ROIs for state: %s\n', restState);
end




%% ================================================================
%  Step 4: Extracting, Inspecting, and Visualizing 1/f Slopes
% ================================================================
clear; clc;

output_dir = '/Users/sanjalavasudev/Desktop/FOOOF nsci 270/WELCH/EO/ROI';

% List of regional ROIs
roi_names = {'LPFC_broad', 'LPFC', 'RPFC', 'CZ', 'Fz', 'P3', 'P4', 'central_left', 'central_right'};

slopes = [];
valid_labels = {};

fprintf('\n========================================\n');
fprintf('  Aperiodic 1/f Slopes Across Regional ROIs \n');
fprintf('========================================\n');

for i = 1:numel(roi_names)
    nm = roi_names{i};
    
    % Locate the fit_residuals file generated in Step 3A
    resFiles = dir(fullfile(output_dir, sprintf('*_%s_fit_residuals.mat', nm)));
    
    if ~isempty(resFiles)
        filePath = fullfile(resFiles(1).folder, resFiles(1).name);
        data = load(filePath);
        
        if isfield(data, 'res_struct')
            rs = data.res_struct;
            f = rs.freqs(:);
            ap = rs.ap_fit(:);
            
            % Compute the 1/f slope from the aperiodic fit in log10(frequency) space
            % FOOOF equation: ap_fit = offset - slope * log10(freq)
            p = polyfit(log10(f), ap, 1);
            slope_val = -p(1); % Negative slope represents the 1/f exponent/slope
            
            fprintf('ROI: %-15s | 1/f Slope: %.3f\n', nm, slope_val);
            
            % Collect values for the bar chart plot
            slopes(end+1) = slope_val;
            valid_labels{end+1} = nm;
        else
            fprintf('ROI: %-15s | res_struct not found.\n', nm);
        end
    else
        fprintf('ROI: %-15s | File not found.\n', nm);
    end
end
fprintf('========================================\n');

% --- Generate and Save Regional Slope Bar Chart ---
if ~isempty(slopes)
    fig = figure('Color', 'w', 'Units', 'normalized', 'Position', [0.2 0.2 0.6 0.4]);
    bar(slopes, 'FaceColor', [0.2 0.4 0.7], 'EdgeColor', 'none');
    
    set(gca, 'XTick', 1:numel(valid_labels), 'XTickLabel', valid_labels, 'XTickLabelRotation', 30);
    ylabel('Aperiodic 1/f Slope (\chi)');
    title('SSDE P101: Regional Aperiodic Slopes');
    grid on; box off;
    
    % Save figure files directly to your ROI folder
    fig_png_path = fullfile(output_dir, 'SSDE_P101_Regional_Slopes_BarChart.png');
    fig_fig_path = fullfile(output_dir, 'SSDE_P101_Regional_Slopes_BarChart.fig');
    
    saveas(fig, fig_png_path);
    savefig(fig, fig_fig_path);
    
    fprintf('\nSaved bar chart figure successfully to:\n  %s\n', fig_png_path);
end




%% ================================================================
%  Step 5: Extracting & Converting Beta Sub-Bands to Interpretable Units
% ================================================================
clear; clc; close all;

output_dir = '/Users/sanjalavasudev/Desktop/FOOOF nsci 270/WELCH/EO/ROI';
roi_names = {'LPFC_broad', 'LPFC', 'RPFC', 'CZ', 'Fz', 'P3', 'P4', 'central_left', 'central_right'};

extracted_data = struct();

fprintf('\n=========================================================================\n');
fprintf('  Extracting & Converting Beta Sub-Bands (Linear Units: uV^2/Hz)       \n');
fprintf('=========================================================================\n');

for i = 1:numel(roi_names)
    nm = roi_names{i};
    resFiles = dir(fullfile(output_dir, sprintf('*_%s_fit_residuals.mat', nm)));

    if ~isempty(resFiles)
        filePath = fullfile(resFiles(1).folder, resFiles(1).name);
        data = load(filePath);
        rs = data.res_struct;

        % Quality Check: Verify model fit quality (R^2 > 0.90)
        r2 = rs.r_squared;
        if r2 < 0.90
            fprintf('WARNING: ROI %s failed QC (R^2 = %.3f). Skipping.\n', nm, r2);
            continue;
        end

        % 1. Compute 1/f slope from aperiodic fit in log-log space
        f = rs.freqs(:);
        ap = rs.ap_fit(:);
        p = polyfit(log10(f), ap, 1);
        slope_val = -p(1);

        % 2. Isolate periodic spectrum (flattened residual in log10 space)
        fooofed_spec = rs.fooofed_spectrum(:);
        periodic_log = fooofed_spec - ap;

        % 3. Extract mean periodic power across three beta sub-bands (Log scale)
        mask1 = (f >= 13) & (f <= 16); % Low Beta
        mask2 = (f > 16) & (f <= 20);  % Mid Beta
        mask3 = (f > 20) & (f <= 30);  % High Beta

        b1_log = mean(periodic_log(mask1));
        b2_log = mean(periodic_log(mask2));
        b3_log = mean(periodic_log(mask3));

        % Convert from log10 back to linear power (uV^2/Hz) for interpretability
        % Formula: Linear Power = 10^(log_power)
        b1_linear = 10^(b1_log);
        b2_linear = 10^(b2_log);
        b3_linear = 10^(b3_log);

        if isnan(slope_val) || isnan(b1_linear) || isnan(b2_linear) || isnan(b3_linear)
            fprintf('WARNING: ROI %s contains NaN values. Skipping.\n', nm);
            continue;
        end

        fprintf('ROI: %-15s | Slope: %.3f | Low (13-16Hz): %.4f | Mid (16-20Hz): %.4f | High (20-30Hz): %.4f\n', ...
            nm, slope_val, b1_linear, b2_linear, b3_linear);

        % Store structured data
        extracted_data.(nm).r2 = r2;
        extracted_data.(nm).slope = slope_val;
        extracted_data.(nm).low_beta_linear = b1_linear;
        extracted_data.(nm).mid_beta_linear = b2_linear;
        extracted_data.(nm).high_beta_linear = b3_linear;
    else
        fprintf('ROI: %-15s | File not found.\n', nm);
    end
end
fprintf('=========================================================================\n');

% --- Build Clean Table & Export to CSV ---
roi_fields = fieldnames(extracted_data);
PaperTable = table();

for k = 1:numel(roi_fields)
    r = roi_fields{k};
    PaperTable.ROI{k,1} = r;
    PaperTable.Model_R2(k,1) = extracted_data.(r).r2;
    PaperTable.Aperiodic_Slope(k,1) = extracted_data.(r).slope;
    PaperTable.Low_Beta_Linear_uV2Hz(k,1) = extracted_data.(r).low_beta_linear;
    PaperTable.Mid_Beta_Linear_uV2Hz(k,1) = extracted_data.(r).mid_beta_linear;
    PaperTable.High_Beta_Linear_uV2Hz(k,1) = extracted_data.(r).high_beta_linear;
end

% Save to your folder as a CSV file ready for your paper
csv_path = fullfile(output_dir, 'SSDE_P101_Interpretable_Beta_Power.csv');
writetable(PaperTable, csv_path);

fprintf('\nInterpretable data table exported successfully to:\n  %s\n', csv_path);



%% ================================================================
%  Step 6: Correlation Analysis & Plotting (100% Base MATLAB)
% ================================================================
clear; clc; close all;

output_dir = '/Users/sanjalavasudev/Desktop/FOOOF nsci 270/WELCH/EO/ROI';
roi_names = {'LPFC_broad', 'LPFC', 'RPFC', 'CZ', 'Fz', 'P3', 'P4', 'central_left', 'central_right'};

slopes = [];
b1_linear = [];
b2_linear = [];
b3_linear = [];
valid_labels = {};

fprintf('\n===================================================\n');
fprintf('  Extracting Data & Running Base corrcoef Analysis \n');
fprintf('===================================================\n');

for i = 1:numel(roi_names)
    nm = roi_names{i};
    resFiles = dir(fullfile(output_dir, sprintf('*_%s_fit_residuals.mat', nm)));
    
    if ~isempty(resFiles)
        filePath = fullfile(resFiles(1).folder, resFiles(1).name);
        data = load(filePath);
        rs = data.res_struct;
        
        % Quality Check: Verify model fit quality (R^2 > 0.90)
        if rs.r_squared < 0.90
            continue; 
        end
        
        % 1. Compute 1/f slope from aperiodic fit in log-log space
        f = rs.freqs(:);
        ap = rs.ap_fit(:);
        p = polyfit(log10(f), ap, 1);
        slope_val = -p(1);
        
        % 2. Isolate periodic spectrum & convert to linear power (uV^2/Hz)
        fooofed_spec = rs.fooofed_spectrum(:);
        periodic_log = fooofed_spec - ap;
        
        b1_val = 10^(mean(periodic_log((f >= 13) & (f <= 16))));
        b2_val = 10^(mean(periodic_log((f > 16) & (f <= 20))));
        b3_val = 10^(mean(periodic_log((f > 20) & (f <= 30))));
        
        slopes(end+1) = slope_val;
        b1_linear(end+1) = b1_val;
        b2_linear(end+1) = b2_val;
        b3_linear(end+1) = b3_val;
        valid_labels{end+1} = nm;
    end
end

% --- Compute Pearson Correlations using Base MATLAB corrcoef ---
[R1, P1] = corrcoef(slopes(:), b1_linear(:));
[R2, P2] = corrcoef(slopes(:), b2_linear(:));
[R3, P3] = corrcoef(slopes(:), b3_linear(:));

r1 = R1(1,2); p1 = P1(1,2);
r2 = R2(1,2); p2 = P2(1,2);
r3 = R3(1,2); p3 = P3(1,2);

% --- Multiple Comparisons Correction (Bonferroni) ---
num_tests = 3; % Low, Mid, High Beta sub-bands
p_vals = [p1, p2, p3];
p_bonf = min(p_vals * num_tests, 1.0);

fprintf('\n===================================================\n');
fprintf('  Statistical Correlation & Significance Report    \n');
fprintf('===================================================\n');
fprintf('Low Beta  (13-16 Hz) vs Slope: r = %.3f | p = %.4f (Bonferroni p = %.4f)\n', r1, p1, p_bonf(1));
fprintf('Mid Beta  (16-20 Hz) vs Slope: r = %.3f | p = %.4f (Bonferroni p = %.4f)\n', r2, p2, p_bonf(2));
fprintf('High Beta (20-30 Hz) vs Slope: r = %.3f | p = %.4f (Bonferroni p = %.4f)\n', r3, p3, p_bonf(3));
fprintf('---------------------------------------------------\n');

% --- Generate Publication-Grade 3-Panel Figure with Base Fit Lines ---
fig = figure('Color', 'w', 'Units', 'inches', 'Position', [1 1 13 4.5]);
clf(fig);

c_low = [0.15 0.35 0.65];   % Slate Blue
c_mid = [0.10 0.55 0.45];   % Deep Teal
c_high = [0.75 0.35 0.20];  % Warm Rust

x_line = linspace(min(slopes), max(slopes), 100);

% Subplot 1: Low Beta
subplot(1, 3, 1);
scatter(slopes, b1_linear, 80, c_low, 'filled', 'MarkerFaceAlpha', 0.85);
hold on;
p_fit1 = polyfit(slopes, b1_linear, 1);
plot(x_line, polyval(p_fit1, x_line), 'Color', c_low, 'LineWidth', 2);
xlabel('Aperiodic 1/f Slope (\chi)', 'FontSize', 11, 'FontWeight', 'bold');
ylabel('Periodic Power [\muV^2/Hz]', 'FontSize', 11, 'FontWeight', 'bold');
title(sprintf('Low Beta (13–16 Hz)\nr = %.2f, uncorr p = %.3f\nBonf p = %.3f', r1, p1, p_bonf(1)), 'FontSize', 11);
grid on; box off;

% Subplot 2: Mid Beta
subplot(1, 3, 2);
scatter(slopes, b2_linear, 80, c_mid, 'filled', 'MarkerFaceAlpha', 0.85);
hold on;
p_fit2 = polyfit(slopes, b2_linear, 1);
plot(x_line, polyval(p_fit2, x_line), 'Color', c_mid, 'LineWidth', 2);
xlabel('Aperiodic 1/f Slope (\chi)', 'FontSize', 11, 'FontWeight', 'bold');
ylabel('Periodic Power [\muV^2/Hz]', 'FontSize', 11, 'FontWeight', 'bold');
title(sprintf('Mid Beta (16–20 Hz)\nr = %.2f, uncorr p = %.3f\nBonf p = %.3f', r2, p2, p_bonf(2)), 'FontSize', 11);
grid on; box off;

% Subplot 3: High Beta
subplot(1, 3, 3);
scatter(slopes, b3_linear, 80, c_high, 'filled', 'MarkerFaceAlpha', 0.85);
hold on;
p_fit3 = polyfit(slopes, b3_linear, 1);
plot(x_line, polyval(p_fit3, x_line), 'Color', c_high, 'LineWidth', 2);
xlabel('Aperiodic 1/f Slope (\chi)', 'FontSize', 11, 'FontWeight', 'bold');
ylabel('Periodic Power [\muV^2/Hz]', 'FontSize', 11, 'FontWeight', 'bold');
title(sprintf('High Beta (20–30 Hz)\nr = %.2f, uncorr p = %.3f\nBonf p = %.3f', r3, p3, p_bonf(3)), 'FontSize', 11);
grid on; box off;

sgtitle('SSDE P101: Regional Aperiodic Slope vs. Periodic Beta Power (Bonferroni Corrected)', 'FontSize', 13, 'FontWeight', 'bold');

% Save Figure files directly to your ROI folder
saveas(fig, fullfile(output_dir, 'SSDE_P101_Slope_vs_Beta_Correlation_Bonferroni.png'));
savefig(fig, fullfile(output_dir, 'SSDE_P101_Slope_vs_Beta_Correlation_Bonferroni.fig'));

fprintf('\nCorrelation analysis and figure saved successfully to:\n  %s\n\n', output_dir);




%% ================================================================
%  Step 6B: Native Robust Regression, Spearman Ranks & ROI Labeling
% ================================================================
clear; clc; close all;

output_dir = '/Users/sanjalavasudev/Desktop/FOOOF nsci 270/WELCH/EO/ROI';
roi_names = {'LPFC_broad', 'LPFC', 'RPFC', 'CZ', 'Fz', 'P3', 'P4', 'central_left', 'central_right'};

slopes = [];
b1_linear = [];
b2_linear = [];
b3_linear = [];
valid_labels = {};

fprintf('\n===================================================\n');
fprintf('  Extracting Data & Running Native Robust Analysis \n');
fprintf('===================================================\n');

for i = 1:numel(roi_names)
    nm = roi_names{i};
    resFiles = dir(fullfile(output_dir, sprintf('*_%s_fit_residuals.mat', nm)));
    
    if ~isempty(resFiles)
        filePath = fullfile(resFiles(1).folder, resFiles(1).name);
        data = load(filePath);
        rs = data.res_struct;
        
        % Quality Check: Verify model fit quality (R^2 > 0.90)
        if rs.r_squared < 0.90
            continue; 
        end
        
        % 1. Compute 1/f slope
        f = rs.freqs(:);
        ap = rs.ap_fit(:);
        p = polyfit(log10(f), ap, 1);
        slope_val = -p(1);
        
        % 2. Isolate periodic spectrum & convert to linear power (uV^2/Hz)
        fooofed_spec = rs.fooofed_spectrum(:);
        periodic_log = fooofed_spec - ap;
        
        b1_val = 10^(mean(periodic_log((f >= 13) & (f <= 16))));
        b2_val = 10^(mean(periodic_log((f > 16) & (f <= 20))));
        b3_val = 10^(mean(periodic_log((f > 20) & (f <= 30))));
        
        slopes(end+1) = slope_val;
        b1_linear(end+1) = b1_val;
        b2_linear(end+1) = b2_val;
        b3_linear(end+1) = b3_val;
        valid_labels{end+1} = nm;
    end
end

% --- Native Spearman Rank Correlations (Zero Toolboxes) ---
[~, idx_s] = sort(slopes(:)); [~, r_slopes] = sort(idx_s);
[~, idx_b1] = sort(b1_linear(:)); [~, r_b1] = sort(idx_b1);
[~, idx_b2] = sort(b2_linear(:)); [~, r_b2] = sort(idx_b2);
[~, idx_b3] = sort(b3_linear(:)); [~, r_b3] = sort(idx_b3);

[R_sp1, P_sp1] = corrcoef(r_slopes, r_b1);
[R_sp2, P_sp2] = corrcoef(r_slopes, r_b2);
[R_sp3, P_sp3] = corrcoef(r_slopes, r_b3); % Fixed reference here

rho1 = R_sp1(1,2); sp1 = P_sp1(1,2);
rho2 = R_sp2(1,2); sp2 = P_sp2(1,2);
rho3 = R_sp3(1,2); sp3 = P_sp3(1,2);

fprintf('\n===================================================\n');
fprintf('  Spearman Rank Correlation Report (Native Math)   \n');
fprintf('===================================================\n');
fprintf('Low Beta  (13-16 Hz) vs Slope: rho = %.3f | p = %.4f\n', rho1, sp1);
fprintf('Mid Beta  (16-20 Hz) vs Slope: rho = %.3f | p = %.4f\n', rho2, sp2);
fprintf('High Beta (20-30 Hz) vs Slope: rho = %.3f | p = %.4f\n', rho3, sp3);
fprintf('---------------------------------------------------\n');

% --- Inline Robust Fitting (IRLS with Bisquare weights) ---
x_data = slopes(:);

% Low Beta Robust Fit
y1 = b1_linear(:);
p = polyfit(x_data, y1, 1); rob_s1 = p(1); rob_i1 = p(2);
for iter = 1:10
    res = y1 - (rob_s1 * x_data + rob_i1);
    s = max(1.4826 * median(abs(res - median(res))), 1e-6);
    u = res / (4.685 * s); w = (abs(u) < 1) .* (1 - u.^2).^2;
    beta = ([x_data, ones(size(x_data))]' * diag(w) * [x_data, ones(size(x_data))]) \ ([x_data, ones(size(x_data))]' * diag(w) * y1);
    rob_s1 = beta(1); rob_i1 = beta(2);
end

% Mid Beta Robust Fit
y2 = b2_linear(:);
p = polyfit(x_data, y2, 1); rob_s2 = p(1); rob_i2 = p(2);
for iter = 1:10
    res = y2 - (rob_s2 * x_data + rob_i2);
    s = max(1.4826 * median(abs(res - median(res))), 1e-6);
    u = res / (4.685 * s); w = (abs(u) < 1) .* (1 - u.^2).^2;
    beta = ([x_data, ones(size(x_data))]' * diag(w) * [x_data, ones(size(x_data))]) \ ([x_data, ones(size(x_data))]' * diag(w) * y2);
    rob_s2 = beta(1); rob_i2 = beta(2);
end

% High Beta Robust Fit
y3 = b3_linear(:);
p = polyfit(x_data, y3, 1); rob_s3 = p(1); rob_i3 = p(2);
for iter = 1:10
    res = y3 - (rob_s3 * x_data + rob_i3);
    s = max(1.4826 * median(abs(res - median(res))), 1e-6);
    u = res / (4.685 * s); w = (abs(u) < 1) .* (1 - u.^2).^2;
    beta = ([x_data, ones(size(x_data))]' * diag(w) * [x_data, ones(size(x_data))]) \ ([x_data, ones(size(x_data))]' * diag(w) * y3);
    rob_s3 = beta(1); rob_i3 = beta(2);
end

% --- Generate Labeled, Robust Fit Figure ---
fig = figure('Color', 'w', 'Units', 'inches', 'Position', [1 1 14 5]);
clf(fig);

c_low = [0.15 0.35 0.65];   
c_mid = [0.10 0.55 0.45];   
c_high = [0.75 0.35 0.20];  

x_line = linspace(min(slopes), max(slopes), 100);

% Subplot 1: Low Beta
subplot(1, 3, 1);
scatter(slopes, b1_linear, 80, c_low, 'filled', 'MarkerFaceAlpha', 0.85);
hold on;
plot(x_line, rob_s1 * x_line + rob_i1, 'Color', c_low, 'LineWidth', 2, 'LineStyle', '--');
text(slopes + 0.005, b1_linear, valid_labels, 'FontSize', 8, 'Color', 'k');
xlabel('Aperiodic 1/f Slope (\chi)', 'FontSize', 11, 'FontWeight', 'bold');
ylabel('Periodic Power [\muV^2/Hz]', 'FontSize', 11, 'FontWeight', 'bold');
title(sprintf('Low Beta (13–16 Hz)\nSpearman \\rho = %.2f, p = %.3f', rho1, sp1), 'FontSize', 11);
grid on; box off;

% Subplot 2: Mid Beta
subplot(1, 3, 2);
scatter(slopes, b2_linear, 80, c_mid, 'filled', 'MarkerFaceAlpha', 0.85);
hold on;
plot(x_line, rob_s2 * x_line + rob_i2, 'Color', c_mid, 'LineWidth', 2, 'LineStyle', '--');
text(slopes + 0.005, b2_linear, valid_labels, 'FontSize', 8, 'Color', 'k');
xlabel('Aperiodic 1/f Slope (\chi)', 'FontSize', 11, 'FontWeight', 'bold');
ylabel('Periodic Power [\muV^2/Hz]', 'FontSize', 11, 'FontWeight', 'bold');
title(sprintf('Mid Beta (16–20 Hz)\nSpearman \\rho = %.2f, p = %.3f', rho2, sp2), 'FontSize', 11);
grid on; box off;

% Subplot 3: High Beta
subplot(1, 3, 3);
scatter(slopes, b3_linear, 80, c_high, 'filled', 'MarkerFaceAlpha', 0.85);
hold on;
plot(x_line, rob_s3 * x_line + rob_i3, 'Color', c_high, 'LineWidth', 2, 'LineStyle', '--');
text(slopes + 0.005, b3_linear, valid_labels, 'FontSize', 8, 'Color', 'k');
xlabel('Aperiodic 1/f Slope (\chi)', 'FontSize', 11, 'FontWeight', 'bold');
ylabel('Periodic Power [\muV^2/Hz]', 'FontSize', 11, 'FontWeight', 'bold');
title(sprintf('High Beta (20–30 Hz)\nSpearman \\rho = %.2f, p = %.3f', rho3, sp3), 'FontSize', 11);
grid on; box off;

sgtitle('SSDE P101: Native Robust Fit, Spearman Ranks & ROI Labels', 'FontSize', 13, 'FontWeight', 'bold');

% Save Figure
saveas(fig, fullfile(output_dir, 'SSDE_P101_Native_Robust_Labeled.png'));
savefig(fig, fullfile(output_dir, 'SSDE_P101_Native_Robust_Labeled.fig'));

fprintf('\nAnalysis complete! Robust labeled figure saved to:\n  %s\n\n', output_dir);
