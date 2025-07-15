clear

% Set condition flag: Pro = 1 for (+) TFIIS, (+) RNase H, or (+) RPA condition
Pro = 0;

%% -------------------- Setup Parameters --------------------
% Grey plotting colors
grey_color_05 = [0.5,0.5,0.5];
grey_color_08 = [0.8,0.8,0.8];

% Thermal energy constant (in pN·nm)
params.kBT = 4.11;

% dsDNA mechanical parameters
params.dsDNA.Lp = 42;       % Persistence length (nm)
params.dsDNA.K0 = 1200;     % Stretch modulus (pN)
params.dsDNA.a = 0.338;     % Rise per base pair (nm/bp)

% RNA:DNA hybrid parameters (from Xiang)
params.dsRDNA.Lp = 45;      % Persistence length (nm)
params.dsRDNA.K0 = 700;     % Stretch modulus (pN)
params.dsRDNA.a = 0.30;     % Rise per base pair (nm/bp)

% ssDNA mechanical parameters
params.ssDNA.Lp = 0.848;    % Persistence length (nm)
params.ssDNA.K0 = 504.16;   % Stretch modulus (pN)
params.ssDNA.a = 0.546;     % Rise per base pair (nm/bp)

% Set known DNA values
params.n_plus1 = 1596; % bp.  This is the location of the +1 on the trunk.
params.n_arms = 8330; % bp.  The number of bp in the arms (combined)

% Data file
params.trace_filename = "Figure 3 and 4_Sample Data.dat";


%% -------------------- Generate Theory Curves --------------------

% Force range for theoretical curves (pN)
t_force = 0:0.1:60;

% Normalized extension (ext/Lc) for each DNA type
t_ds_ext = x_MMS(t_force, params.kBT, params.dsDNA.Lp, params.dsDNA.K0);
t_ds_RD_ext = x_MMS(t_force, params.kBT, params.dsRDNA.Lp, params.dsRDNA.K0);
t_ss_ext = x_FJC(t_force, params.kBT, params.ssDNA.Lp, params.ssDNA.K0);

%% -------------------- Load and Process Data --------------------

% Load experimental trace data
data = readtable(params.trace_filename);

% Load theoretical model data
theory = readtable("HeadOn_Template_Theory");

% Calculate ssDNA extension from theory
theory.ssDNAExtension_nm_ = theory.extension_nm_ - ...
    x_MMS(theory.Force_pN_, params.kBT, params.dsDNA.Lp, params.dsDNA.K0) * ...
    params.dsDNA.a * params.n_arms;

% Calculate j-index (normalized ssDNA extension)
theory.j_index_ = theory.ssDNAExtension_nm_ ./ ...
    (x_FJC(theory.Force_pN_, params.kBT, params.ssDNA.Lp, params.ssDNA.K0) * ...
    params.ssDNA.a) / 2;

% Normalize time to start at zero
time = data.Time - data.Time(1);

% Extract relevant columns from data
force = data.Force;
ext_ds = data.Extensionnm;
ssext = data.ssExt;
jindex = data.Jindex;
step = data.Step;

%% -------------------- Define Step Events --------------------

% Define step indices and labels based on condition
if Pro == 1
    step_list = [3, 6, 7, 10, 14];
    step_names = [...
        "Initial unzip", ...
        "Backtrack", ...
        "Reduce force to 1 pN", ...
        "Final unzip"];
else
    step_list = [3, 6, 7, 9, 12];
    step_names = [...
        "Initial unzip", ...
        "Backtrack", ...
        "Reduce force to 1 pN", ...
        "Final unzip"];
end

%% -------------------- Format File Title --------------------

% Extract and clean up file title for labeling
file_name = params.trace_filename;
file_title = extractAfter(file_name, 'Data_');
file_title = convertStringsToChars(file_title);


%% Plot Extension Without Arms
% Create and clear figure 1 for plotting
figure(1);
clf

% Subplot 1: Force vs Time (Step Identification
subplot(3,1,1);
title(params.trace_filename, "Interpreter","none")
hold on

% Plot force vs time for each step interval
for k = 1:(length(step_list)-1)
    slice = data.Step >= step_list(k) & data.Step < step_list(k+1);
    plot(time(slice),force(slice));
end

% Axis labels and legend
ylim([-5,50])
ylabel("Force (pN)")
xlabel("Time (s)")
legend(step_names,'Location','eastoutside')

% Subplot 2 & 3: Force vs Extension (Main Analysis
subplot(3,1,[2 3]);
hold on


for k = 1:(length(step_list)-1)
    slice = data.Step >= step_list(k) & data.Step < step_list(k+1);

    % Step 1: Initial unzip — determine Pol II encounter point
    if k == 1
        slice = slice & force > 5;
        [ext_pt1, force_pt] = get_pt(ssext(slice), force(slice), 22);
        params.n_chase = get_n_ss(ext_pt1, force_pt, params) / 2;
    end

    % Step 3: Backtrack — calculate backtrack distance
    if k == 3
        [ext_pt, force_pt] = get_pt(ssext(slice), force(slice), 20);
        params.backtrack = get_n_ss(ext_pt, force_pt, params) / 2 - params.n_chase;
        params.backtracknm = ext_pt - ext_pt1;

        delta_ext = ext_pt - 2 * params.n_chase * ...
            x_FJC(force_pt, params.kBT, params.ssDNA.Lp, params.ssDNA.K0) * params.ssDNA.a;
        params.backtrack_hybrid = get_n_hybrid(delta_ext, force_pt, params) / 2;
    end

    % Step 4: Final unzip — determine post-replication extension
    if k == 4
        slice = slice & force > 0;
        [ext_pt2, force_pt2] = get_pt(ssext(slice), force(slice), 20);
        params.post_replication_nm = ext_pt - ext_pt2;
        params.post_replication_ssDNA = get_n_ss(ext_pt2, force_pt2, params) / 2;
    end

    % Plot force vs extension for this step
    plot(ssext(slice), force(slice));
end

% Axis formatting
xlim([-100, 3000])
ylim([0,50])
xlabel("Extension with arms removed (nm)")
ylabel("Force (pN)")

% Overlay Theory Curves and Annotations 
% Unzipping theory curve
plot(theory.ssDNAExtension_nm_, ...
    theory.Force_pN_, ...
    "LineWidth", ...
    2, ...
    "Color", ...
    grey_color_05);

% Plot the +1 theory curve
plot(t_ss_ext*params.n_plus1*2*params.ssDNA.a, ...
    t_force, ...
    '--', ...
    "LineWidth", ...
    2, ...
    "Color", ...
    grey_color_05);

% Chased Pol II position
plot(t_ss_ext * params.n_chase * 2 * params.ssDNA.a, t_force, ...
    '--', "LineWidth", 2, "Color", default_MATLAB_colors(1));
annotation('textbox', [0.65, 0.4, 0.25, 0.1], ...
    'String', sprintf("Chased position: %0.0f bp\nTranscribed distance: %0.0f nt", ...
    params.n_chase, params.n_plus1 - params.n_chase), ...
    'FitBoxToText', 'on', ...
    'BackgroundColor', 'white', ...
    'EdgeColor', 'black', ...
    'FontSize', 12, ...
    'Color', default_MATLAB_colors(1));


% Backtracked Pol II (no hybrid)
plot(t_ss_ext * (params.n_chase * 2 + params.backtrack * 2) * params.ssDNA.a, ...
    t_force, '--', "LineWidth", 2, "Color", default_MATLAB_colors(3));

% Backtracked Pol II (full hybrid)
plot(t_ss_ext * (params.n_chase * 2 + params.backtrack_hybrid) * params.ssDNA.a + ...
    t_ds_RD_ext * params.backtrack_hybrid * params.dsRDNA.a, ...
    t_force, '--', "LineWidth", 2, ...
    "Color", lighten_plot_color(default_MATLAB_colors(3), 0.5));
annotation('textbox', [0.65, 0.34, 0.25, 0.1], ...
    'String', sprintf("Backtracked Distance: %0.1f nm\nBacktracked (no hybrid): %0.0f nt\nBacktracked (full hybrid): %0.0f nt", ...
    params.backtracknm, params.backtrack, params.backtrack_hybrid), ...
    'FitBoxToText', 'on', ...
    'BackgroundColor', 'white', ...
    'EdgeColor', 'black', ...
    'FontSize', 12, ...
    'Color', default_MATLAB_colors(3));

% Post-replication unzip position
plot(t_ss_ext * (params.post_replication_ssDNA * 2) * params.ssDNA.a, ...
    t_force, '--', "LineWidth", 2, ...
    "Color", lighten_plot_color(default_MATLAB_colors(6), 0.5));
annotation('textbox', [0.65, 0.48, 0.25, 0.1], ...
    'String', sprintf("2nd Unzip position: %0.0f bp\nTranscription from +1: %0.0f bp\nBacktrack recovery: %0.1f nm", ...
    params.post_replication_ssDNA, ...
    params.n_plus1 - params.post_replication_ssDNA, ...
    params.post_replication_nm), ...
    'FitBoxToText', 'on', ...
    'BackgroundColor', 'white', ...
    'EdgeColor', 'black', ...
    'FontSize', 12, ...
    'Color', default_MATLAB_colors(6));

% Set font size for all elements
set(findall(gcf, '-property', 'FontSize'), 'FontSize', 12);

%% Optional: Plot Total Extension (with arms) 
% % Create and clear figure 5 for plotting
% figure(5);
% clf
% hold on
% 
% % Loop through each step and plot extension vs force
% for k = 1:(length(step_list)-1)
%     slice = data.Step >= step_list(k) & data.Step < step_list(k+1);
%     % For steps 1, 3, and 4, only include data with force > 0
%     if (k == 1) || (k == 3) || (k == 4)
%         slice = slice & force > 0;
%     end
% 
%     % Plot the sliced ext vs force data for each step
%     plot(ext_ds(slice),force(slice));
% end
% 
% % Axis formatting
% xlim([0, 5000])
% xlabel("Extension with arms (nm)")
% ylabel("Force (pN)")

%% Plot Step 4: Unzipping After Backtrack 

% Create and clear figure 2 for plotting
figure(2)
clf

% Slice data for Step 4 where force > 0 pN
k = 4;
slice = data.Step >= step_list(k) & data.Step < step_list(k+1);
slice = slice & force > 0;
hold on

% Plot experimental data: extension vs force
plot(ext_ds(slice), ...
    force(slice), ...
    "Color", ...
    default_MATLAB_colors(6));

% Format Axes and Labels 
title('2nd Unzip - Step 4')
xlim([-100, 6000])
ylim([0,50])
xlabel("Extension with arms removed (nm)")
ylabel("Force (pN)")

smoothed_ds_ext = movmean(ext_ds(slice),50);
smoothed_force = movmean(force(slice),50);

% -------------------- Overlay Theory Curves --------------------
% Unzipping theory curve (based on ssDNA model)
plot(theory.extension_nm_, ...
    theory.Force_pN_, ...
    "LineWidth", ...
    2, ...
    "Color", ...
    grey_color_05);

% dsDNA arms theory curve
plot(t_ds_ext*params.n_arms*params.dsDNA.a, ...
    t_force, ...
    '--', ...
    "LineWidth", ...
    2, ...
    "Color", ...
    grey_color_05);

% Transcription start site (+1) theory curve
plot(t_ss_ext*params.n_plus1*2*params.ssDNA.a + t_ds_ext*params.n_arms*params.dsDNA.a, ...
    t_force, ...
    '--', ...
    "LineWidth", ...
    2, ...
    "Color", ...
    grey_color_08);

index = find(smoothed_ds_ext > 2750,1,"first");

ext_interest = smoothed_ds_ext(index);
force_interest = smoothed_force(index);

plot(ext_interest,force_interest,'r.','MarkerSize',25)

% Set font size for all text in the figure
set(findall(gcf,'-property','FontSize'),'FontSize',14)

% Add legend for clarity
legend('Unzip After Backtrack', ...
    'Unzipping Theory', ...
    'dsDNA y-arm WLC',...
    'Transcription Start Site', ...
    'Extenion = 2750')

% Add annotation for transcribed distance
text(50,45,sprintf("Force at Ext = 2750 nm: %0.2f pN", force_interest));

%% Measure the length of Hybrid Produced During the Force Clamp
% Create and clear figure 3 for plotting
figure(3);
clf

% Step 2: Force clamp phase
k = 2;

% Slice data for Step 2
slice = data.Step >= step_list(k) & data.Step < step_list(k+1);

% Create top subplot for extension
ax1 = subplot(3,1,[1 2]);
title('Force Clamp - Step 2')
hold on

% Smooth extension and force data using moving average
smoothed_ext = movmean(ssext(slice),50);
smoothed_force = movmean(force(slice),50);

% Plot raw extension vs time
plot(time(slice), ...
    ssext(slice), ...
    "Color", ...
    lighten_plot_color(default_MATLAB_colors(k),0.7));

% Plot smoothed extension vs time
plot(time(slice), ...
    smoothed_ext, ...
    "Color", ...
    default_MATLAB_colors(k));

% Determine start index (optional: based on force threshold)
index_start = 1;
time2 = time(slice);
time_start = time2(index_start);

% Measure extension range during clamp
ext_start = min(smoothed_ext);  % Minimum extension (start)
ext_end = max(smoothed_ext);    % Maximum extension (end)

% Mean force during clamp
mean_force = mean(smoothed_force(time2>time_start));

% Add annotation with backtrack measurements
ylabel("Extension with arms removed (nm)")
xlabel("Time (s)")

% Position annotation in top-left corner of plot
NW = [min(xlim) max(ylim)]+[diff(xlim) -diff(ylim)]*0.15;
text(NW(1),NW(2),sprintf("Backtracked distance (no hybrid): %0.0f bp" + ...
    "\nBacktraced distance (full hybrid): %0.0f bp" + ...
    "\n Backtrack Extension Max - Min: %0.0f nm", ...
    params.backtrack, params.backtrack_hybrid,params.backtracknm))

% Create bottom subplot for force
ax2 = subplot(3,1,3);

% Plot raw force vs time
plot(time(slice), ...
    force(slice), ...
    "Color", ...
    lighten_plot_color(default_MATLAB_colors(k),0.7));

% Plot smoothed force vs time
plot(time(slice), ...
    smoothed_force, ...
    "Color", ...
    default_MATLAB_colors(k));

% Add labels
ylabel("Force (pN)")
xlabel("Time (s)")

% Link x-axes of both subplots for synchronized zoom/pan
linkaxes([ax1,ax2],'x')

% Set font size for all elements
set(findall(gcf,'-property','FontSize'),'FontSize',14)

%% Analyze Force Dip During Rezip (Step 3) 
% Create and clear figure 3 for plotting
figure(4)
clf
hold on

% Step 3: Rezip
k = 3;

% Slice data for Step 3 where force > 0
slice = data.Step >= step_list(k) & data.Step < step_list(k+1);
slice = slice & force > 0;

% Plot extension vs force for Step 3
plot(ext_ds(slice), ...
    force(slice), ...
    "Color", ...
    default_MATLAB_colors(6));

% Format Axes and Labels 
title('Rezip - Step 3')
xlim([-1000, 6000])
ylim([0,50])
xlabel("Extension (nm)")
ylabel("Force (pN)")

% -------------------- Overlay Theory Curves --------------------
% unzipping theory curve
plot(theory.extension_nm_, ...
    theory.Force_pN_, ...
    "LineWidth", ...
    2, ...
    "Color", ...
    grey_color_08);

% Arms-only extension theory
arm_ext = t_ds_ext*params.n_arms*params.dsDNA.a;
plot(arm_ext,t_force, ...
    '--', ...
    "LineWidth", ...
    2, ...
    "Color", ...
    grey_color_08);

% +1 site theory curve (arms + transcribed region
plot(t_ss_ext*params.n_plus1*2*params.ssDNA.a + arm_ext, ...
    t_force, ...
    '--', ...
    "LineWidth", ...
    2, ...
    "Color", ...
    grey_color_08);

% Backtracked Pol II (no hybrid)
plot(t_ss_ext*(params.n_chase *2 + params.backtrack * 2)*params.ssDNA.a + arm_ext, ...
    t_force, ...
    '--', ...
    "LineWidth", ...
    2, ...
    "Color", ...
    default_MATLAB_colors(3));

% Backtracked Pol II (full hybrid)
plot(t_ss_ext*(params.n_chase *2 + params.backtrack_hybrid)*params.ssDNA.a...
    +t_ds_RD_ext*(params.backtrack_hybrid)*params.dsRDNA.a + arm_ext, ...
    t_force, ...
    '--', ...
    "LineWidth", ...
    2, ...
    "Color", ...
    lighten_plot_color(default_MATLAB_colors(3),0.5));


% -------------------- Interactive Zoom and User Input --------------------
% Prompt user to zoom into region of interest
msgbox('Use the zoom tool to zoom into region of interest. Once zoomed, hit Enter key');
zoom on;

% Reset CurrentCharacter to avoid skipping waitfor
set(gcf, 'CurrentCharacter', 'a'); % Set to something other than char(13)
waitfor(gcf, 'CurrentCharacter', char(13)); % Wait for Enter key

% Reset zoom after user input
zoom reset;
zoom off;

% Let user click on the force dip point
[ext_dip, force_dip] = ginput(1);

ext_ds_rezip = ext_ds(slice);
force_rezip = force(slice);
% Compute Euclidean distance from selected point to all data points
distances = sqrt((ext_ds_rezip - ext_dip).^2 + (force_rezip - force_dip).^2);

% Define a threshold radius in (nm, pN) space
threshold_end_pt = 8;   % Adjust this value based on data scale

% Find indices of points within the threshold
idx_within_threshold = distances < threshold_end_pt;

% Compute average extension and force of nearby points
ext_avg = mean(ext_ds_rezip(idx_within_threshold));
force_avg = mean(force_rezip(idx_within_threshold));

% Optional: plot the selected region and average point
plot(ext_avg, force_avg, 'ko', 'MarkerFaceColor', 'y', 'DisplayName', 'Average Point');

% Set font size
set(findall(gcf,'-property','FontSize'),'FontSize',12)

% Format Axes and Labels 
xlim([-1000, 6000])
ylim([0,50])

% Add annotation in top-left corner
NW = [min(xlim) max(ylim)]+[diff(xlim) -diff(ylim)]*0.05;
text(NW(1),NW(2),sprintf("Force Dip Value: %0.2f pN", force_dip))

%% Optional: Save Figures
% saveas(figure(1),[file_title '_1_FullTrace.png'])
% saveas(figure(2),[file_title '_2_2ndUnzip.png'])
% saveas(figure(3),[file_title '_3_BacktrackClamp.png'])
% saveas(figure(4),[file_title '_4_Reanneal.png'])
% saveas(figure(5),[file_title '_5_FullTraceTotalExtension.png'])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function color = default_MATLAB_colors(index)
    if index < 1
        index = 1;
    end
    index = mod(index-1,7)+1;
    plot_colors = {[0 0.4470 0.7410],[0.8500 0.3250 0.0980],[0.9290 0.6940 0.1250],[0.4940 0.1840 0.5560],[0.4660 0.6740 0.1880],[0.3010 0.7450 0.9330],[0.6350 0.0780 0.1840]};
    color = plot_colors{index};
end

function [ext_pt, force_pt] = get_pt(ext, force, f)
    % Returns the mean extension and force near a target force value `f`
    % Inputs:
    %   ext   - extension data
    %   force - force data
    %   f     - target force (pN)
    % Outputs:
    %   ext_pt   - mean extension near force f
    %   force_pt - mean force near force f

    mask = force > (f - 0.5) & force < (f + 0.5);
    ext_pt = mean(ext(mask));
    force_pt = mean(force(mask));
end

function n = get_n_ss(ext, force, params)
    % Calculates number of nucleotides in ssDNA from extension and force
    n = ext ./ (x_FJC(force, params.kBT, params.ssDNA.Lp, params.ssDNA.K0) * params.ssDNA.a);
end

function n = get_n_ds(ext, force, params)
    % Calculates number of nucleotides in dsDNA from extension and force
    n = ext ./ (x_MMS(force, params.kBT, params.dsDNA.Lp, params.dsDNA.K0) * params.dsDNA.a);
end

function n = get_n_hybrid(ext, force, params)
    % Estimates total number of nucleotides in a strand that is half ssDNA and half dsDNA
    % Example: ---=== would have n = 6
    x_ss = x_FJC(force, params.kBT, params.ssDNA.Lp, params.ssDNA.K0) * params.ssDNA.a;
    x_ds = x_MMS(force, params.kBT, params.dsRDNA.Lp, params.dsRDNA.K0) * params.dsRDNA.a;
    n = 2 * ext ./ (x_ss + x_ds);
end

function f = get_n_hybrid_fraction(ext, force, n, params)
    % Returns the fraction of a hybrid strand that is ssDNA
    % Example: --=== with n = 5 → f = 0.6
    x_ss = x_FJC(force, params.kBT, params.ssDNA.Lp, params.ssDNA.K0) * params.ssDNA.a;
    x_ds = x_MMS(force, params.kBT, params.dsDNA.Lp, params.dsDNA.K0) * params.dsDNA.a;
    f = (ext / n - x_ss) ./ (x_ds - x_ss);
end

function f = get_n_hybrid_number(ext, force, n, params)
    % Returns the number of ssDNA nucleotides in a hybrid strand of total length n
    x_ss = x_FJC(force, params.kBT, params.ssDNA.Lp, params.ssDNA.K0) * params.ssDNA.a;
    x_ds = x_MMS(force, params.kBT, params.dsDNA.Lp, params.dsDNA.K0) * params.dsDNA.a;
    f = (ext - n * x_ss) ./ (x_ds - x_ss);
end

function Smooth = Decimate_Average(data, SmoothFactor)
    % Downsamples data by averaging over chunks of size SmoothFactor
    numChunks = floor(length(data) / SmoothFactor);
    reshapedData = reshape(data(1:numChunks * SmoothFactor), SmoothFactor, []);
    Smooth = mean(reshapedData);

    % Handle remaining data points at the end
    if mod(length(data), SmoothFactor) ~= 0
        remainingPoints = data(end - mod(length(data), SmoothFactor) + 1:end);
        Smooth = [Smooth, mean(remainingPoints)];
    end
end

function x = x_FJC(F, kBT, Lp, K0)
    % Computes normalized extension using the FJC model
    A = 2 * Lp / kBT * F;
    x = (1 ./ tanh(A) - 1 ./ A) .* (F / K0 + 1);
end

function x = x_MMS(F, kBT, Lp, K0)
    % Computes normalized extension using the MMS model (Biophys J. 72, 1335)
    c = K0 * Lp / kBT;
    f = F / K0;

    a0 = 1/4 - (1/4 + (c + 1) * f) .* ((1 + f).^2);
    a1 = (1 + f) .* ((1 + f) + 1/2 + 2 * (c + 1) * f);
    a2 = -(2 * (1 + f) + 1/4 + (c + 1) * f);

    x = zeros(size(f));
    for i = 1:length(f)
        r = roots([1, a2(i), a1(i), a0(i)]);
        x(i) = min(real(r));  % Choose smallest real root
    end
end

function plot_color_light = lighten_plot_color(plot_color,percent)
    %Lighten plot color by mixing with white.
    plot_color_light = plot_color + percent*(1-plot_color);
end

