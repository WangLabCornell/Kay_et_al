clear

% Set condition flag: HO = 1 for head-on templates
HO = 1;

%% -------------------- Setup Parameters --------------------
% Grey plotting colors
grey_color_05 = [0.5,0.5,0.5];
grey_color_08 = [0.8,0.8,0.8];

% Get the default blue color
defaultBlue = get(groot, 'defaultAxesColorOrder');
blue = defaultBlue(1, :); % First color is the default blue

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

% Known DNA lengths
if HO == 1
    params.n_plus1 = 1596; % bp.  5' position of the RNA
else
    params.n_plus1 = 738; % bp.  5' position of the RNA
end
params.n_arms = 8330;       % Total base pairs in arms

% Data file path
params.trace_filename = "Figure 1 and 2_Sample Data.dat";
params.n_hybrid = 350;
max_force = [];

% 20% of the vleocity clamp speed
dwell_thresh = 0.05;
width = 0.06;

%% -------------------- Generate Theory Curves --------------------

% Force range for theoretical curves (pN)
t_force = 0:0.1:60;

% Normalized extension (ext/Lc) for each DNA type
t_ds_ext = x_MMS(t_force, params.kBT, params.dsDNA.Lp, params.dsDNA.K0);
t_ds_RD_ext = x_MMS(t_force, params.kBT, params.dsRDNA.Lp, params.dsRDNA.K0);
t_ss_ext = x_FJC(t_force, params.kBT, params.ssDNA.Lp, params.ssDNA.K0);

plot_options.linewidth = 2;

%% -------------------- Load and Process Data --------------------

% Load experimental trace data
data = readtable(params.trace_filename);

% Load theoretical model data
if HO == 1
    theory = readtable("HeadOn_Template_Theory");
else
    theory = readtable("CoDirectional_Template_Theory");
end

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
slice = force > 5;
force = force(slice);
ext_ds = data.Extensionnm(slice);
ssext = data.ssExt(slice);
jindex = data.Jindex(slice);
step = data.Step(slice);

%% -------------------- Format File Title --------------------

% Extract and clean up file title for labeling
file_name = params.trace_filename;
file_title = extractAfter(file_name, 'Data_');
file_title = convertStringsToChars(file_title);

%% Plot force vs extension with max force and sliding

% Create and clear figure 1 for plotting
fig1 = figure(1);
clf

% Plot raw unzipping data
plot(ssext,force,'.','MarkerSize',10)
hold on

% Plot smoothed unzipping trace
plot(movmean(ssext,20), movmean(force,20), 'Color', blue, 'LineWidth', 0.5)

% Plot theoretical ssDNA curve
plot(theory.ssDNAExtension_nm_, theory.Force_pN_, "LineWidth", plot_options.linewidth)

%%% Identify First Encounter with Pol II
thrsh = 20;                  % Force threshold (pN)
min_consecutive = 10;        % Minimum consecutive points above threshold
count = 0;

% Step 1: Find first sustained force deviation
for i = 1:length(force)
    if force(i) > thrsh
        count = count + 1;
        if count == min_consecutive
            idx_19 = i - min_consecutive + 1;
            break;
        end
    else
        count = 0;
    end
end

% Step 2: Find last point just below threshold
tolerance = 0.05;
idx_18 = find(abs(force(1:idx_19) - (thrsh - 3.5)) < tolerance, 1, 'last');

% Step 3: Find first point just above threshold
idx_20 = find(abs(force(1:idx_19) - (thrsh - 2.5)) < tolerance, 1, 'last');
% if ~isempty(idx_20)
%     idx_20 = idx_20 + idx_19;
% end

% Compute average extension and force at Pol II encounter
ext_pt = mean(ssext(idx_18:idx_20));
force_pt = mean(force(idx_18:idx_20));

% Plot Pol II encounter point
plot(ext_pt,force_pt,'r.','MarkerSize',25)

% Find closest data point to Pol II encounter
[~, index] = min(sqrt((ssext - ext_pt).^2 + (force - force_pt).^2));
selected_time_1 = time(index);
jindex_1 = jindex(index);

% Determine Pol II location on DNA
trace.n_last_aligned = get_n_ss(ext_pt,force_pt, params)/2; 


% Calculate transcribed distance based on orientation
if HO == 1
    trace.n_transcribed = params.n_plus1 - trace.n_last_aligned; % Head-On
else 
    trace.n_transcribed = trace.n_last_aligned - params.n_plus1; % Co-Directional
end

trace.n_setZero = trace.n_last_aligned;

%%% User Zoom Interaction to Select End of Sliding

msgbox('Use the zoom tool to zoom into region of interest, the data region before the return to baseline/the end of sliding. Once zoomed hit enter key')
zoom on;
% Reset CurrentCharacter to avoid skipping waitfor
set(gcf, 'CurrentCharacter', 'a');% Set to something other than char(13)
% Press Enter to get out of the zoom mode
waitfor(gcf, 'CurrentCharacter', char(13));
zoom reset
zoom off

% User selects populated area directly before return to baseline
[ext_pt1,force_pt1] = ginput(1);

% Find nearby points within threshold
distances1 = sqrt((ssext - ext_pt1).^2 + (force - force_pt1).^2);
threshold_end_pt = 2;
idx_within_threshold = distances1 < threshold_end_pt;

% Compute average of nearby points
ext_avg = mean(ssext(idx_within_threshold));
force_avg = mean(force(idx_within_threshold));

% Plot selected end point
plot(ext_avg, force_avg, 'r.', 'MarkerSize', 25)

% Find closest data point to selected end
[~, index1] = min(sqrt((ssext - ext_avg).^2 + (force - force_avg).^2));
selected_time_2 = time(index1);
jindex_2 = jindex(index1);

% Determine final Pol II location
trace.n_setEnd = get_n_ss(ext_avg, force_avg, params) / 2;

% Estimate hybridized region at end of sliding
if trace.n_setEnd > params.n_plus1
    ext_pt_TSS = params.n_plus1 * 2 * params.ssDNA.a * x_FJC(force_pt, params.kBT, params.ssDNA.Lp, params.ssDNA.K0);
    trace.n_setEndHybrid = get_n_hybrid(ext_pt_TSS - ...
        2 * trace.n_setZero * x_FJC(force_pt, params.kBT, params.ssDNA.Lp, params.ssDNA.K0) * params.ssDNA.a, ...
        force_pt, params) / 2;
else
    trace.n_setEndHybrid = get_n_hybrid(ext_pt - ...
        2 * trace.n_setZero * x_FJC(force_pt, params.kBT, params.ssDNA.Lp, params.ssDNA.K0) * params.ssDNA.a, ...
        force_pt, params) / 2;
end

% Compute sliding distances
n_Slide = trace.n_setEnd - trace.n_setZero;
n_SlideHybrid = trace.n_setEndHybrid;

%%% Plot Annotations and Highlights
% Plot the +1 theory curve
plot(t_ss_ext*params.n_plus1*2*params.ssDNA.a, ...
    t_force, ...
    '--', ...
    "LineWidth", ...
    2, ...
    "Color", ...
    grey_color_05);

% Plot Pol II transcribed position
plot(t_ss_ext*trace.n_last_aligned*params.ssDNA.a*2, ...
    t_force,'--', ...
    "LineWidth",2, ...
    "Color", ...
    default_MATLAB_colors(1));

% Check if max force is empty
TF = isempty(max_force);
if TF == 1
    max_force = max(force);
end

% Highlight max force
max_index = find(force == max_force);
scatter(ssext(max_index), ...
    max_force, ...
    'r', ...
    'o', ...
    'LineWidth', ...
    3)

% Plot Pol II encounter and end positions
plot(t_ss_ext*trace.n_setZero*params.ssDNA.a*2, ...
    t_force, ...
    '--', ...
    "LineWidth", ...
    2, ...
    "Color", ...
    default_MATLAB_colors(2));
plot(t_ss_ext*trace.n_setEnd*params.ssDNA.a*2, ...
    t_force, ...
    '--', ...
    "LineWidth", ...
    2, ...
    "Color", ...
    default_MATLAB_colors(3));

% Add legend
legend(["Unzip Data", ...
    "Unzip Theory", ...
    "+1 Location", ...
    "Pol II Transcribed Position", ...
    "Pol II Ending Position", ...
    "Displaced Force", ...
    "Slide Start", ...
    "Slide End" ...
    ], ...
    "Location", ...
    "northeast")


% Add annotation message
message = sprintf("+1 location: %d\nTranscribed distance: % 0.1f\nDisplacement Force: % 0.2f", ...
    params.n_plus1, ...
    trace.n_transcribed, ...
    max_force);
text(100, 50, message)


% Final plot formatting

xlim([-100, 3000])
ylim([0, 70])
ylabel("Force (pN)")
xlabel("Extension (nm) no arms")
title(file_title, "Interpreter", "none")

% Set all the fonts in the figure to size 16
set(findall(gcf, '-property', 'FontSize'), 'FontSize', 16)

% Optional: Save figure
% saveas(fig1,file_title + "_Sliding.png") % Save sliding plot 

% %% Optional: Additional View Points of Sliding Selection
% 
% % Define a slice of data where force is above 5 pN
% slice = force > 5;
% 
% % Create and clear figure 3
% fig3 = figure(3);
% clf
% 
% % Subplot 1: Extension vs Force (Main Panel)
% subplot(2, 2, [1 3])
% 
% % Plot raw ssDNA extension vs force
% plot(ssext(slice), force(slice), '.', 'MarkerSize', 10)
% hold on
% 
% % Plot smoothed extension vs force
% plot(movmean(ssext(slice), 20), movmean(force(slice), 20), ...
%     'Color', blue, 'LineWidth', 0.5)
% 
% % Plot theoretical ssDNA curve
% plot(theory.ssDNAExtension_nm_, theory.Force_pN_, ...
%     'LineWidth', plot_options.linewidth)
% 
% % Highlight Pol II encounter and end points
% plot(ext_pt, force_pt, 'r.', 'MarkerSize', 25)
% plot(ext_avg, force_avg, 'r.', 'MarkerSize', 25)
% 
% % Plot +1 theory curve
% plot(t_ss_ext * params.n_plus1 * 2 * params.ssDNA.a, t_force, ...
%     '--', 'LineWidth', 2, 'Color', grey_color_05)
% 
% % Plot Pol II transcribed position
% plot(t_ss_ext * trace.n_last_aligned * params.ssDNA.a * 2, t_force, ...
%     '--', 'LineWidth', 2, 'Color', default_MATLAB_colors(1))
% 
% % Highlight max force point
% scatter(ssext(max_index), max_force, 'r', 'o', 'LineWidth', 3)
% 
% % Plot Pol II sliding start and end positions
% plot(t_ss_ext * trace.n_setZero * params.ssDNA.a * 2, t_force, ...
%     '--', 'LineWidth', 2, 'Color', default_MATLAB_colors(2))
% plot(t_ss_ext * trace.n_setEnd * params.ssDNA.a * 2, t_force, ...
%     '--', 'LineWidth', 2, 'Color', default_MATLAB_colors(3))
% 
% % Add legend
% legend({ ...
%     "Unzip Data", ...
%     "Unzip Theory", ...
%     "+1 Location", ...
%     "Pol II Transcribed Position", ...
%     "Displaced Force", ...
%     "Slide Start", ...
%     "Slide End" ...
%     }, 'Location', 'northeast')
% 
% % Add annotation message
% message = sprintf("+1 location: %d\nTranscribed distance: %.1f\nDisplacement Force: %.2f", ...
%     params.n_plus1, trace.n_transcribed, max_force);
% text(100, 50, message)
% 
% % Set axis labels and limits
% xlim([-100, 3000])
% ylim([0, 70])
% xlabel("Extension (nm) no arms")
% ylabel("Force (pN)")
% title(file_title, "Interpreter", "none")
% 
% % Subplot 2: Time vs Force
% subplot(2, 2, 2)
% plot(time(slice), force(slice))
% hold on
% plot(selected_time_1, force_pt, 'r.', 'MarkerSize', 25)
% plot(selected_time_2, force_pt1, 'r.', 'MarkerSize', 25)
% xlabel('Time (s)')
% ylabel('Force (pN)')
% title(file_title, "Interpreter", "none")
% 
% % Subplot 3: Time vs J-index
% subplot(2, 2, 4)
% plot(time(slice), jindex(slice))
% hold on
% plot(selected_time_1, jindex_1, 'r.', 'MarkerSize', 25)
% plot(selected_time_2, jindex_2, 'r.', 'MarkerSize', 25)
% xlabel('Time (s)')
% ylabel('J-index (bp)')
% 
% % Set font size for all elements
% set(findall(gcf, '-property', 'FontSize'), 'FontSize', 16)
% 
% % Optional: Save figure
% % saveas(fig3, file_title + "_Rate.png")

%% Sliding Analysis - plot dwell time histogram to remove stretching until first slide behavior 
% Sliding Rate
time_sliding = time(index:index1);
jindex_sliding = jindex(index:index1);
dt = time_sliding(2) - time_sliding(1);

% Calculate the velocity and position
[smoothed_J, velocity] = calculate_velocity_position(time_sliding, jindex_sliding, width);

% Create and clear figure 4 for plotting
fig4 = figure(4);
clf

% Define histogram bin edges based on the range of smoothed_J values
bins = round(min(smoothed_J)) - 1 : round(max(smoothed_J)) + 1;

% Plot a horizontal histogram of smoothed_J values
h = histogram(smoothed_J,bins,'Orientation', 'horizontal');

% Convert histogram bin counts to time by multiplying with time step (dt
h.BinCounts = h.BinCounts*dt;

% Identify pause regions in the histogram using a threshold
pauses = get_pauses(h.BinCounts,bins,dwell_thresh);

% Close the figure after extracting pause information
close(fig4);

% Initialize a logical index array to mark regions with pausing
index2 = false(size(smoothed_J));

% Loop through each pause and mark the corresponding time intervals
for i = 1:length(pauses)
    min_index = find((smoothed_J > pauses(i).start),1);
    max_index = find((smoothed_J < (pauses(i).end)+1),1,'last');
    index2 = index2 | ((time_sliding > time_sliding(min_index)) & (time_sliding < time_sliding(max_index)));
end

% Find the index for the first detected pause
first_pause = find(index2 == 1, 1, 'first');

% Check if the first pause is preceded by a significant jump in smoothed_J
if smoothed_J(first_pause) > (smoothed_J(1) + 10)
    % If so, ignore the pause and set end index to 1
    first_pause_end = 1;
else
    % Refine the pause region
    first_pause = find(index2, 1, 'first');
    pause_segment = index2(first_pause:end);
    first_pause_end = find(pause_segment == 0, 1, 'first');
    first_pause_end = first_pause + first_pause_end - 1;

    % Check for short gap followed by another pause
    post_pause = index2(first_pause_end:end);
    next_pause = find(post_pause == 1, 1, 'first');
    if ~isempty(next_pause)
        next_pause = next_pause + first_pause_end - 1;
        if (next_pause - first_pause_end) < 50
            extended_pause = index2(next_pause:end);
            pause_end = find(extended_pause == 0, 1, 'first');
            if ~isempty(pause_end)
                first_pause_end = next_pause + pause_end - 1;
            end
         end
    end
end

% Calculate sliding metrics
dJ = smoothed_J(end) - smoothed_J(first_pause_end); % Total J index change
dT = time_sliding(end) - time_sliding(first_pause_end); % Total time
slide_rate_smoothed = dJ / dT; % Average slide rate

% Compute rate of change (derivative) of smoothed_J
rate_of_change = diff(smoothed_J(first_pause_end:end-1)) ./ ...
                diff(time_sliding(first_pause_end:end-1));

% Compute average rate of change
average_rate = mean(rate_of_change);

%% Plot j vs time during sliding

% Create and clear figure 5 for plotting
fig5 = figure(5);
clf

% Set the title of the plot using the file title from the trace structure
title(file_title, "Interpreter","none")

% Plot the raw sliding data (J index over time)
plot(time_sliding,jindex_sliding)
hold on

% Precompute index ranges
all_data = first_pause_end:length(time_sliding);
pause_1 = 1:first_pause_end;

% Define segments and colors
% Highlight the sliding phase (after the first pause ends) in green
% Highlight the initial stretching phase (before the first pause ends) in red
plot_segments = {
    all_data, 'Green';
    pause_1, 'Red'
};

% Loop through and plot
for i = 1:size(plot_segments, 1)
    idx = plot_segments{i, 1};
    color = plot_segments{i, 2};
    plot(time_sliding(idx), smoothed_J(idx), '.', 'MarkerSize', 10, 'Color', color)
end

% Mark the end of the trace with another large yellow dot
plot(time_sliding(first_pause_end), ...
    smoothed_J(first_pause_end), ...
    'y.', ...
    'MarkerSize', ...
    40)

% Mark the start of sliding with a large yellow dot
plot(time_sliding(end), ...
    smoothed_J(end), ...
    'y.', ...
    'MarkerSize', ...
    40)


% Create a message string with sliding metrics and display it on the plot
message = sprintf( ...
    'Sliding Time (s): %0.1f\nSliding Distance no hybrid (bp): %0.1f', ...
    dT, dJ);
text(time_sliding(1), max(smoothed_J)*0.95, message, 'VerticalAlignment', 'top')

% Label the axes
ylabel('J Index (bp)')
xlabel('Time (s)')

% Add legend
legend({'Raw J Index', 'Sliding Phase', 'Stretching Phase', 'Start/End Markers'}, 'Location', 'best')

% Set font size for all text in the figure
set(findall(gca,'-property','FontSize'),'FontSize',16)

% Optional: Save figure
% saveas(fig5, file_title + "_FinalRate.png")

%% Hybrid Size Analysis – Figure 2

plot_shift = 150; % Number of points to shift for visual clarity

% Define threshold for detecting transition to hybrid region
j_thrsh = trace.n_last_aligned;
min_consecutive = 10;
count = 0;

% Step 1: Find the index where extension exceeds threshold for N consecutive points
for i = 1:length(force)
    if ssext(i) > j_thrsh
        count = count + 1;
        if count == min_consecutive
            j_idx = i - min_consecutive + 1;
            break;
        end
     else
        count = 0;
    end
end


% Step 2: Extract force and extension data after the transition
force_after = force(j_idx:end);
extension_after = ssext(j_idx:end);


% Step 3: Compute hybrid extension and subtract from total extension
extension_hybrid = x_MMS(force_after, params.kBT, params.dsRDNA.Lp, params.dsRDNA.K0) ...
                    * params.dsRDNA.a * params.n_hybrid;

% Step 4: Compute extension due to ssDNA only
extension_fromssDNA = extension_after - extension_hybrid;


% Step 5: Calculate j-index after Pol II using ssDNA model
jindex_after = ((extension_fromssDNA) ./ ...
                (x_FJC(force_after, params.kBT, params.ssDNA.Lp, params.ssDNA.K0) * params.ssDNA.a) ...
                + params.n_hybrid) / 2;


%%% Plotting – Figure 2
% Create and clear figure 2 for plotting
fig2 = figure(2);
clf

% Plot j-index vs force before Pol II encounter
plot(jindex(1:j_idx+plot_shift), ...
    force(1:j_idx+plot_shift), ...
    "LineWidth", ...
    plot_options.linewidth)
hold on

% Plot j-index vs force after Pol II encounter
plot(jindex_after(plot_shift:end), ...
    force_after(plot_shift:end), ...
    "LineWidth", ...
    plot_options.linewidth)

% Plot theoretical j-index vs force curve
plot(theory.j_index_, ...
    theory.Force_pN_, ...
    "LineWidth", ...
    plot_options.linewidth)

% Annotate +1 location and transcribed location
xline(params.n_plus1, '--') % +1 transcription start site
xline(trace.n_last_aligned, '--', 'color', 'red') % Pol II encounter site

% Add Legend
legend(["Unzip prior to polII", ...
    "Unzip after polII", ...
    "Theory", "+1 Location", ...
    "Transcribed location" ...
    ], ...
    "Location", ...
    "northeast")

% Add annotation message
message = sprintf("+1 location: %d\nTranscribed distance: %d\nHybrid size: %d", params.n_plus1,trace.n_transcribed,params.n_hybrid);
text(100,50,message)

% Set plot dimensions and labels
xlim([-100, 3200])
ylim([0,60])
ylabel("Force (pN)")
xlabel("Number of basepairs unzipped")
title(file_title, "Interpreter","none")

% Set font size for all elements
set(findall(gcf,'-property','FontSize'),'FontSize',16)

% Optional: Save figure
% saveas(fig2,file_title + "_JIndex.png") % Save jindex plot

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

function pauses = get_pauses(dwell_hist,bins,dwell_thresh)

    ind = dwell_hist > dwell_thresh;

    % Label regions with unique value
    [labeledVector, numRegions] = bwlabel(ind);

    pauses = [];
    for i = 1:numRegions
        bins_selected = bins(labeledVector == i);
        pauses(i).start = min(bins_selected);
        pauses(i).end = max(bins_selected);
    end

end
