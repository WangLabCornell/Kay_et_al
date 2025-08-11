clear
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
params.hyDNA.Lp = 45;      % Persistence length (nm)
params.hyDNA.K0 = 700;     % Stretch modulus (pN)
params.hyDNA.a = 0.30;     % Rise per base pair (nm/bp)

% ssDNA mechanical parameters
params.ssDNA.Lp = 0.848;    % Persistence length (nm)
params.ssDNA.K0 = 504.16;   % Stretch modulus (pN)
params.ssDNA.a = 0.546;     % Rise per base pair (nm/bp)

% Set known DNA values
params.n_plus1 = 1596; % bp.  This is the location of the +1 on the trunk.
params.n_arms = 8330; % bp.  The number of bp in the arms (combined)

% Data file
params.trace_filename = "Figure 5_Sample Data.dat";

%% -------------------- Generate Theory Curves --------------------

% Force range for theoretical curves (pN)
t_force = 0:0.1:60;

% Normalized extension (ext/Lc) for each DNA type
t_ds_ext = x_MMS(t_force, params.kBT, params.dsDNA.Lp, params.dsDNA.K0);
t_hy_ext = x_MMS(t_force, params.kBT, params.hyDNA.Lp, params.hyDNA.K0);
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
% Define step indices
step_list = [3,6,7,8,9,12];
step_names = [...
    "Initial unzip",...
    "Backtrack",...
    "Reduce force to 2 pN",...
    "Replication at 1 pN",...
    "Final unzip"];
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

% Subplot 2 & 3: Force vs Extension (Main Analysis)
subplot(3,1,[2 3]);
hold on


% Subplot 2 & 3: Force vs Extension (Main Analysis)
subplot(3,1,[2 3]);
hold on

% Loop through each step interval and analyze force-extension behavior
for k = 1:(length(step_list)-1)
    slice = data.Step >= step_list(k) & data.Step < step_list(k+1);

    % Step 1: Identify chased PolII position
    if k == 1
        slice = slice & force > 5;
        [ext_pt, force_pt] = get_pt(ssext(slice), force(slice), 20);
        params.n_chase = get_n_ss(ext_pt, force_pt, params) / 2;
    end

    % Step 3: Identify backtracking behavior
    if k == 3
        [ext_pt, force_pt] = get_pt(ssext(slice), force(slice), 20);
        params.backtrack = get_n_ss(ext_pt, force_pt, params) / 2 - params.n_chase;

        % Calculate extension due to backtracked hybrid
        delta_ext = ext_pt - 2 * params.n_chase * ...
            x_FJC(force_pt, params.kBT, params.ssDNA.Lp, params.ssDNA.K0) * params.ssDNA.a;
        params.backtrack_hybrid = get_n_hybrid(delta_ext, force_pt, params) / 2;
    end

    % Step 5: Analyze post-replication ssDNA and replication fraction
    if k == 5
        slice = slice & force > 0.5 & ssext < 1800;
        [ext_pt, force_pt] = get_pt(ssext(slice), force(slice), 20);
        params.post_replication_ssDNA = get_n_ss(ext_pt, force_pt, params) / 2;

        % Estimate replicated fraction (no further backtracking)
        delta_ext = ext_pt ...
             - params.n_chase * x_FJC(force_pt, params.kBT, params.ssDNA.Lp, params.ssDNA.K0) * params.ssDNA.a ...
            - params.backtrack_hybrid * x_FJC(force_pt, params.kBT, params.ssDNA.Lp, params.ssDNA.K0) * params.ssDNA.a ...
            - params.backtrack_hybrid * x_MMS(force_pt, params.kBT, params.hyDNA.Lp, params.hyDNA.K0) * params.hyDNA.a;

        params.fraction_replicated = get_n_hybrid_fraction(delta_ext, force_pt, params.n_chase, params);
        params.fraction_replicated_number = get_n_hybrid_number(delta_ext, force_pt, params.n_chase, params);

        % Estimate replicated fraction assuming full backtrack
        delta_ext = ext_pt ...
            - params.n_plus1 * x_FJC(force_pt, params.kBT, params.ssDNA.Lp, params.ssDNA.K0) * params.ssDNA.a ...
            - (params.n_plus1 - params.n_chase) * x_MMS(force_pt, params.kBT, params.dsDNA.Lp, params.dsDNA.K0) * params.dsDNA.a;

        params.fraction_replicated_full_backtrack = get_n_hybrid_fraction(delta_ext, force_pt, params.n_chase, params);

        % Store max force and corresponding extension
        max_force = max(force(slice));
        max_index = find(force(slice) == max_force);
        ext2 = ssext(slice);
        max_force_ext = ext2(max_index);
    end

    % Plot force vs extension for current step
    plot(ssext(slice), force(slice));
end

% Format force-extension plot
xlim([-100, 3000])
ylim([0,max_force+10])
xlabel("Extension with arms removed (nm)")
ylabel("Force (pN)")

% -----------Add Theory Curves to Force-Extension Plot------------------
% Plot ssDNA unzipping theory curve
plot(theory.ssDNAExtension_nm_,theory.Force_pN_,"LineWidth",2, "Color",grey_color_05);
hold on

% Plot +1 theory curve (full transcription)
plot(t_ss_ext*params.n_plus1*2*params.ssDNA.a, ...
    t_force, ...
    '--', ...
    "LineWidth", ...
    2, ...
    "Color", ...
    grey_color_08);



% Plot chased PolII position (partial transcription)
plot(t_ss_ext*params.n_chase*params.ssDNA.a*2, ...
    t_force, ...
    '--', ...
    "LineWidth", ...
    2, ...
    "Color", ...
    default_MATLAB_colors(1));

annotation('textbox', [0.6, 0.4, 0.25, 0.1], ...
    'String', sprintf("Chased position: %0.0f bp\nTranscribed distance: %0.0f nt", ...
    params.n_chase, params.n_plus1 - params.n_chase), ...
    'FitBoxToText', 'on', ...
    'BackgroundColor', 'white', ...
    'EdgeColor', 'black', ...
    'FontSize', 12, ...
    'Color', default_MATLAB_colors(1));

% Plot Backtracked PolII Positions Using Two Models
% 1) Model: RNA does NOT hybridize (pure ssDNA backtrack)
plot(t_ss_ext*(params.n_chase *2 + params.backtrack * 2)*params.ssDNA.a, ...
    t_force, ...
    '--', ...
    "LineWidth", ...
    2, ...
    "Color", ...
    default_MATLAB_colors(3));

% 2) Model: RNA fully hybridizes (ssDNA + hybrid DNA)
plot(t_ss_ext*(params.n_chase *2 + params.backtrack_hybrid)*params.ssDNA.a...
    +t_hy_ext*(params.backtrack_hybrid)*params.hyDNA.a, ...
    t_force, ...
    '--', ...
    "LineWidth", ...
    2, ...
    "Color", ...
    lighten_plot_color(default_MATLAB_colors(3),0.5));

% Annotate backtracked distances for both models
annotation('textbox', [0.6, 0.34, 0.25, 0.1], ...
    'String', sprintf("Backtracked distance (no hybrid): %0.0f nt\nBacktraced distance (full hybrid): %0.0f nt", ...
    params.backtrack, params.backtrack_hybrid), ...
    'FitBoxToText', 'on', ...
    'BackgroundColor', 'white', ...
    'EdgeColor', 'black', ...
    'FontSize', 12, ...
    'Color', default_MATLAB_colors(3));


% Theory curve assuming PolII does not continue to backtracks AND 
% The DNAP fully replicates
plot(t_ss_ext*(params.n_chase + params.backtrack_hybrid)*params.ssDNA.a...
    +t_hy_ext*(params.backtrack_hybrid)*params.hyDNA.a...
    +t_ds_ext*(params.n_chase)*params.dsDNA.a, t_force,...
    '--', "LineWidth",2, "Color",grey_color_08);

% Theory curve assuming PolII does not continue to backtrack AND 
% uses "fraction replicated" to fit data
plot(t_ss_ext*(params.n_chase + params.backtrack_hybrid + params.n_chase* (1- params.fraction_replicated))*params.ssDNA.a...
    +t_hy_ext*(params.backtrack_hybrid)*params.hyDNA.a...
    +t_ds_ext*(params.n_chase*params.fraction_replicated)*params.dsDNA.a, ...
    t_force,...
    '--', "LineWidth",2, "Color",default_MATLAB_colors(5));

annotation('textbox', [0.6, 0.49, 0.25, 0.1], ...
    'String', sprintf("DNA availble to replciate: %0.0f bp", ...
    params.n_chase), ...
    'FitBoxToText', 'on', ...
    'BackgroundColor', 'white', ...
    'EdgeColor', 'black', ...
    'FontSize', 12, ...
    'Color', default_MATLAB_colors(5));

% Plot the disruption force
scatter(max_force_ext,max_force,'o','r','LineWidth',3);
text(200,40,sprintf("MaxForce: %0.1f pN", max_force));

% Set all the fonts in the figure to size 14
set(findall(gcf,'-property','FontSize'),'FontSize',14)

% Optional: save plot
% saveas(fig1,params.file_title + "_Steps_MaxForce.png")

%% Measure the dsDNA produced during the force clamp.
% Create and clear figure 1 for plotting
fig2 = figure(2);
clf

% Select the step interval for analysis (Step 4)
k = 4;
slice = data.Step >= step_list(k) & data.Step < step_list(k+1);


% Extract extension and force data for this interval
ext_slice = ssext(slice);
force_slice = force(slice);

% Create upper subplot for extension analysis
ax1 = subplot(3,1,[1 2]);
title(params.trace_filename, "Interpreter","none")
hold on

% Decimate to ~100 Hz (reduce noise and data size)
params.Decimate = 10;
desimate_ext = Decimate_Average(ext_slice, params.Decimate);
desimate_force = Decimate_Average(force_slice, params.Decimate);

% Smooth with moving average to ~1 Hz
params.Smoothing = 20;
smooth_ext = movmean(desimate_ext, params.Smoothing);
smooth_force = movmean(desimate_force, params.Smoothing);
mean_force = mean(smooth_force);


% Find the start of the force clamp (force < 1.4 pN)
index_start = find(smooth_force < 1.4, 1, "first");

% Prepare time vectors
time2 = time(slice);
time2_decimate = Decimate_Average(time2, params.Decimate);
time2_smooth = movmean(time2_decimate, params.Smoothing);
time_start = time2_smooth(index_start);

% Calculate mean extension during early clamp (0.5–3 s after start)
ext_mean = mean(smooth_ext(time2_smooth > (time_start + 0.5) & ...
                            time2_smooth < (time_start + 3)));

% Estimate ssDNA extension at 1 pN
ss_Ext = ext_mean - ...
    (params.backtrack_hybrid * x_MMS(mean(smooth_force), params.kBT, ...
    params.hyDNA.Lp, params.hyDNA.K0) * params.hyDNA.a);
ssEXT_1pN = ss_Ext / (2 * params.n_chase + params.backtrack_hybrid);
ssEXT_1pN_calc = 0.0356;    % Empirical or theoretical estimate

% Convert Extension to dsDNA Estimates

% --- 100 Hz (decimated) dsDNA conversion ---
ext_chase_decimate = desimate_ext ...
    - params.n_chase * ssEXT_1pN_calc ...
    - params.backtrack_hybrid * ssEXT_1pN_calc ...
    - params.backtrack_hybrid * x_MMS(desimate_force, params.kBT, ...
        params.hyDNA.Lp, params.hyDNA.K0) * params.hyDNA.a;

decimate_dsDNAConverted = (ext_chase_decimate - params.n_chase * ssEXT_1pN_calc) ./ ...
    (-ssEXT_1pN_calc + x_MMS(desimate_force, params.kBT, ...
    params.dsDNA.Lp, params.dsDNA.K0) * params.dsDNA.a);

% --- 1 Hz (smoothed) dsDNA conversion ---
ext_chase_smooth = smooth_ext ...
    - params.n_chase * ssEXT_1pN_calc ...
    - params.backtrack_hybrid * ssEXT_1pN_calc ...
    - params.backtrack_hybrid * x_MMS(mean_force, params.kBT, ...
        params.hyDNA.Lp, params.hyDNA.K0) * params.hyDNA.a;

smooth_dsDNAConverted = (ext_chase_smooth - params.n_chase * ssEXT_1pN_calc) ./ ...
    (-ssEXT_1pN_calc + x_MMS(smooth_force, params.kBT, ...
    params.dsDNA.Lp, params.dsDNA.K0) * params.dsDNA.a);

% Prepare Time Vectors and Compute Replication Change
% Time vectors aligned to start of trace
time_slice = time(slice) - time(1);
time_slice_decimate = Decimate_Average(time_slice, params.Decimate);
time_slice_smooth = movmean(time_slice_decimate, params.Smoothing);

% Compute dsDNA at start (0.5–3 s after clamp begins)
dsConverted_start = mean(smooth_dsDNAConverted(time2_smooth > (time_start + 0.5) & ...
                                                time2_smooth < (time_start + 3)));

% Compute dsDNA at end (last 10 s of trace)
dsConverted_end = mean(smooth_dsDNAConverted(time2_smooth > (time2_smooth(end) - 10)));

% Total change in dsDNA (replication progress)
ext_change = dsConverted_end - dsConverted_start;

% Plot decimated dsDNA conversion (faint line)
plot(time_slice_decimate, ...
    decimate_dsDNAConverted, ...
    "Color", ...
    lighten_plot_color(default_MATLAB_colors(k),0.7));

% Plot smoothed dsDNA conversion (bold line)
plot(time_slice_smooth, ...
    smooth_dsDNAConverted, ...
    "Color", ...
    default_MATLAB_colors(k));


% Highlight key reference lines
yline(params.n_chase + dsConverted_start, '--', "Color", "Green", 'LineWidth', 3);          % Expected replication
yline(dsConverted_start, '--', "Color", "Black", 'LineWidth', 3);       % Start of clamp
yline(dsConverted_end, '--', "Color", "Black", 'LineWidth', 3);         % End of clamp


% Set y-axis limits
ylim([-500, 2500]);

% Add annotation text
text(60, 1800, ...
    sprintf("Force: %0.1f pN\nTranscribed Distance: %0.0f bp\n# of bp Replicated: %0.0f bp", ...
    mean_force, params.n_plus1 - params.n_chase, ext_change), ...
    "Color", default_MATLAB_colors(k));

text(180, params.n_chase+100, 'Expected', ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');

text(80, dsConverted_start, 'Data Start - BP Replicated', ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');

text(80, dsConverted_end, 'Data End - BP Replicated', ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');

% Plot Force Trace (Bottom Subplot)
ax2 = subplot(3,1,3);

% Plot decimated force
plot(time_slice_decimate, ...
     desimate_force, ...
     "Color", lighten_plot_color(default_MATLAB_colors(k), 0.7));
hold on

% Plot smoothed force
plot(time_slice_smooth, ...
     smooth_force, ...
     "Color", default_MATLAB_colors(k));

% Format axes
ylim([0, 2]);
ylabel("Force (pN)");
xlabel("Time (s)");

% Link x-axes of both subplots
linkaxes([ax1, ax2], 'x');

% Set font size for all elements
set(findall(gcf, '-property', 'FontSize'), 'FontSize', 14);

% Optional: Save the figure (uncomment to use)
% saveas(fig2, params.file_title + "_DNAPActivity.png");


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
    x_ds = x_MMS(force, params.kBT, params.hyDNA.Lp, params.hyDNA.K0) * params.hyDNA.a;
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
