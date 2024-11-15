clear

plot_options.linewidth = 2;

% ds, RNA-DNA hybrid, and ssDNA parameters
params.kBT = 4.09; % pNnm
params.dsDNA.Lp = 42; % nm
params.dsDNA.K0 = 1200; % pN
params.dsDNA.a = 0.338; % nm/bp

params.hyDNA.Lp = 50; % nm
params.hyDNA.K0 = 700; % pN
params.hyDNA.a = 0.3; % nm/bp

params.ssDNA.Lp = 0.848; % nm
params.ssDNA.K0 = 504.16; % pN
params.ssDNA.a = 0.546; % nm/bp

% Template specific information
params.n_plus1 = 0; % bp.  5' position of the RNA (HO 1596. CD 738)
params.n_arms = 0; % bp.  The number of bp in the arms (combined) (8330 bp of dsDNA)

params.smoothing = 50;

params.path = "/"; % folder name
params.trace_filename = ".dat"; % trace name
params.shift = 0; % shift param
params.stretch = 1; % stretch param
hybridshift = 0;

% Load data
data = readtable(params.path + params.trace_filename);
theory = readtable(".dat"); % Theoretically predicited unzipping file

theory.ssDNAExtension_nm_ = theory.extension_nm_ - x_MMS(theory.Force_pN_,params.kBT,params.dsDNA.Lp,params.dsDNA.K0)*params.dsDNA.a*params.n_arms;
theory.j_index_ = (theory.ssDNAExtension_nm_)./(x_FJC(theory.Force_pN_,params.kBT, params.ssDNA.Lp, params.ssDNA.K0)*params.ssDNA.a)/2;

params.n_transcribed = params.n_last_aligned - params.n_plus1;

% Extract out data from data file
time = data.Time;
time = time - time(1);
force = data.F_Scaled;
slice = force>5;
force = force(slice);
% Remove dsDNA y-arm extension from base extension to allow for shift and
% stretch alignment to the theory
ssext_raw = data.Extensionnm(slice) - x_MMS(force,params.kBT,params.dsDNA.Lp,params.dsDNA.K0)*params.dsDNA.a*params.n_arms;
ssext = params.stretch*(ssext_raw - params.shift);
% Conversion of "ssDNA" into number of base pairs unzipped - conversion
% only works for data preceeding encounter with protein
jindex = (ssext)./(x_FJC(force,params.kBT, params.ssDNA.Lp, params.ssDNA.K0)*params.ssDNA.a)/2;
step = data.Step;


% Steps with extra backtracking
step_list = [3,6,7,8,9,12];
step_names = [...
    "Initial unzip",...
    "Backtrack",...
    "Reduce force to 2 pN",...
    "Replication at 1 pN",...
    "Final unzip"];

file_name = params.trace_filename;
file_title = extractBefore(file_name,'_Converted');
file_title = convertStringsToChars(file_title);
params.file_title = file_title;

% Plot the unzipping with a shift
fig1 = figure(1)
clf

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot the force vs time to help identify steps.
subplot(3,1,1);
title(params.trace_filename, "Interpreter","none")
hold on
for k = 1:(length(step_list)-1)
    slice = data.Step >= step_list(k) & data.Step < step_list(k+1);
    plot(time(slice),force(slice));
end
ylim([-5,35])
ylabel("Force (pN)")
xlabel("Time (s)")
legend(step_names)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot the force vs extension
subplot(3,1,[2 3]);
hold on
for k = 1:(length(step_list)-1)
    slice = data.Step >= step_list(k) & data.Step < step_list(k+1);
    if k==1
        slice = slice & force > 5;
        % Find the position of the PolII first encountered on template
        [ext_pt,force_pt] = get_pt(ext(slice),force(slice),20); % get point from data
        params.n_chase = get_n_ss(ext_pt,force_pt, params)/2; 
    end

    if k==3
        [ext_pt,force_pt] = get_pt(ext(slice),force(slice),20); % get point from data
        params.backtrack = get_n_ss(ext_pt,force_pt, params)/2 - params.n_chase;

        delta_ext = ext_pt - 2 * params.n_chase * x_FJC(force_pt,params.kBT, params.ssDNA.Lp, params.ssDNA.K0)*params.ssDNA.a;
        params.backtrack_hybrid = get_n_hybrid(delta_ext,force_pt, params)/2;
%         params.n_backtrack_hybrid = get_n_hybrid(ext(slice),force(slice),20, params)
        
    end

    if k==5
        slice = slice & force > 0.5;
        slice = slice & ext < 1800;
        
        [ext_pt,force_pt] = get_pt(ext(slice),force(slice),20); % get point from data

       
        params.post_replication_ssDNA = get_n_ss(ext_pt,force_pt, params)/2; 
            
        %Fit the data assuming no further backtracking to see what fraction
        %is replicated.
        delta_ext = ext_pt...
            - params.n_chase * x_FJC(force_pt,params.kBT, params.ssDNA.Lp, params.ssDNA.K0)*params.ssDNA.a...
            - params.backtrack_hybrid * x_FJC(force_pt,params.kBT, params.ssDNA.Lp, params.ssDNA.K0)*params.ssDNA.a...
            - params.backtrack_hybrid * x_MMS(force_pt,params.kBT, params.hyDNA.Lp, params.hyDNA.K0)*params.hyDNA.a;
        params.fraction_replicated = get_n_hybrid_fraction(delta_ext,force_pt, params.n_chase, params);
        params.fraction_replicated_number = get_n_hybrid_number(delta_ext,force_pt, params.n_chase, params);
        


        delta_ext = ext_pt...
            - params.n_plus1 * x_FJC(force_pt,params.kBT, params.ssDNA.Lp, params.ssDNA.K0)*params.ssDNA.a...
            - (params.n_plus1 - params.n_chase) * x_MMS(force_pt,params.kBT, params.dsDNA.Lp, params.dsDNA.K0)*params.dsDNA.a;
        params.fraction_replicated_full_backtrack = get_n_hybrid_fraction(delta_ext,force_pt, params.n_chase, params);
        
        max_force = max(force(slice));
        max_index = find(force(slice) == max_force);
        ext2 = ext(slice);
        max_force_ext = ext2(max_index);

    end
    
    plot(ext(slice),force(slice));
end


% Find the position of the PolII
xlim([-100, 3000])
xline(0,'--',"LineWidth",2)
ylim([0,max_force+10])
xlabel("Extension with arms removed (nm)")
ylabel("Force (pN)")

% Add theory curves to unzipping plot
grey_color = [0.5,0.5,0.5];

% Plot the unzipping theory curve
plot(theory.ssDNAExtension_nm_,theory.Force_pN_,"LineWidth",2, "Color",grey_color);
hold on
t_force = 0:0.1:40; % Array of force values used to plot theory curves
t_ss_ext = x_FJC(t_force,params.kBT, params.ssDNA.Lp, params.ssDNA.K0); % ext/Lc (normalized extension) of ssDNA used to plot theory curves
t_hy_ext = x_MMS(t_force,params.kBT, params.hyDNA.Lp, params.hyDNA.K0); % ext/Lc (normalized extension) of dsDNA used to plot theory curves
t_ds_ext = x_MMS(t_force,params.kBT, params.dsDNA.Lp, params.dsDNA.K0); % ext/Lc (normalized extension) of dsDNA used to plot theory curves


% Plot the +1 theory curve
plot(t_ss_ext*params.n_plus1*2*params.ssDNA.a, t_force,'--', "LineWidth",2, "Color",grey_color);

% Plot the position chased PolII
plot(t_ss_ext*params.n_chase*params.ssDNA.a*2, t_force,'--', "LineWidth",2, "Color",default_MATLAB_colors(1));
text(x_FJC(30,params.kBT, params.ssDNA.Lp, params.ssDNA.K0)*params.ssDNA.a*params.n_chase*2+20,50,...
    sprintf("Chased position: %0.0f bp\nTranscribed distance: %0.0f nt", params.n_chase, params.n_plus1 - params.n_chase),...
    "Color", default_MATLAB_colors(1));
% Plot the position of the backtracked PolII using two models:
% 1) The RNA does not hybridize
% 2) The RNA fully hybridizes.
plot(t_ss_ext*(params.n_chase *2 + params.backtrack * 2)*params.ssDNA.a, t_force,'--', "LineWidth",2, "Color",default_MATLAB_colors(3));
plot(t_ss_ext*(params.n_chase *2 + params.backtrack_hybrid)*params.ssDNA.a...
    +t_hy_ext*(params.backtrack_hybrid)*params.hyDNA.a, t_force,'--', "LineWidth",2, "Color",lighten_plot_color(default_MATLAB_colors(3),0.5));
text(x_FJC(25,params.kBT, params.ssDNA.Lp, params.ssDNA.K0)*params.ssDNA.a*(params.n_chase *2 + params.backtrack * 2)+20,45,...
    sprintf("Backtracked distance (no hybrid): %0.0f nt\nBacktraced distance (full hybrid): %0.0f nt", params.backtrack, params.backtrack_hybrid),...
    "Color", default_MATLAB_colors(3));

% Theory curve assuming PolII does not continue to backtracks AND 
% The DNAP fully replicates
plot(t_ss_ext*(params.n_chase + params.backtrack_hybrid)*params.ssDNA.a...
    +t_hy_ext*(params.backtrack_hybrid)*params.hyDNA.a...
    +t_ds_ext*(params.n_chase)*params.dsDNA.a, t_force,...
    '--', "LineWidth",2, "Color",grey_color);
% +t_ds_ext*(params.n_chase + params.backtrack_hybrid)*params.dsDNA.a, t_force,...

% Theory curve assuming PolII does not continue to backtrack AND 
% uses "fraction replicated" to fit data
plot(t_ss_ext*(params.n_chase + params.backtrack_hybrid + params.n_chase* (1- params.fraction_replicated))*params.ssDNA.a...
    +t_hy_ext*(params.backtrack_hybrid)*params.hyDNA.a...
    +t_ds_ext*(params.n_chase*params.fraction_replicated)*params.dsDNA.a, t_force,...
    '--', "LineWidth",2, "Color","Green");
text(1500,10,...
    sprintf("fraction replicated: %0.0f bp\nfraction replicated(2): %0.0f bp\nDNA availble to replciate: %0.0f bp",...
    params.n_chase*params.fraction_replicated,params.fraction_replicated_number,params.n_chase),...
    "Color","Green");

% Plot the disruption force
scatter(max_force_ext,max_force,'o','r','LineWidth',3);
text(200,40,sprintf("MaxForce: %0.1f pN", max_force));

text()

% Set all the fonts in the figure to size 14
set(findall(gcf,'-property','FontSize'),'FontSize',14)

% Measure the dsDNA produced during the force clamp.

fig2 = figure(2);
clf
k=4;
slice = data.Step >= step_list(k) & data.Step < step_list(k+1);
% slice = slice & force > 0.5;
ext_slice = ext(slice);
force_slice = force(slice);
ax1 = subplot(3,1,[1 2]);
title(params.trace_filename, "Interpreter","none")
hold on

% decimate to 100 hz
params.Decimate = 10;
desimate_ext = Decimate_Average(ext_slice,params.Decimate);
desimate_force = Decimate_Average(force_slice,params.Decimate);

% smooth with movmean to 1 Hz
params.Smoothing = 20;
smooth_ext = movmean(desimate_ext,params.Smoothing);
smooth_force = movmean(desimate_force,params.Smoothing);
mean_force = mean(smooth_force);

index_start = find(smooth_force < 1.4, 1, "first");
time2 = time(slice);
time2_decimate = Decimate_Average(time2,params.Decimate);
time2_smooth = movmean(time2_decimate,params.Smoothing);
time_start = time2_smooth(index_start);

ext_mean = mean(smooth_ext(time2_smooth > (time_start+0.5) & time2_smooth < (time_start + 3))) ;
ss_Ext = ext_mean - (params.backtrack_hybrid * x_MMS(mean(smooth_force),params.kBT, params.hyDNA.Lp, params.hyDNA.K0)*params.hyDNA.a);
ssEXT_1pN = ss_Ext / (2*params.n_chase + params.backtrack_hybrid);
ssEXT_1pN_calc = 0.0356;

%100 Hz - ds conversion
ext_chase_decimate = desimate_ext...
    - params.n_chase * ssEXT_1pN_calc...
    - params.backtrack_hybrid * ssEXT_1pN_calc...
    - params.backtrack_hybrid * x_MMS(desimate_force,params.kBT, params.hyDNA.Lp, params.hyDNA.K0)*params.hyDNA.a;

decimate_dsDNAConverted = (ext_chase_decimate - params.n_chase*(ssEXT_1pN_calc))./(-ssEXT_1pN_calc+...
    x_MMS(desimate_force,params.kBT, params.dsDNA.Lp, params.dsDNA.K0)*params.dsDNA.a);

%1 Hz - ds conversion
ext_chase_smooth = smooth_ext...
    - params.n_chase * ssEXT_1pN_calc...
    - params.backtrack_hybrid * ssEXT_1pN_calc...
    - params.backtrack_hybrid * x_MMS(mean_force,params.kBT, params.hyDNA.Lp, params.hyDNA.K0)*params.hyDNA.a;

smooth_dsDNAConverted = (ext_chase_smooth - params.n_chase*(ssEXT_1pN_calc))./(-ssEXT_1pN_calc+...
    x_MMS(smooth_force,params.kBT, params.dsDNA.Lp, params.dsDNA.K0)*params.dsDNA.a);


time_slice = time(slice);
time_slice = time_slice - time_slice(1);
time_slice_decimate = Decimate_Average(time_slice,params.Decimate);
time_slice_smooth = movmean(time_slice_decimate,params.Smoothing);

% The START extension is the first 0.5 to 3 seconds of the trace
dsConverted_start = mean(smooth_dsDNAConverted(time2_smooth > (time_start+0.5) & time2_smooth < (time_start + 3))) ;

% The end extension is the last 10 seconds of the trace
dsConverted_end = mean(smooth_dsDNAConverted(time2_smooth > (time2_smooth(end) - 10))); % Find the extension in the 

ext_change = dsConverted_end - dsConverted_start;

plot(time_slice_decimate,decimate_dsDNAConverted, "Color", lighten_plot_color(default_MATLAB_colors(k),0.7));
plot(time_slice_smooth,smooth_dsDNAConverted, "Color", default_MATLAB_colors(k));
yline(params.n_chase, '--', "Color", "Green",'LineWidth',3)
yline(dsConverted_start, '--', "Color", "Black",'LineWidth',3)
yline(dsConverted_end, '--', "Color", "Black",'LineWidth',3)
ylim([-500,2500])
text(60,-250,...
    sprintf("Force: %0.1f pN\nTranscribed Distance %0.0f bp\n# of bp Replicated: %0.0f bp",mean_force,params.n_plus1 - params.n_chase,ext_change),...
    "Color", default_MATLAB_colors(k))
text(20, params.n_chase, 'Expected', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
text(60, dsConverted_start, 'Data Start - BP Replcated', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
text(60, dsConverted_end, 'Data End - BP Replciated', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');

ax2 = subplot(3,1,3);
plot(time_slice_decimate,desimate_force, "Color", lighten_plot_color(default_MATLAB_colors(k),0.7));
hold on
plot(time_slice_smooth,smooth_force, "Color", default_MATLAB_colors(k));

ylim([0,2])
ylabel("Force (pN)")
xlabel("Time (s)")

linkaxes([ax1,ax2],'x')
set(findall(gcf,'-property','FontSize'),'FontSize',14)
saveas(fig2,params.file_title + "_DNAPActivity.png")

time_slice_smooth = time_slice_smooth';
smooth_dsDNAConverted = smooth_dsDNAConverted';
smooth_ext = smooth_ext';


index1 = find(time_slice_smooth>0.5,1);
index2 = find(time_slice_smooth>3,1);
zero = mean(smooth_dsDNAConverted(index1:index2));
smooth_dsDNAConverted2 = smooth_dsDNAConverted - zero;

Lagzero_ext = (zero * (-ssEXT_1pN_calc+...
    x_MMS(mean(force_slice),params.kBT, params.dsDNA.Lp, params.dsDNA.K0)*params.dsDNA.a)) +...
    params.n_chase*(ssEXT_1pN_calc);
zero_ext = Lagzero_ext...
    + params.n_chase * ssEXT_1pN_calc...
    + params.backtrack_hybrid * ssEXT_1pN_calc...
    + params.backtrack_hybrid * x_MMS(mean(force_slice),params.kBT, params.hyDNA.Lp, params.hyDNA.K0)*params.hyDNA.a;

smooth_ext2 = smooth_ext - zero_ext;
fig3 = figure(3);
clf
plot(time_slice_smooth,smooth_dsDNAConverted2)
hold on
yline(params.n_chase, '--', "Color", "Green",'LineWidth',3)
yline(dsConverted_start - zero, '--', "Color", "Black",'LineWidth',3)
yline(dsConverted_end - zero, '--', "Color", "Black",'LineWidth',3)

annotation('textbox', [0.2, 0.8, 0.1, 0.1], 'String', sprintf("Force: %0.1f pN\nTranscribed Distance %0.0f bp\n# of bp Replicated: %0.0f bp",mean_force,params.n_plus1 - params.n_chase,ext_change));
text(20, params.n_chase, 'Expected', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
text(60, dsConverted_start - zero, 'Data Start - BP Replcated', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
text(60, dsConverted_end - zero, 'Data End - BP Replciated', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');

yyaxis right
plot(time_slice_smooth,smooth_ext2)

title(params.trace_filename, "Interpreter","none")

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ext_pt, force_pt] = get_pt(ext,force,f)
    %get data around the force
    ext_pt = mean(ext(force>(f-0.5) & force<(f+0.5)));
    force_pt = mean(force(force>(f-0.5) & force<(f+0.5)));
end

% function y = ssEXT_1pN_calc(Force)
% y = -0.0165*Force + 0.0499;
% end


function n = get_n_ss(ext,force, params)
    % Get the number of nucleotides in a strand of ssDNA
    n = ext./(x_FJC(force,params.kBT, params.ssDNA.Lp, params.ssDNA.K0)*params.ssDNA.a);
end

function n = get_n_ds(ext,force, params)
    % Get the number of nucleotides in a strand of dsDNA
    n = ext./(x_MMS(force,params.kBT, params.dsDNA.Lp, params.dsDNA.K0)*params.dsDNA.a);
end

function n = get_n_hybrid(ext,force, params)
    % Get the TOTAL number of nucleotides in a strand that is half ssDNA and half
    % dsDNA
    % Example: ---=== would have n = 6
    n = 2 * ext/(x_FJC(force,params.kBT, params.ssDNA.Lp, params.ssDNA.K0)*params.ssDNA.a...
        + x_MMS(force,params.kBT, params.hyDNA.Lp, params.hyDNA.K0)*params.hyDNA.a);
end

function f = get_n_hybrid_fraction(ext,force, n, params)
    % Returns the fraction f of a DNA strand that is ss given the total
    % length is n.
    % Example:  --=== with n = 5, f = 0.6
    x_ss = x_FJC(force,params.kBT, params.ssDNA.Lp, params.ssDNA.K0)*params.ssDNA.a; %ssDNA length/bp at force
    x_ds = x_MMS(force,params.kBT, params.dsDNA.Lp, params.dsDNA.K0)*params.dsDNA.a; %dsDNA length/bp at force
    f = (ext/n - x_ss)./(x_ds - x_ss);
end

function f = get_n_hybrid_number(ext,force, n, params)
    % Returns the fraction f of a DNA strand that is ss given the total
    % length is n.
    % Example:  --=== with n = 5, f = 0.6
    x_ss = x_FJC(force,params.kBT, params.ssDNA.Lp, params.ssDNA.K0)*params.ssDNA.a; %ssDNA length/bp at force
    x_ds = x_MMS(force,params.kBT, params.dsDNA.Lp, params.dsDNA.K0)*params.dsDNA.a; %dsDNA length/bp at force
    f = (ext - n*x_ss)./(x_ds - x_ss);
end

function Smooth = Decimate_Average(data,SmoothFactor)
    % Calculate the number of complete chunks of size N
    numChunks = floor(length(data) / SmoothFactor);

    % Reshape only the complete chunks
    reshapedData = reshape(data(1:numChunks*SmoothFactor), SmoothFactor, []);


    % Calculate the mean of each column
    Smooth = mean(reshapedData);

    % If the length of data is not a multiple of N, handle the remaining points
    if mod(length(data), SmoothFactor) ~= 0
        remainingPoints = data(end-mod(length(data), SmoothFactor)+1:end);
        Smooth = [Smooth, mean(remainingPoints)];
    end
end


function plot_guideline(ext,force,n_ss,n_ds,params)
    

end

function x = x_FJC(F,kBT,Lp,K0)

% def extension_FJC(force,Lp,K0,kbT):
%     if isinstance(force, np.ndarray):
%         return np.array([extension_FJC_eval(f,Lp,K0,kbT) for f in force])
%     if isinstance(force, int) or isinstance(force, float):
%         return extension_FJC_eval(force,Lp,K0,kbT)
%     else:
%         raise ValueError("Input must be numeric")
% 
% def extension_FJC_eval(force,Lp,K0,kbT):
%     '''
%     Formula is from MDW, et al, Biophys. J. 72, 1335 (1997).
%     See p. 1342.
%     '''
%     if force:
%         A = 2*Lp/kbT*force
%         return (1/np.tanh(A) - 1/A ) * (force/K0 + 1)
%     else:
%         return 0

    A = 2*Lp/kBT*F;
    x = (1./tanh(A) - 1./A) .* (F/K0 + 1);
end

function x = x_MMS(F,kBT,Lp,K0)
    % Formula is from MDW, et al, Biophys. J. 72, 1335 (1997).
    % See p. 1342.

    % Returns the normalized extension x := extension / L_0
    c = K0*Lp/kBT;
    f = F/K0;
    
    % a0 + a1 * x + a2 * x**2 + x**3 = 0
    a0 = 1/4-(1/4+(c+1)*f).*((1+f).^2);
    a1 = (1+f).*((1+f)+1/2+2*(c+1)*f);
    a2 = -(2*(1+f)+1/4+(c+1)*f);

    x = zeros(size(a0));
    for i = 1:length(x)
        r = roots([1,a2(i),a1(i),a0(i)]);
        % This function has three roots for force < 0.176 pN and then three real roots after that point
        % The correct root is always the one with the smallest real part.
        x(i) = min(real(r));
    end
end

function plot_color_light = lighten_plot_color(plot_color,percent)
    %Lighten plot color by mixing with white.
    plot_color_light = plot_color + percent*(1-plot_color);
end

function color = default_MATLAB_colors(index)
    if index < 1
        index = 1;
    end
    index = mod(index-1,7)+1;
    plot_colors = {[0 0.4470 0.7410],[0.8500 0.3250 0.0980],[0.9290 0.6940 0.1250],[0.4940 0.1840 0.5560],[0.4660 0.6740 0.1880],[0.3010 0.7450 0.9330],[0.6350 0.0780 0.1840]};
    color = plot_colors{index};
end
