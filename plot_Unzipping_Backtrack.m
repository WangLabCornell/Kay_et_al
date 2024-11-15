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
shift = 0;

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

broke_through = 20; 

% Steps with extra backtracking
step_list = [3,6,7,8,12];
step_names = [...
    "Initial unzip",...
    "Backtrack",...
    "Reduce force to 1 pN",...
    "Final unzip"];

file_name = params.trace_filename;
file_title = extractBefore(file_name,'_Converted');
file_title = convertStringsToChars(file_title);

% Plot the unzipping with a shift
figure(1);
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
ylim([-5,50])
ylabel("Force (pN)")
xlabel("Time (s)")
legend(step_names,'Location','eastoutside')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot the force vs extension
subplot(3,1,[2 3]);
hold on
for k = 1:(length(step_list)-1)
    slice = data.Step >= step_list(k) & data.Step < step_list(k+1);
    if k==1
        slice = slice & force > 5;
        % Find the position of the PolII first encountered on template
        [ext_pt1,force_pt] = get_pt(ext(slice),force(slice),22); % get point from data
        params.n_chase = get_n_ss(ext_pt1,force_pt, params)/2; 
    end

    if k==3
        [ext_pt,force_pt] = get_pt(ext(slice),force(slice),broke_through); % get point from data
        params.backtrack = get_n_ss(ext_pt,force_pt, params)/2 - params.n_chase;
        
        params.backtracknm = ext_pt - ext_pt1;

        delta_ext = ext_pt - 2 * params.n_chase * x_FJC(force_pt,params.kBT, params.ssDNA.Lp, params.ssDNA.K0)*params.ssDNA.a;
        params.backtrack_hybrid = get_n_hybrid(delta_ext,force_pt, params)/2;
        %params.n_backtrack_hybrid = get_n_hybrid(ext(slice),force(slice),20, params)
        
    end

    if k==4
        slice = slice & force > 0;
        [ext_pt2,force_pt2] = get_pt(ext(slice),force(slice),20); % get point from data

        params.post_replication_nm = ext_pt - ext_pt2;
        params.post_replication_ssDNA = get_n_ss(ext_pt2,force_pt2, params)/2; 
            
    end
    
    plot(ext(slice),force(slice));
end


% Find the position of the PolII

xlim([-100, 3000])
xline(0,'--',"LineWidth",2)
ylim([0,50])
xlabel("Extension with arms removed (nm)")
ylabel("Force (pN)")

% Add theory curves to unzipping plot
grey_color = [0.5,0.5,0.5];

% Plot the unzipping theory curve
plot(theory.ssDNAExtension_nm_,theory.Force_pN_,"LineWidth",2, "Color",grey_color);
hold on
t_force = 0:0.1:40; % Array of force values used to plot theory curves
t_ss_ext = x_FJC(t_force,params.kBT, params.ssDNA.Lp, params.ssDNA.K0); % ext/Lc (normalized extension) of ssDNA used to plot theory curves
t_ds_ext = x_MMS(t_force,params.kBT, params.dsDNA.Lp, params.dsDNA.K0); % ext/Lc (normalized extension) of dsDNA used to plot theory curves
t_ds_RD_ext = x_MMS(t_force,params.kBT, params.dsRDNA.Lp, params.dsRDNA.K0); % ext/Lc (normalized extension) of dsDNA used to plot theory curves
% Plot the +1 theory curve
plot(t_ss_ext*params.n_plus1*2*params.ssDNA.a, t_force,'--', "LineWidth",2, "Color",grey_color);


% Plot the position chased PolII
plot(t_ss_ext*params.n_chase*params.ssDNA.a*2, t_force,'--', "LineWidth",2, "Color",default_MATLAB_colors(1));
text(x_FJC(45,params.kBT, params.ssDNA.Lp, params.ssDNA.K0)*params.ssDNA.a*params.n_chase*2+20,45,...
    sprintf("Chased position: %0.0f bp\nTranscribed distance: %0.0f nt", params.n_chase, params.n_plus1 - params.n_chase),...
    "Color", default_MATLAB_colors(1));
% Plot the position of the backtracked PolII using two models:
% 1) The RNA does not hybridize
% 2) The RNA fully hybridizes.
plot(t_ss_ext*(params.n_chase *2 + params.backtrack * 2)*params.ssDNA.a, t_force,'--', "LineWidth",2, "Color",default_MATLAB_colors(3));
plot(t_ss_ext*(params.n_chase *2 + params.backtrack_hybrid)*params.ssDNA.a...
    +t_ds_RD_ext*(params.backtrack_hybrid)*params.dsRDNA.a, t_force,'--', "LineWidth",2, "Color",lighten_plot_color(default_MATLAB_colors(3),0.5));
text(x_FJC(40,params.kBT, params.ssDNA.Lp, params.ssDNA.K0)*params.ssDNA.a*(params.n_chase *2 + params.backtrack * 2)+20,40,...
    sprintf("Backtracked Distance: %0.1f nm\nBacktracked distance (no hybrid): %0.0f nt\nBacktraced distance (full hybrid): %0.0f nt", params.backtracknm,params.backtrack, params.backtrack_hybrid),...
    "Color", default_MATLAB_colors(3));

plot(t_ss_ext*(params.post_replication_ssDNA*2)*params.ssDNA.a, t_force,'--', "LineWidth",2, "Color",lighten_plot_color(default_MATLAB_colors(6),0.5));
text(x_FJC(50,params.kBT, params.ssDNA.Lp, params.ssDNA.K0)*params.ssDNA.a*(params.post_replication_ssDNA *2)+20,50,...
    sprintf("2nd Unzip position: %0.0f bp\n Transcription from +1: %0.0f bp\n Transcription from BackTrack: %0.1f nm", params.post_replication_ssDNA,params.n_plus1 - params.post_replication_ssDNA,params.post_replication_nm),...
    "Color", default_MATLAB_colors(6));

text()

% Set all the fonts in the figure to size 14
set(findall(gcf,'-property','FontSize'),'FontSize',12)

figure(2)
clf
k=4;
slice = data.Step >= step_list(k) & data.Step < step_list(k+1);
slice = slice & force > 1;
hold on
plot(ext(slice),force(slice), "Color", default_MATLAB_colors(6));
title(params.trace_filename, "Interpreter","none")

xlim([-100, 3000])
xline(0,'--',"LineWidth",2)
ylim([0,50])
xlabel("Extension with arms removed (nm)")
ylabel("Force (pN)")


% Add theory curve to unzipping plot
grey_color = [0.8,0.8,0.8];
plot(theory.ssDNAExtension_nm_,theory.Force_pN_,"LineWidth",2, "Color",grey_color);
hold on
t_force = 0:0.1:40; % Array of force values used to plot theory curves
t_ss_ext = x_FJC(t_force,params.kBT, params.ssDNA.Lp, params.ssDNA.K0); % ext/Lc (normalized extension) of ssDNA used to plot theory curves
t_ds_ext = x_MMS(t_force,params.kBT, params.dsDNA.Lp, params.dsDNA.K0); % ext/Lc (normalized extension) of dsDNA used to plot theory curves
plot(t_ss_ext*params.n_plus1*2*params.ssDNA.a, t_force,'--', "LineWidth",2, "Color",grey_color);

grey_color = [0.5,0.5,0.5];


plot(theory.ssDNAExtension_nm_- shift,theory.Force_pN_,"LineWidth",2, "Color",grey_color);
%plot(t_ss_ext*params.n_plus1*2*params.ssDNA.a-shift, t_force,'--', "LineWidth",2, "Color",grey_color);

params.shift_unzip = shift/(x_FJC(15,params.kBT, params.ssDNA.Lp, params.ssDNA.K0)*params.ssDNA.a-x_MMS(15,params.kBT, params.dsRDNA.Lp, params.dsRDNA.K0)*params.dsRDNA.a)

set(findall(gcf,'-property','FontSize'),'FontSize',14)
legend('Unzip After Backtrack','Ext = 0','Unzipping Theory','Transcription Start Site - Unshifted','Shifted Unzipping Theory')
text(50,45,sprintf("Transcribed Distance: %0.0f bp\nShift Size: %0.0f bp", params.n_plus1 - params.n_chase, params.shift_unzip));

% Measure the hybrid produced during the force clamp.
figure(3);
clf
k=2;
slice = data.Step >= step_list(k) & data.Step < step_list(k+1);
%slice = slice & force > 20;
ax1 = subplot(3,1,[1 2]);
title(params.trace_filename, "Interpreter","none")
hold on

smoothed_ext = movmean(ext(slice),50);
smoothed_force = movmean(force(slice),50);
plot(time(slice),ext(slice), "Color", lighten_plot_color(default_MATLAB_colors(k),0.7));
plot(time(slice),smoothed_ext, "Color", default_MATLAB_colors(k));

%index_start = find(smoothed_force > 22, 1, "first");
index_start = 1;
time2 = time(slice);
time_start = time2(index_start);
% The start extension between 1 and 3 seconds after the force drops
ext_start = min(smoothed_ext) ;

% The end extension is the last 2 seconds of the trace
ext_end = max(smoothed_ext); % Find the extension in the 


mean_force = mean(smoothed_force(time2>time_start));

% yline(ext_start, '--', "Color", "Black")
% yline(ext_end, '--', "Color", "Black")
%ylim([0,900])
ylabel("Extension with arms removed (nm)")
xlabel("Time (s)")
NW = [min(xlim) max(ylim)]+[diff(xlim) -diff(ylim)]*0.05;
text(NW(1),NW(2),sprintf("Backtracked distance (no hybrid): %0.0f bp\nBacktraced distance (full hybrid): %0.0f bp\n Backtrack Extension Max - Min: %0.0f nm", params.backtrack, params.backtrack_hybrid,ext_end - ext_start))


ext_change = ext_end - ext_start;

ax2 = subplot(3,1,3);
plot(time(slice),force(slice), "Color", lighten_plot_color(default_MATLAB_colors(k),0.7));
plot(time(slice),smoothed_force, "Color", default_MATLAB_colors(k));
%ylim([-0.5,2])
ylabel("Force (pN)")
xlabel("Time (s)")

linkaxes([ax1,ax2],'x')
set(findall(gcf,'-property','FontSize'),'FontSize',14)


figure(4)
clf
k=3;
slice = data.Step >= step_list(k) & data.Step < step_list(k+1);
slice = slice & force > 0;
hold on
plot(ext(slice),force(slice), "Color", default_MATLAB_colors(6));
title(params.trace_filename, "Interpreter","none")

xlim([-1000, 3000])
xline(0,'--',"LineWidth",2)
ylim([0,50])
xlabel("Extension with arms removed (nm)")
ylabel("Force (pN)")

% Add theory curve to unzipping plot
grey_color = [0.8,0.8,0.8];
plot(theory.ssDNAExtension_nm_,theory.Force_pN_,"LineWidth",2, "Color",grey_color);
hold on
t_force = 0:0.1:40; % Array of force values used to plot theory curves
t_ss_ext = x_FJC(t_force,params.kBT, params.ssDNA.Lp, params.ssDNA.K0); % ext/Lc (normalized extension) of ssDNA used to plot theory curves
t_ds_ext = x_MMS(t_force,params.kBT, params.dsDNA.Lp, params.dsDNA.K0); % ext/Lc (normalized extension) of dsDNA used to plot theory curves
plot(t_ss_ext*params.n_plus1*2*params.ssDNA.a, t_force,'--', "LineWidth",2, "Color",grey_color);


% Plot the position of the backtracked PolII using two models:
% 1) The RNA does not hybridize
% 2) The RNA fully hybridizes.
plot(t_ss_ext*(params.n_chase *2 + params.backtrack * 2)*params.ssDNA.a, t_force,'--', "LineWidth",2, "Color",default_MATLAB_colors(3));
plot(t_ss_ext*(params.n_chase *2 + params.backtrack_hybrid)*params.ssDNA.a...
    +t_ds_RD_ext*(params.backtrack_hybrid)*params.dsRDNA.a, t_force,'--', "LineWidth",2, "Color",lighten_plot_color(default_MATLAB_colors(3),0.5));

force_dip] = ginput(1);

set(findall(gcf,'-property','FontSize'),'FontSize',14)

NW = [min(xlim) max(ylim)]+[diff(xlim) -diff(ylim)]*0.05;
text(NW(1),NW(2),sprintf("Extension Change: %0.2f nm\nBacktracked distance (no hybrid): %0.2f bp\n Force Dip Value: %0.2f pN", params.backtracknm, params.backtrack,force_dip))

figure(1)
text(NW(1),NW(2),sprintf("Force Dip Value: %0.2f pN",force_dip))
set(findall(gcf,'-property','FontSize'),'FontSize',12)

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
        + x_MMS(force,params.kBT, params.dsDNA.Lp, params.dsDNA.K0)*params.dsDNA.a);
end

function f = get_n_hybrid_fraction(ext,force, n, params)
    % Returns the fraction f of a DNA strand that is ss given the total
    % length is n.
    % Example:  --=== with n = 5, f = 0.6
    x_ss = x_FJC(force,params.kBT, params.ssDNA.Lp, params.ssDNA.K0)*params.ssDNA.a; %ssDNA length/bp at force
    x_ds = x_MMS(force,params.kBT, params.dsDNA.Lp, params.dsDNA.K0)*params.dsDNA.a; %dsDNA length/bp at force
    f = (ext/n - x_ss)./(x_ds - x_ss);
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