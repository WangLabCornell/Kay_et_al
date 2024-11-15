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
params.n_hybrid = 0; % shift size
params.n_last_aligned = 0; % last aligned location
max_force = [];
plot_shift = 100;

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

index = find(jindex > params.n_last_aligned,1,"first");    
force_after = force(index:end);
extension_after = ssext(index:end);
extension_hybrid = x_MMS(force_after,params.kBT, params.hyDNA.Lp, params.hyDNA.K0)*params.hyDNA.a*params.n_hybrid;
jindex_after = ((extension_after - extension_hybrid)./(x_FJC(force_after,params.kBT, params.ssDNA.Lp, params.ssDNA.K0)*params.ssDNA.a) + params.n_hybrid)/2;

TF = isempty(max_force);
if TF == 1
    max_force = max(force_after);
end


% Steps with extra backtracking
step_list = [3,6,];
step_names = ["Initial unzip"];

file_name = params.trace_filename;
file_title = extractBefore(file_name,'_Converted');
file_title = convertStringsToChars(file_title);
params.file_title = file_title;


fig2 = figure(2);
clf
plot(jindex(1:index+plot_shift),force(1:index+plot_shift), "LineWidth",plot_options.linewidth)
hold on
plot(jindex_after(plot_shift:end),force_after(plot_shift:end), "LineWidth",plot_options.linewidth)
plot(theory.j_index_, theory.Force_pN_, "LineWidth",plot_options.linewidth)

xline(params.n_plus1, '--')
xline(params.n_last_aligned, '--', 'color', 'red')


legend(["Unzip prior to polII", "Unzip after polII","Theory", "+1 Location", "Transcribed location"], "Location", "northeast")

message = sprintf("+1 location: %d\nTranscribed distance: %d\nHybrid size: %d\nDisruption Force: % 0.2f", params.n_plus1,params.n_transcribed,params.n_hybrid,max_force);
text(100,50,message)

xlim([-100, 2600])
ylim([0,60])
ylabel("Force (pN)")
xlabel("Number of basepairs unzipped")
title(params.file_title, "Interpreter","none")
% Set all the fonts in the figure to size 14
set(findall(gcf,'-property','FontSize'),'FontSize',16)
saveas(fig2,params.file_title + "_JIndex.png")

function color = default_MATLAB_colors(index)
    if index < 1
        index = 1;
    end
    index = mod(index-1,7)+1;
    plot_colors = {[0 0.4470 0.7410],[0.8500 0.3250 0.0980],[0.9290 0.6940 0.1250],[0.4940 0.1840 0.5560],[0.4660 0.6740 0.1880],[0.3010 0.7450 0.9330],[0.6350 0.0780 0.1840]};
    color = plot_colors{index};
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