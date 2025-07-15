function [position_smoothed,velocity] = calculate_velocity_position(time,position,width)
% Author: ji33
% Version: 1.0
% Date: 230411
%
%The velocity function calculates the velocity from the position vs time
%data.
% 
% Position is the position data from which to calculate the velocity
% Time is the time array
% Assumption 1: Time array is assumed to be uniformly sampled and increasing.
% Width is the standard deviation (sigma) of the Guassian weight function
% used in the fit. Units of the width parameter should match the units of
% the time array.

dt = time(2)-time(1); %Get the dt between data points.  Uses assumptions 1.
sample_factor = 3; %How many sigma on each side to include in the sampling region.  Default is three.
width_index = ceil(width*sample_factor/dt); %Get the half width in array index units
%Width_index must be at least 1.  This results in 3 data points being used
%for the velocity fit. 
if width_index < 1 
    width_index = 1;
end

%Generate the weight array to be used in fitting.
x_values = (-width_index:width_index)*dt; %X values to calculate the weights
weights = gaussian(x_values,0,width); %Array of weights

position_length = length(position); %Length of the input poisition data
velocity = zeros(size(position)); %Initialize the output array.
position_smoothed = zeros(size(position));

%%%
%Each point is indepent of the next.  I made this a parfor loop.  MATLAB
%takes a little while to launch the parrallel processing environment.  Then
%the calculation is more efficient.  This could be changed back to a
%regular loop if desired.
parfor i = 1:position_length
    %Slice a region of the position and time arrays for fitting.
    %start_index is the index to start the slicing.
    start_index = i - width_index;
    %Fix the slice start_index when near the start of the data
    if start_index < 1
        start_index = 1;
        weight_start_index = width_index - i + 2;
    else
        weight_start_index = 1;
    end
    
    %end_index is the index to end the slicing.
    end_index = i + width_index;
    %Fix the slice end_index when near the end of the data
    if end_index > position_length
        end_index = position_length;
        weight_end_index = width_index + 1 - (i - position_length);
    else
        weight_end_index = 2*width_index + 1;
    end
    
    %Send the sliced data to the MATLAB fit function for linear fit with
    %weights
    result = fit(time(start_index:end_index),position(start_index:end_index),'poly1','Weight',weights(weight_start_index:weight_end_index));
    %Extract the linear fit parameter and store it in the velocity variable
    velocity(i) = result.p1;
    position_smoothed(i) = result(time(i));
end

end

function y = gaussian(x,mu,sigma)
% Calculation a normal distribuition.    
y = exp(-((x-mu)/sigma).^2/2)./(sigma*sqrt(2*pi));
end