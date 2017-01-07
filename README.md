## Summary

Add clickable interactivity to an existing scatter, line or image plot. The user supplies additional data specific to each data point on the reference plot. This data is displayed in a separate suplot ('data window') when the user clicks points on the reference plot. Multiple 'windows' can be added to the same plot.

## Installation

Download, unzip (or clone) and add the folder to Matlab's search path (right-clcik, 'add to path') and you're good to go!

## Example 

**See the attached test_plot_interaction.m for more examples**

Visualise mu*sigma of a 2d gaussian with sigma values varying on the x and mu on the y axes
close all
clearvars

murange = 0:0.02:1;
sigrange = 0:0.01:1;

[sig,mu] = meshgrid(sigrange,murange);
ellipt = mu.*sig;

imagesc(ellipt)
xlabel('Sig. 1')
ylabel('Sig. 2')
title('Ellipticity of gaussian')

xrange = 0:0.01:1;
interaction_data = cell(length(sigrange),length(murange));
for i = 1:length(sigrange)
    for j = 1:length(murange)
        interaction_data{i,j} = [xrange;normpdf(xrange,murange(j),sigrange(i))];
    end
end

add_plot_interaction_new(gca,interaction_data,'o-')
xlabel('X')
ylabel('Y')