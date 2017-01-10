%% For example, generate a random set of mu and sigma values and plot a
% scatter plot. The interaction data will show plots of a gaussian
% distribution with the parameters corresponding to the datapoints on the
% left.
clearvars
close all
plot_data = rand(10,2); plot_data(:,1) = plot_data(:,1) - 0.5;
plot(plot_data(:,1),plot_data(:,2),'o')
xlabel('mu'); ylabel('sig')

 x = -1:0.01:1;

for i = 1:size(plot_data,1)
    interaction_data{i} = [x(:),reshape(normpdf(x,plot_data(i,1),plot_data(i,2)),[],1)];
end

add_plot_interaction(gca,interaction_data)

xlabel('Interaction Term X')
ylabel('Interaction Term Y')

%% Visualise mu*sigma of a 2d gaussian with sigma values varying on the x and mu on the y axes
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

add_plot_interaction(gca,interaction_data)
xlabel('X')
ylabel('Y')

%% Visualise ellipticity of a 2d gaussian varying sigma 1 and sigma 2
close all
clearvars

mu = 0.5;

sigrange1 = 0:0.01:1;
sigrange2 = 0:0.01:2;

[sig1,sig2] = meshgrid(sigrange1,sigrange2);
ellipt = sig1.*sig2;

imagesc(ellipt)
xlabel('Sig. 1')
ylabel('Sig. 2')
title('Ellipticity of gaussian')

xrange = 0:0.01:1;
yrange = 0:0.01:2;
[xm,ym] = meshgrid(xrange,yrange);

interaction_data = cell(length(sigrange1),length(sigrange2));
for i = 1:length(sigrange1)
    for j = 1:length(sigrange2)
        interaction_data{i,j} = normpdf(xm,mu,sigrange1(i)).*normpdf(ym,mu,sigrange2(j));
    end
end

add_plot_interaction(gca,interaction_data)
xlabel('X')
ylabel('Y')

%% Add multiple interaction windows
clearvars
close all
plot_data = rand(10,2); plot_data(:,1) = plot_data(:,1) - 0.5;
plot(plot_data(:,1),plot_data(:,2),'o');
xlabel('mu'); ylabel('sig')

title('Reference Window')

 x = -1:0.01:1;

% First window
for i = 1:size(plot_data,1)
    interaction_data1{i} = [x(:),reshape(normpdf(x,plot_data(i,1),plot_data(i,2)),[],1)];
end

add_plot_interaction(gca,interaction_data1)

xlabel('Interaction Term X')
ylabel('Interaction Term Y')

title('Interaction Window 1')

% Second window
[xm,ym] = meshgrid(x);

for i = 1:size(plot_data,1)
    interaction_data2{i} = normpdf(xm,plot_data(i,1),plot_data(i,2)).*normpdf(ym,plot_data(i,1),plot_data(i,2));
end

add_plot_interaction(gca,interaction_data2)

xlabel('Interaction Term X')
ylabel('Interaction Term Y')

title('Interaction window 2')

suptitle('Interaction figure with multiple viewing windows')
