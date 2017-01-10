%% Add interactivity to any plot
%
% Talfan Evans (November 2016)
%
% For any existing plot, this function adds a subplot window which displays
% data corresponding to the points in the original plot. The rightmost plot
% (which we term 'interaction data') updates to display data corresponding
% to clickable datapoints on the left (the original data).
%
% Example:
%
% Visualise mu*sigma of a 2d gaussian with sigma values varying on the x and mu on the y axes
% close all
% clearvars
%
% murange = 0:0.02:1;
% sigrange = 0:0.01:1;
%
% [sig,mu] = meshgrid(sigrange,murange);
% ellipt = mu.*sig;
%
% imagesc(ellipt)
% xlabel('Sig. 1')
% ylabel('Sig. 2')
% title('Ellipticity of gaussian')
%
% xrange = 0:0.01:1;
% interaction_data = cell(length(sigrange),length(murange));
% for i = 1:length(sigrange)
%     for j = 1:length(murange)
%         interaction_data{i,j} = [xrange;normpdf(xrange,murange(j),sigrange(i))];
%     end
% end
%
% add_plot_interaction_new(gca,interaction_data,'o-')
% xlabel('X')
% ylabel('Y')

function add_plot_interaction(axis_handle,interaction_data,varargin)

% Check whether the existing plot already contains an interaction
if any(strcmpi(get(allchild(axis_handle),'Tag'),'current_point'))
    delete(axis_handle.Children(strcmpi(get(allchild(axis_handle),'Tag'),'current_point')))
end

%% Default values
opts = [];

%% Check that the axis_handle provided is compatible
[axis_handle,opts] = parse_axis_handle(axis_handle,opts);

%% Get the data stored in the existing axis_handle
opts = get_existing_data(axis_handle,opts);

%% Check that the interaction data is compatible
[interaction_data,opts,axis_handle] = parse_interaction_data(interaction_data,opts,axis_handle);
 
%% Copy existing figure into a new subplot figure
[hax,opts] = copy_to_new_canvas(axis_handle,interaction_data,opts,varargin{:});

%% Assign the callback function to the new axis
function_cell = cell(1,opts.Nsubplts-1);
for i = 1:opts.Nsubplts-1
    function_cell{i} = eval(['@(src,event)fig_callback(src,event,hax,',num2str(i+1),')']);
end

set(hax(1),'ButtonDownFcn',...
    @(src,event)(cellfun(@(x)feval(x,src,event),function_cell)))

%% Delete original figure
delete(axis_handle.Parent)                                                 % Delete the original figure

end

%% Make sure the axis_handle is in the required format and determine its
% format, if the format is not explicitly supplied
function [axis_handle,opts] = parse_axis_handle(axis_handle,opts)

% If the handle provided is to the figure, rather than axis of an existing
% interaction figure, find the axis handle ofthe original axis
if strcmpi(axis_handle.Tag,'interactive_figure')
    axis_handle = axis_handle.Children(...
        strcmpi(get(axis_handle.Children,'Tag'),...        % Find the indices of the existing interaction subplots
            'interactive_figure_base'));    
% If the handle provided is a subplot of an already existing subfigure, get
% the original axis_handle
elseif strcmpi(axis_handle.Parent.Tag,'interactive_figure')
    axis_handle = axis_handle.Parent.Children(...
        strcmpi(get(axis_handle.Parent.Children,'Tag'),...        % Find the indices of the existing interaction subplots
            'interactive_figure_base'));
% Determine whether the supplied axis contains coordinate data (lines,
% points etc.), matrix data (e.g. an image of pixels) or bars (e.g. a
% histogram)
% If a figure handle has been supplied and that figure only has one child
% axis, use this axis as the axis_handle. Otherwise, return error.
elseif strcmpi(class(axis_handle),'matlab.graphics.axis.Axes')
    % Everything is in order...
elseif strcmpi(class(axis_handle),'matlab.ui.Figure')
    if length(axis_handle.Children)==1
        if strcmpi(axis_handle.Children(1),'matlab.graphics.axis.Axes')
            axis_handle = axis_handle.Children(1);
        else
            error('Need to supply an axis handle, not a figure handle.')
        end
    else
        error('Need to supply an axis handle, not a figure handle.')
    end
% Sometimes gca grabs the objects of the axis, not the axis. If this is the
% case, take the parent as the axis_handle.
elseif strcmpi(class(axis_handle.Parent),'matlab.graphics.axis.Axes')
    axis_handle = axis_handle.Parent;
% Otherwise throw an error
else
    error('Need to supply an axis handle')
end

% Determine the type data being displayed. 'Points' data is data 
% supplied as two x and y vectors. 'matrix' data is data supplied as a 
% matrix, e.g. pixels in an image. 'bar' data is in x and y vector format,
% but is dealt with differently later on.
if     strcmpi(class(axis_handle.Children),'matlab.graphics.chart.primitive.Line')
    axis_handle.UserData.plot_type = 'points';
elseif strcmpi(class(axis_handle.Children),'matlab.graphics.chart.Scatter')
    axis_handle.UserData.plot_type = 'points';
elseif strcmpi(class(axis_handle.Children),'matlab.graphics.chart.primitive.Bar')
    axis_handle.UserData.plot_type = 'bar';
elseif strcmpi(class(axis_handle.Children),'matlab.graphics.primitive.Image')
    axis_handle.UserData.plot_type = 'image';
end

% Check whether the existing plot already contains an interaction
if any(strcmpi(get(allchild(axis_handle),'Tag'),'current_point'))
    delete(axis_handle.Children(strcmpi(get(allchild(axis_handle),'Tag'),'current_point')))
end

end

%% Get the data stored in the existing axis_handle
function opts = get_existing_data(axis_handle,opts)
% Find the number of distinct datasets represented on the figure. Plots can
% have many lines overlaid, each with its own set of points. Data will be
% stored in a [1,Ndata] cell array.
opts.Ndata = length(axis_handle.Children);

    switch lower(axis_handle.UserData.plot_type)
        case {'points','bar'}
            % The data for points and bar types are stored in thhe same way
            for n = 1:opts.Ndata
                axis_handle.UserData.plot_data{n} = [axis_handle.Children(n).XData(:),axis_handle.Children(n).YData(:)];
            end
        case 'image'
            % Treat image data differently. Since it doesn't contain
            % explicit points, we will assign coordinates to each of the
            % pixels. Need to generate the x and y coordinate vectors as
            % they are not stored in image handles unless they were
            % explicitly specified when displaying the image
            xvec = linspace(axis_handle.Children(1).XData(1),axis_handle.Children(1).XData(end),size(axis_handle.Children(1).CData,2));
            yvec = linspace(axis_handle.Children(1).YData(1),axis_handle.Children(1).YData(end),size(axis_handle.Children(1).CData,1));
            axis_handle.UserData.plot_data{1}.X = xvec(:);
            axis_handle.UserData.plot_data{1}.Y = yvec(:);
            axis_handle.UserData.plot_data{1}.Z = axis_handle.Children(1).CData;
    end
    
end

%% Make sure the interaction_data is in the required format and determine its
% format, if the format is not explicitly supplied
function [interaction_data,opts,axis_handle] = parse_interaction_data(interaction_data,opts,axis_handle)
switch lower(axis_handle.UserData.plot_type)
    case {'points','bar'}
        % 'points' and 'bar' type figures can contain mutliple
        % datasets. The interaction data needs to provide a
        % {Ndata}{Npoints} cell array, where Npoints is the number of
        % points corresponding to each dataset.
        
        % Firstly, if there is only one dataset, the user may provide
        % the interaction_data as a single layer cell array. If this is
        % the case, put it in the right format.
        if (opts.Ndata==1) && (length(interaction_data)>1)
            if (length(interaction_data)==size(axis_handle.UserData.plot_data{1},1))
                interaction_data = {interaction_data};
            else
                error('The interaction data must have the same dimension as the plot data')
            end
            % If there is more than one data set in the plot, check that
            % the interaction data has the right dimensions
        elseif (opts.Ndata>1)
            if length(interaction_data)==opts.Ndata
                for n = 1:opts.Ndata
                    if length(interaction_data{n})~=size(axis_handle.UserData.plot_data{n},1)
                        error('Number of data points in data set %i does not match number of points in the existing plot. Interaction data must be provided as cell array of dimensions {NData}{NPoints}, where NData is the number of distinct datasets in the existing plot and NPoints the number of datapoints in each dataset',n)
                    end
                end
            else
                error('Interaction data must be provided as cell array of dimensions {NData}{NPoints}, where NData is the number of distinct datasets in the existing plot and NPoints the number of datapoints in each dataset')
            end
        end
        
        % Store information about the subplots in the axis handle
        if ~isfield(axis_handle.UserData,'subplot_type')
            axis_handle.UserData.subplot_type = {};
        end
        
        % If interaction data has smallest dimension greater than 2, assume its an
        % image        
        if min(size(interaction_data{1}{1}))>2
            axis_handle.UserData.subplot_type{end+1} = 'image';
            % If interaction data has one dimension equal to 2, assume it's x and y data
        elseif min(size(interaction_data{1}{1}))==2
            axis_handle.UserData.subplot_type{end+1} = 'points';
        else
            error('Interaction data must either be supplied as a [2,L] vector of coordinates or an image array')
        end
    case 'image'
        % Image interaction data needs to be provided in a cell array
        % that is the same dimension as the image (one cell for each
        % pixel)
        if     all(size(interaction_data)==[length(axis_handle.UserData.plot_data{1}.Y),length(axis_handle.UserData.plot_data{1}.X)])~=1
            if all(size(interaction_data)==[length(axis_handle.UserData.plot_data{1}.X),length(axis_handle.UserData.plot_data{1}.Y)])==1
                interaction_data = interaction_data';
            else
                error('Interaction data for image object must be provided as a {N_Y,N_X} cell array')
            end
        end
        
        % Store information about the subplots in the axis handle
        if ~isfield(axis_handle.UserData,'subplot_type')
            axis_handle.UserData.subplot_type = {};
        end
        
        % If interaction data has smallest dimension greater than 2, assume its an
        % image
        if min(size(interaction_data{1,1}))>2
            axis_handle.UserData.subplot_type{end+1} = 'image';
            % If interaction data has one dimension equal to 2, assume it's x and y data
        elseif any(min(size(interaction_data{1,1}))==[1,2])
            axis_handle.UserData.subplot_type{end+1} = 'points';
        else
            error('Interaction data must either be supplied as a [2,L] vector of coordinates or an image array')
        end
        
        % Reshape into a nested cell array so as to be in the same format
        % as the points data, which saves on writing flexible code later on
        in_data = cell(1,size(interaction_data,1));
        for k=1:size(interaction_data,1)
            in_data{k} = interaction_data(k,:);
        end
        interaction_data = in_data; clear in_data;
        
end

end

%% Copy the existing figure object to a new canvas
function [hax,opts] = copy_to_new_canvas(axis_handle,interaction_data,opts,varargin)

if strcmpi(axis_handle.Parent.Tag,'interactive_figure')
   N = length(axis_handle.Parent.Children) + 1; 
else
    N = 2;
end
opts.Nsubplts = N;
M = 1; p =1;

hf1=figure;                                                                % Create a new canvas
hf1.Tag = 'interactive_figure';                                            % Allow the figure to be identified as an interactive figure later on
s1=subplot(M,N,p); p=p+1;                                                  % Place an empty subplot where the existing plot will be copied to
pos1=get(s1,'Position');                                                   % We don't need the actual subplot canvas, we should use Matlab's default placement to help position our plot
delete(s1);                                                                % We have the position of Matlab's default placement, now we can get rid of the actual canvas

hax(1)=copyobj(axis_handle,hf1);                                           % Copy the existing figure to the new canvas
set(hax(1), 'Position', pos1);                                             % Place it where the old subplot was
set(hax(1),'Tag','interactive_figure_base')                                % Tag this sbplot as th original figure

hold on; plot(0,0,'ok','MarkerSize',10);                                   % Add a marker to indicate the current datapoint
set(hax(1).Children(1),'Visible','off')                                    % Make the initial points invisible until the user makes the first click
set(hax(1).Children(1),'Tag','current_point')                              % Need to refer to this later
set(allchild(hax(1)),'HitTest','off')                                      % Set the data lines so as not to get in the way of mouse-clicks

child_ind = find(strcmpi(get(axis_handle.Parent.Children,'Tag'),...        % Find the indices of the existing interaction subplots
    'interaction_data'))';
child_ind = fliplr(child_ind);                                             % Matlab orders the axes in reverse
for pl = 1:length(child_ind)                                               % Copy over pre-existing interaction plots
    s1=subplot(M,N,p); p=p+1;                                              % Place an empty subplot where the existing plot will be copied to
    pos1=get(s1,'Position');                                               % We don't need the actual subplot canvas, we should use Matlab's default placement to help position our plot
    delete(s1);                                                            % We have the position of Matlab's default placement, now we can get rid of the actual canvas
    hax(pl+1)=copyobj(axis_handle.Parent.Children(child_ind(pl)),hf1);     % Copy the existing figure to the new canvas
    set(hax(pl+1),'Position', pos1);                                       % Place it where the old subplot was
    set(allchild(hax(pl+1)),'Visible','off')                               % Make th figure initially invisible
end

hax(N)=subplot(M,N,p); p=p+1; axis equal                                   % Create second subplot to hold the interaction data

switch lower(axis_handle.UserData.subplot_type{N-1})
    case {'points','bar'}
        if min(size(interaction_data{1}{1}))==1                            % If x vector is not explicitly included
            plot(interaction_data{1}{1},varargin{:})
        else
            plot(interaction_data{1}{1}(:,1),...
                interaction_data{1}{1}(:,2),varargin{:})                   % Plot some example data but make it invisible before any specific point is clicked
        end
    case 'image'
        imagesc(interaction_data{1}{1},varargin{:}); axis equal            % There can be only one datasets for image data
        xlim([1,size(interaction_data{1}{1},2)])
        ylim([1,size(interaction_data{1}{1},1)])
end

set(allchild(hax(N)),'Visible','off')                                      % We can't just leave a blank figure as we use the set(gca,'XData',xdata) method later
set(hax(N),'Tag','interaction_data')                                       % Identify as a subplot displaying interaction data
if isfield(hax(1).UserData,'interaction_data')
    hax(1).UserData.interaction_data{end+1} = interaction_data;            % Store the interaction data in the figure handle
else
    hax(1).UserData.interaction_data{1} = interaction_data;
end

end

%% The callback function which updates the interaction data according to the
% current data-point on the original axis
function fig_callback(~,~,hax,hax_number) %#ok<DEFNU>

HAX = hax(hax_number);                                                     % Set the correct subplot

interaction_data = hax(1).UserData.interaction_data{hax_number-1};         % Get the interaction data, which is stored in the axis handle's UserData field

set(allchild(HAX),'Visible','on')                                          % Set the lines to be visible since we turned them off earlier

c = get(gca,'CurrentPoint'); c=[c(1,1),c(1,2)];                            % Get the position of the mouse-click

% Get the indices of the corresponding cell
switch lower(hax(1).UserData.plot_type)
    case {'points'}
        for i = 1:length(hax(1).UserData.plot_data)
            [mincell(i),Ind2(i)] = min(sqrt( (hax(1).UserData.plot_data{i}(:,1)-c(1)).^2 + ...% Calculate the closest corresponding data-point
                (hax(1).UserData.plot_data{i}(:,2)-c(2)).^2));   %#ok<AGROW>
        end
        
        [~,Ind1] = min(mincell);                                         % Find index of minimum point from all datasets within the vector containing the indicies of minimum points from all datasets
        Ind2 = Ind2(Ind1);                                                 % Find index of minimum points within the datasets containing the minimum point
        
        mrkr = hax(1).UserData.plot_data{Ind1}(Ind2,:);                             % Find coordinates of closest datapoint
    case 'image'        
        Ind1 = limRound(c(2),1,length(interaction_data)); % y
        Ind2 = limRound(c(1),1,length(interaction_data{1})); % x
        
        mrkr = [Ind2,Ind1];
end

% Plot the cell
switch lower(hax(1).UserData.subplot_type{hax_number-1})
            case 'points'
                ch = get(gca,'Children');                                          % Set a marker to indicate that the
                set(ch(1),'XData',mrkr(1),'YData',mrkr(2),'Visible','on')          % datapoint is currently selected
                                
                if size(interaction_data{Ind1}{Ind2},2)>2                         % Transpose the interaction data if necessary
                    interaction_data{Ind1}{Ind2} = interaction_data{Ind1}{Ind2}';
                end
                
                if min(size(interaction_data{Ind1}{Ind2}))==1
                    HAX.Children.XData = 1:length(interaction_data{Ind1}{Ind2}); % Generate an xvector if it was not explicitly defined
                    HAX.Children.YData = interaction_data{Ind1}{Ind2};
                else
                    HAX.Children.XData = interaction_data{Ind1}{Ind2}(:,1);          % Update the plot on the right to correspond
                    HAX.Children.YData = interaction_data{Ind1}{Ind2}(:,2);          % to the relevant data point
                end
            case 'bar'
                % Need to add this functionality
            case 'image'
                ch = get(gca,'Children');                                          % Set a marker to indicate that the
                set(ch(1),'XData',mrkr(1),'YData',mrkr(2),'Visible','on')                % datapoint is currently selected
                
                HAX.Children.CData = interaction_data{Ind1}{Ind2};                % Update the plot on the right to correspond to the relevant point
end
 
end

%% Round number within limits
function x = limRound(x,LB,UB)

x = round(x);

x(x<LB) = LB;
x(x>UB) = UB;

end