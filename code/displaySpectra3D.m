function displaySpectra3D(X,ppm,isbackground,database_gissmo,annotate_matrix,delete_annotation,delete_width)

%{
    Chris Esselman 9.21.24
    
    displaySpectra3D
        Create a 3D stacked plot. Initially for visualizing fraction
        library, but can be used for visualizing more than one spectrum.
        Zooming in easiest when clipping is turned off. Can also use the
        xlim function after this function to get a more specific zoom. 
        
        displaySpectra3D(X,ppm)
        xlim([0.5 2])

    Inputs:
    
    X               Data matrix ex/ 140x32k double

    ppm             ppm values corresponding to data matrix ex/ 1x32k

    Optional Inputs:

    isbackground        True or false. Will color the last spectrum
                        black

    database_gissmo     struct containing ppm ranges of peaks in GISSMO
    
    annotate_matrix     matrix containing three columns of strings
                        1. Metabolite name
                        2. Lower fraction number
                        3. Higher fraction number

    delete_annotation   boolean. True if want to delete instead of color

    delete_width        double. Change the width of points to be deleted
                        ex/50
        

%}

% So I can pass optional arguments
arguments
    X double
    ppm (1,:) double
    isbackground logical = false
    database_gissmo = 53
    annotate_matrix = 54
    delete_annotation = false
    delete_width = 0
end

% Change the range of values to visualize
X_shrunk = X;
ppm_shrunk = ppm;

% Make so values between 0 and 1
X_shrunk=(X_shrunk-min(X_shrunk(:)))/(max(X_shrunk(:))-min(X_shrunk(:)));
% Find out how many fracs there are
rows = size(X_shrunk,1);
% Make matrix for plot3 function. Matrix with rows number of fractions and
% columns the size of ppm. Values are repeated ppm vector
for i = 1:rows
    x(i,:) = ppm_shrunk;
end
%Values that hold the height of the matrix
y = X_shrunk;
%So the fracs spread out
for i = 1:rows
    z(i,:) = (zeros(1,length(ppm_shrunk)) + i);
end
%Color the plot. Could change this to be a different color map
colors = winter(rows);
fig = figure();
ax1 = axes(fig);
hold on
%Plot each fraction
if ~delete_annotation && ~isbackground
    for i = 1:rows - 1
        plot3(x(i,:),z(i,:),y(i,:),'Color',colors(i,:),'LineWidth',0.5)
    end
    plot3(x(rows,:),z(rows,:),y(rows,:),'Color',colors(i,:),'LineWidth',0.5)
elseif ~delete_annotation && isbackground
    for i = 1:rows - 1
        plot3(x(i,:),z(i,:),y(i,:),'Color',colors(i,:),'LineWidth',0.5)
    end
    plot3(x(rows,:),z(rows,:),y(rows,:),'Color',"k",'LineWidth',0.5)
end
% Do the annotation
% First see if input arguments are correct
if isstruct(database_gissmo)
    if ~isstring(annotate_matrix)
        error('Error. To annotate need both database_gissmo and annotate_matrix')
    else
        struct_annotate2 = struct;
        for i = 1:size(annotate_matrix,1)
            struct_annotate2(i).metab = annotate_matrix(i,1);
            struct_annotate2(i).lower_frac = str2double(annotate_matrix(i,2));
            struct_annotate2(i).high_frac = str2double(annotate_matrix(i,3));
        end
        if delete_annotation
            for k = 1:size(struct_annotate2,2)
                for i = 1:size(database_gissmo,2)
                    if strcmp(struct_annotate2(k).metab,database_gissmo(i).Metab) == 1
                        index_metab = i;
                        break;
                    end
                end
                low_frac = struct_annotate2(k).lower_frac;
                high_frac = struct_annotate2(k).high_frac;
                for i = 1:2:length(database_gissmo(index_metab).ppm_range)
                    idx1 = matchPPMs(database_gissmo(index_metab).ppm_range(i),ppm_shrunk);
                    idx2 = matchPPMs(database_gissmo(index_metab).ppm_range(i+1),ppm_shrunk);
                    if (idx1 - delete_width) > 0
                        idx1 = idx1 - delete_width;
                    end
                    if (idx2 ~= 1) && ((idx2 + delete_width)  < length(ppm_shrunk))
                        idx2 = idx2 + delete_width;
                    end
                    X_shrunk(low_frac:high_frac,idx1:idx2) = NaN;
                end
            end
            y = X_shrunk;
            for i = 1:rows - 1
                plot3(x(i,:),z(i,:),y(i,:),'Color',colors(i,:),'LineWidth',0.5)
            end
            plot3(x(rows,:),z(rows,:),y(rows,:),'Color',"k",'LineWidth',0.5)
        else
            %Make a nxn NaN matrix and change all the columns and rows to
            %values. These will be colored red
            annotate_matrix_X = NaN(size(X_shrunk,1),size(X_shrunk,2));
            for k = 1:size(struct_annotate2,2)
                for i = 1:size(database_gissmo,2)
                    if strcmp(struct_annotate2(k).metab,database_gissmo(i).Metab) == 1
                        index_metab = i;
                        break;
                    end
                end
                low_frac = struct_annotate2(k).lower_frac;
                high_frac = struct_annotate2(k).high_frac;
                for i = 1:2:length(database_gissmo(index_metab).ppm_range)
                    idx1 = matchPPMs(database_gissmo(index_metab).ppm_range(i),ppm_shrunk);
                    idx2 = matchPPMs(database_gissmo(index_metab).ppm_range(i+1),ppm_shrunk);
                    annotate_matrix_X(low_frac:high_frac,idx1:idx2) = X_shrunk(low_frac:high_frac,idx1:idx2);
                end
            end
            % Plot on the red traces
            y = annotate_matrix_X;
            for i = 1:size(annotate_matrix_X,1)
                plot3(x(i,:),z(i,:),y(i,:),'Color',"r",'LineWidth',0.5)
            end
        end
    end
end
hold off
% Make struct to hold to color over
set(gca,'XDir','rev');
view(0,45)
originalZLim = zlim(ax1);
ax1.Clipping = "off";
clipping_variable = false;
% Create a slider control
slider = uicontrol('Style', 'slider', 'Min', 0.01, 'Max', 50, 'Value', 1, ...
    'Position', [100, 50, 220, 20], ...
    'Callback', @(src, event) adjustVerticalZoom(ax1, src, originalZLim));
% slider = uicontrol('Style', 'slider', 'Min', 0.01, 'Max', 100, 'Value', 1, ...
%                    'Position', [100, 50, 220, 20], ...
%                    'Callback', @(src, event) adjustVerticalZoom(ax, src, originalZLim));

% Create a text label for the slider
sliderLabel = uicontrol('Style', 'text', 'Position', [100, 70, 120, 20], ...
    'String', 'Vertical Zoom');

% Create a button to enable the selection mode
selectionButton = uicontrol('Style', 'pushbutton', 'String', 'Change Clipping', ...
    'Position', [100, 100, 120, 20], ...
    'Callback', @(src, event) selectRegion(ax1));

%% function defs
    function adjustVerticalZoom(ax, slider, ~)
        % Get the current slider value
        zoomFactor = slider.Value;

        % Adjust the upper limit of the z-axis based on the zoom factor
        zData = get(ax.Children, 'ZData');
        if iscell(zData)
            zData = cell2mat(zData);
        end
        zMax = max(zData(:));
        baseline = -zMax / zoomFactor/3; % Ensure the baseline is always visible
        zlim(ax, [baseline, zMax / zoomFactor]);
    end
    function selectRegion(ax)
        if clipping_variable
            ax.Clipping = "off";
            clipping_variable = false;
        else
            ax.Clipping = "on";
            clipping_variable = true;
        end
    end
end
