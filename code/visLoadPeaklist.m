function [ppm_list] = visLoadPeaklist(XALN,loadings,ppm,corr_cutoff,show_on_original)
%{
    Chris Esselman 4.19.23
    Edited log
        Chris Esselman 5.18.23- added functionality for full resolution
        spectra
        Chirs Esselman 5.31.23- Automatically detect if peak picked or full
        resolution
        Chris Esselman 9.20.24- Combined Showpeaks and autopick (explained
        below) into one function- visLoadPeaklist
    VisLoad_autopick
    
        An added functionality to VisLoadings1D where ppm values of high
        correlation peaks are returned as an array.

    VisLoad_Showpeaks
    
        A functionality to VisLoadings1D where the peaks above the absolute
        value of the correlation cutoff are circled on an average of the
        normalized spectra instead of covariance figure. Meant to actually
        show where the peaks are on the orginal spectra. This function is a
        mix of VisLoad_auto and Displaypeak1D

 
    Inputs:
 
        XALN            Data matrix

        loadings        PCA or PLSDA loadings for one PC (ex/PCA.loadings(1,:)) 

        ppm             ppm values corresponding to XALN

    Optional Inputs:
    
        corr_cutoff     Value of desired correlation cutoff (ex/0.8). If
                        not specified, default 0.8

        show_on_original True or false. Will default to false
        

    Outputs:

        ppm_list        Array of ppm values corresponding to peaks with 
                        correlation above absolute value of corr_cutoff

        figure of visual loadings

Future work would be to have some way of telling if the ppm value was
positively or negatively correlated
%}

arguments
    XALN double
    loadings (1,:) double
    ppm (1,:) double
    corr_cutoff = 0.8
    show_on_original = false
end

if ~show_on_original
    %Sort ppm values from smallest to greatest
    %feature_sort = sorted ppm values, order = sorted indices of input
    [feature_sort,order]=sort(ppm);

    %Not actual correlation values(I wish it was)
    %Copies ordered loadings into variable corr
    corr=loadings(order);

    %Calculate covariance
    covar=loadings(order).*std(XALN(:,order));

    %create range from negative of largest loading to largest loading
    range=[-1*max(abs(corr)),max(abs(corr))];

    %Create colormap
    %Jet is just a specific type of map. Could also try turbo(100)?
    cmap=jet(100);

    %Loop works by creating a NaN matrix of 100 by number of ppm values;
    %k = Negative of largest loading to largest loading with 100 equal
    %intervals;
    %For each row of the NaN matrix, first find loadings larger than k, second
    %copy the covariance values of the largers loadings indices onto the row.

    %Will see later that if number is copied it will be colored with a new
    %color. If it remains NaN, it will not be recolored, so if every row of
    %matrix column is copied, the peak at that ppm value will be dark red.
    %If every row of matrix column is NaN, peak at that ppm value will be dark
    %blue.
    lines=NaN(size(cmap,1),size(corr,2));
    ind=1;
    for k=range(1):(range(2)-range(1))/size(cmap,1):range(2)
        lines(ind,find(corr>k))=covar(find(corr>k));
        ind=ind+1;
    end

    %Create the ppm list
    %First loop finds the index in the lines matrix where the correlation
    %coloring will be above or equal to the corr_cutoff
    peak_list_index = 1;
    select_correlation = (-1:2/100:1)';
    correlation_index = 0;
    for i = 51:size(select_correlation,1)
        if select_correlation(i) >= corr_cutoff
            correlation_index = i;
            break
        end
    end

    %Loops through columns of lines matrix.
    %First if statement, finds all ppm values where the NaN is located. If high
    %corr_cutoff, means that these values are highly negatively correlated to
    %creating the seperation on the PC. (Dark blue)

    %Second if statemnt, finds all ppm values that are not NaN. Means ppm
    %values are highly correlated to causing the seperation. (Dark red)

    %Included overarching if statement that will distinguish points from local
    %maxima for full resolution spectra. Chris Esselman 5.18.23
    peak_pick_cutoff = 5000;
    if length(ppm) < peak_pick_cutoff
        ppm_list = [];
        for i = 1:size(lines,2)
            dummy_variable = lines(((101-correlation_index) + 1),i);
            if isnan(dummy_variable) == 1
                ppm_list(peak_list_index) = feature_sort(i);
                peak_list_index = peak_list_index + 1;
            end
            dummy_variable = lines(correlation_index,i);
            if isnan(dummy_variable) == 0
                ppm_list(peak_list_index) = feature_sort(i);
                peak_list_index = peak_list_index + 1;
            end
        end
    elseif length(ppm) > peak_pick_cutoff
        ppm_list = [];
        X_ave = mean(XALN,1);
        for i = 2:(size(lines,2)-1)
            dummy_variable = lines(((101-correlation_index) + 1),i);
            if (isnan(dummy_variable) == 1) && (X_ave(1,i)>X_ave(1,i-1))...
                    && (X_ave(1,i)>X_ave(i+1))
                ppm_list(peak_list_index) = feature_sort(i);
                peak_list_index = peak_list_index + 1;
            end
            dummy_variable = lines(correlation_index,i);
            if (isnan(dummy_variable) == 0) && (X_ave(1,i)>X_ave(1,i-1))...
                    && (X_ave(1,i)>X_ave(i+1))
                ppm_list(peak_list_index) = feature_sort(i);
                peak_list_index = peak_list_index + 1;
            end
        end
    end

    %As mentioned above, plotting works by first making a plot that is a
    %uniform dark blue. Will continually redraw over the plot and peaks at
    %certain ppm values will change color based on if that value is NaN or not
    %in the lines matrix

    figure
    %subplot(3,2,6);
    plot(feature_sort,lines(1,:),'Color',cmap(1,:))
    hold on
    for k=2:size(cmap,1)
        %subplot(3,2,6);
        plot(feature_sort,lines(k,:),'Color',cmap(k,:));
    end
    xlabel('Chemical Shift')
    ylabel(['Loadings coefficients * std(data)'])

    colormap(cmap);
    t=colorbar;
    set(get(t,'ylabel'),'String', 'Loadings Coefficients');
    caxis([-1 1])
    set(gca,'XDir','rev');
    hold off
else
    %Sort ppm values from smallest to greatest
    %feature_sort = sorted ppm values, order = sorted indices of input
    [feature_sort,order]=sort(ppm);

    %Not actual correlation values(I wish it was)
    %Copies ordered loadings into variable corr
    corr=loadings(order);

    %Calculate covariance
    covar=loadings(order).*std(XALN(:,order));

    %create range from negative of largest loading to largest loading
    range=[-1*max(abs(corr)),max(abs(corr))];

    %Create colormap
    %Jet is just a specific type of map. Could also try turbo(100)?
    cmap=jet(100);

    %Loop works by creating a NaN matrix of 100 by number of ppm values;
    %k = Negative of largest loading to largest loading with 100 equal
    %intervals;
    %For each row of the NaN matrix, first find loadings larger than k, second
    %copy the covariance values of the largers loadings indices onto the row.

    %Will see later that if number is copied it will be colored with a new
    %color. If it remains NaN, it will not be recolored, so if every row of
    %matrix column is copied, the peak at that ppm value will be dark red.
    %If every row of matrix column is NaN, peak at that ppm value will be dark
    %blue.
    lines=NaN(size(cmap,1),size(corr,2));
    ind=1;
    for k=range(1):(range(2)-range(1))/size(cmap,1):range(2)
        lines(ind,find(corr>k))=covar(find(corr>k));
        ind=ind+1;
    end

    %Create the ppm list
    %First loop finds the index in the lines matrix where the correlation
    %coloring will be above or equal to the corr_cutoff
    peak_list_index = 1;
    select_correlation = (-1:2/100:1)';
    correlation_index = 0;
    for i = 51:size(select_correlation,1)
        if select_correlation(i) >= corr_cutoff
            correlation_index = i;
            break
        end
    end

    %Loops through columns of lines matrix.
    %First if statement, finds all ppm values where the NaN is located. If high
    %corr_cutoff, means that these values are highly negatively correlated to
    %creating the seperation on the PC. (Dark blue)

    %Second if statemnt, finds all ppm values that are not NaN. Means ppm
    %values are highly correlated to causing the seperation. (Dark red)

    %Included overarching if statement that will distinguish points from local
    %maxima for full resolution spectra. Chris Esselman 5.18.23

    %Set cutoff to know if full resolution or peak picked spectra
    peak_pick_cutoff = 5000;

    %List to know how to color the peak circles
    red_list = [];
    red_list_index = 1;
    blue_list = [];
    blue_list_index = 1;

    if length(ppm) < peak_pick_cutoff
        ppm_list = [];
        for i = 1:size(lines,2)
            dummy_variable = lines(((101-correlation_index) + 1),i);
            if isnan(dummy_variable) == 1
                ppm_list(peak_list_index) = feature_sort(i);
                blue_list(blue_list_index) = feature_sort(i);
                peak_list_index = peak_list_index + 1;
                blue_list_index = blue_list_index + 1;
            end
            dummy_variable = lines(correlation_index,i);
            if isnan(dummy_variable) == 0
                ppm_list(peak_list_index) = feature_sort(i);
                red_list(red_list_index) = feature_sort(i);
                peak_list_index = peak_list_index + 1;
                red_list_index = red_list_index + 1;
            end
        end
    elseif length(ppm) > peak_pick_cutoff
        ppm_list = [];
        X_ave = mean(XALN,1);
        for i = 2:(size(lines,2)-1)
            dummy_variable = lines(((101-correlation_index) + 1),i);
            if (isnan(dummy_variable) == 1) && (X_ave(1,i)>X_ave(1,i-1))...
                    && (X_ave(1,i)>X_ave(i+1))
                ppm_list(peak_list_index) = feature_sort(i);
                blue_list(blue_list_index) = feature_sort(i);
                peak_list_index = peak_list_index + 1;
                blue_list_index = blue_list_index + 1;
            end
            dummy_variable = lines(correlation_index,i);
            if (isnan(dummy_variable) == 0) && (X_ave(1,i)>X_ave(1,i-1))...
                    && (X_ave(1,i)>X_ave(i+1))
                ppm_list(peak_list_index) = feature_sort(i);
                red_list(red_list_index) = feature_sort(i);
                peak_list_index = peak_list_index + 1;
                red_list_index = red_list_index + 1;
            end
        end
    end

    %Create the figure that shows peaks above absolute value of correlation
    %cutoff

    %Create average spectra from normalized spectra
    X_ave = mean(XALN,1);

    %Plot the ppm on the X and the intensities on y
    figure, hold on;
    plot(ppm,X_ave)

    %Create circles around important peaks or ppms. Blue circles for high
    %negative correlation and red circles for high positive correlation
    if blue_list_index ~= 1
        for k=1:length(blue_list)
            [~,b(k)]=min(abs(ppm-blue_list(k)));
        end
        scatter(blue_list,X_ave(:,b),'blue');
    end

    if red_list_index ~= 1
        for k=1:length(red_list)
            [~,r(k)]=min(abs(ppm-red_list(k)));
        end
        scatter(red_list,X_ave(:,r),'red');
    end
    set(gca,'XDir','reverse');
end
end