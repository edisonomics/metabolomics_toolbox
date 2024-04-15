%% About this workflow:
% ----------------------------------------------------------------------------
% Workflow type: NMR metabolomics-toolbox tutorial
% Workflow ID: matlab_workflow1_complete.nan.m
% Experiment type:1D-1H
% Steps coverd: From loading data to PCA
% Time estimate to complete: 2 days
% (For a summarized version that needs less computational time, use:
% nan_matlab_workflow1_summarized.m)
% ----------------------------------------------------------------------------

%% Acknowledgements:
% ----------------------------------------------------------------------------
% This workflow was originally developed by Olatomiwa O Bifarin at UGA.
% For further details, please see Bifarin et al. (J. Proteome Res., 2021).
% The original workflow was modified for the purpose of this tutorial.
% ----------------------------------------------------------------------------

%% Research Question: What are the NMR features separating the groups?

%% Script's architecture
% ------------------------------------------------------------------------
%  SECTION 1: ACTIVATE TOOLBOX, SET PATHS, AND LOAD ORIGINAL DATA
%  SECTION 2: REFERENCING
%  SECTION 3  REGION REMOVAL
%  SECTION 4: ALIGNMENT
%  SECTION 5: NORMALIZATION
%  SECTION 6: SCALING
%  SECTION 7: MULTIVARIATE ANALYSIS
%  SECTION 8: SUMMARY
% ------------------------------------------------------------------------

%% SECTION 1: ACTIVATE TOOLBOX, SET PATHS, AND LOAD ORIGINAL DATA

%% Activate metabolomics_toolbox in NMRbox
activate_metabolomics_toolbox
    % ----------------------------------------------------------------
    % This step is ONLY for NMRbox. Skip this step if you are using MATLAB on your laptop.
    % 
    % For more details about this toolbox, go to:
    % https://github.com/edisonomics/metabolomics_toolbox
    % ----------------------------------------------------------------

%% Set paths
     % ----------------------------------------------------------------
     % Note:
     % Here, we set paths to directories. You move between paths during the
     % analysis so that you can structure the workflow (for example, you
     % load your original data from 'paths.data', but you save your processed
     % data in 'paths.results').
     % ----------------------------------------------------------------
paths.project = findCurrentFile()
paths.scripts = [paths.project,'/','scripts']
paths.data = [paths.project,'/','data'] 
paths.spectra = [paths.data,'/','spectra'] 
paths.results = [paths.project,'/','results']

%% Load NMR data into MATLAB (function 'Load1D')
data.spectra=Load1D(paths.spectra,'bruker');
     % ----------------------------------------------------------------
     % Note:
     % Sample 269 (SS054) was removed from the dataset in advance
     % based on the decision made in the original study.
     % ----------------------------------------------------------------

%% Create a data matrix using the data loaded (function 'Setup1D')
[data.X,data.ppm]=Setup1D(data.spectra); 

%% Visualize the orignal spectra (function 'plotr')
 plotr(data.ppm,data.X); set(gca,'FontSize',18)

%% Import a study metadata sheet
cd(paths.data); data.T = readtable('sample_removedSS054.xlsx','Format','auto');

%% Create a vector 'Yvec'
     % ----------------------------------------------------------------
     % Note:
     % 'Yvec' (i.e. Y vector) will be used to select samples of interest in
     % the downstream analysis.
     % In this study, to select samples for the Control group, we use Yvec =
     % 1, but for the Patient group, we use Yvec = 2. 
     % ----------------------------------------------------------------
data.Yvec=data.T.Yvec;

%% SECTION 2: REFERENCING

%% Reference spectra (function 'ref_spectra')
    % ----------------------------------------------------------------
    % Note:
    % (1) After you excute 'ref_spectra', a figure showing picked peaks wil automatically
    % pop up. You need to click on one circle that you want to make 0 ppm.
    % (2) Next, a figure showing referenced spectra will pop up. Zoom in, and make sure
    % they are referenced correctly.
    % ----------------------------------------------------------------
data.spectra1 = ref_spectra(data.spectra,-0.02);

%% Create a new data matrix using the referenced spectra
[data.X1,data.ppm]=Setup1D(data.spectra1);

%% If you want to see the new matrix you just created
displaypeak1D(data.X1,data.ppm,0,data.Yvec);

%% SECTION 3  REGION REMOVAL

%% Remove end and water regions ('remove_region')
    % ----------------------------------------------------------------
    % Note:
    % 'remove_region' removes regions that are not necessary in the downstream
    % analysis.
    % Here, we remove end regions (<-0.5 ppm and >10 ppm), and the water
    % region (4.683-4.893 ppm) 
    % ----------------------------------------------------------------

data.ppmEndMin=-0.5       % Change the value here based on your spectra
data.ppmEndMax=10.0      % Change the value here based on your spectra
data.ppmWaterMin=4.683  % Change the value here based on your spectra
data.ppmWaterMax=4.893 % Change the value here based on your spectra

[~,data.k1]=min(abs(data.ppm-data.ppmEndMin));
data.XR=data.X1(:,data.k1:end);
data.ppmR=data.ppm(data.k1:end);
[~,data.k2]=min(abs(data.ppmR-data.ppmEndMax));
data.XR1=data.XR(:,1:data.k2);
data.ppmR=data.ppmR(1:data.k2);
data.XR2=remove_region(data.XR1,data.ppmR,data.ppmWaterMax,data.ppmWaterMin);
%% Make sure regions were removed properly
displaypeak1D(data.XR2,data.ppmR,0,data.Yvec); set(gca,'FontSize',18)

%% SECTION 4: ALIGNMENT

%% Align spectra (function 'guide_align1D')
    % ----------------------------------------------------------------
    % Note:
    % 'guide_align1D' aligns spectra using available distance_metric
    % (either 'correlation' or 'spearman') and method (either
    % 'CCOW','RAFFT', or 'PAFFT') parameters.
    % ('ICOSHIFT' is also one of the options but does not currently work.)
    %
    % Here, we create a set of aligned spectra using combinations of the parameters.
    %
    % CAUTION!!!: This section needs one overnight (mainly because of CCOW).
    % Execute this subsection, go home, and check the output next morning.
    % ----------------------------------------------------------------

    distance_metric={'correlation', 'spearman'};
    method={'RAFFT', 'PAFFT', 'CCOW'};
    a=1:length(distance_metric);
    b=1:length(method);
    comb=combvec(a,b)';
    for i=1:length(comb)
        data.XAL{i,1}=distance_metric(comb(i,1));
        data.XAL{i,2}=method(comb(i,2));
        tic
        data.XAL{i,3}=guide_align1D(data.XR2,data.ppmR,string(data.XAL{i,1}),string(data.XAL{i,2}));
        data.XAL{i,4}=toc
    end

    clear('distance_metric','method','a','b','comb','i')

%% Compare the aligned spectra
    % ----------------------------------------------------------------
    % Note:
    % Zoom in and decide which alignment you use.
    % (Funciton 'multiAlign_interactive_draft_2' is also abailable as an
    % alternative but we skip it here.)
    % ----------------------------------------------------------------

distance_metric={'spearman'}; % Enter the one to use: 'correlation' or 'spearman'
method={'CCOW'};                   % Enter the one to use: 'RAFFT', 'PAFFT', or 'CCOW'

for i= 1:length(data.XAL)
    a(i,1) = i;
    a(i,2)=strcmp(data.XAL{i,1},distance_metric);
    a(i,3)=strcmp(data.XAL{i,2},method);
end
row=a(a(:,2)==1&a(:,3)==1);
displaypeak1D(data.XAL{row, 3},data.ppmR,0,data.Yvec); set(gca,'FontSize',18); title(strcat(distance_metric,'+',method));

clear('a','i','distance_metric','method', 'row')

%% Decide aligned spectra to use
    % ----------------------------------------------------------------
    % Note:
    % The original study decided to use a combination of spearman and CCOW.
    % ----------------------------------------------------------------

distance_metric={'spearman'}; % Enter the one to use: 'correlation' or 'spearman'
method={'CCOW'};                   % Enter the one to use: 'RAFFT', 'PAFFT', or 'CCOW'

for i= 1:length(data.XAL)
    a(i,1) = i;
    a(i,2)=strcmp(data.XAL{i,1},distance_metric);
    a(i,3)=strcmp(data.XAL{i,2},method);
end

a(a(:,2)==1&a(:,3)==1);

data.distance_metric=string(distance_metric);
data.method=string(method);
data.XAL1=data.XAL{a(a(:,2)==1&a(:,3)==1),3};

clear('a','i','distance_metric','method')

%% If you want to see the finalized spectra
displaypeak1D(data.XAL1,data.ppmR,0,data.Yvec); set(gca,'FontSize',18)

%% SECTION 5: NORMALIZATION

%% Limit the dataset to the real samples
    % ----------------------------------------------------------------
    % Note:
    % Hereafter, we limit the samples to real study samples.
    % This can eliminate the effects of other samples (Internal/External
    % controls, blanks) on normalization.
    %
    % Here, we use only Control and Patient group samples (Yvec=1 and Yvec=2,
    % respectibvely).
    % ----------------------------------------------------------------

%% Overview the data disribution of the original data (function 'normcheck')
normcheck(data.XAL1(data.Yvec==1|data.Yvec==2,:))

%% Normalization (function 'normalize')
    % ----------------------------------------------------------------
    % Note:
    % Several normalization methods ('total', 'PQN', 'quantile',
    % 'intensity', and 'integral') are available. This study used PQN.
    % ----------------------------------------------------------------
data.normalization=['PQN']; % Enter the one to use: 'total', 'PQN', 'quantile', 'intensity', or 'integral
data.XAL1N=normalize(data.XAL1(data.Yvec==1|data.Yvec==2,:),data.ppmR,data.normalization);

%%  Overview the data disribution of the normalized data
normcheck(data.XAL1N)

%% SECTION 6: SCALING

%%  Create copies of the normalized data for the downstream analysis
analysis.XAL1N=data.XAL1N;
analysis.ppmR=data.ppmR;
analysis.Yvec=data.Yvec(data.Yvec==1|data.Yvec==2,:);

%% Scale data (function 'scale')
    % ----------------------------------------------------------------
    % Note:
    % Several scaling methods ('log', 'logoff', 'mc', and 'pareto') are available. 
    % This study used logoff.
    % ----------------------------------------------------------------
analysis.scaling=['logoff'];  % Enter the one to use: 'log', 'logoff', 'mc', and 'pareto'
analysis.XAL1NS=scale(analysis.XAL1N,analysis.scaling)

%% Check scaling for logoff
varcheck(analysis.XAL1NS)

%% SECTION 7: MULTIVARIATE ANALYSIS

%% Run Principal Component Analysis (function 'nipalsPCA')
analysis.PCA2logoff=nipalsPCA(analysis.XAL1NS,10);

%% Check the contribution of each component
plot(analysis.PCA2logoff.variance,'o-k','MarkerFaceColor','k','MarkerSize',6)
title('Contribution of each component');
xlabel('Component');
ylabel('%');

%% Create score plots (function 'VisScores')
VisScores(analysis.XAL1NS,analysis.PCA2logoff,[1 2],'Y',analysis.Yvec,'conf_ellipse',false, 'showlegend', analysis.Yvec);
    % ----------------------------------------------------------------
    % Note:
    % If you want to see, for example, PC1 and PC3, run it with [1 3]
    % instead of [1 2].
    % ----------------------------------------------------------------

%% Create loading plots (function 'VisLoadings1D')
VisLoadings1D(analysis.XAL1NS,analysis.PCA2logoff.loadings(1,:),analysis.ppmR)
    % ----------------------------------------------------------------
    % Note:
    % If you want to see, for example, PC2, use (2,:), instead of (1,:)
    % ----------------------------------------------------------------

%%  SECTION 8: SUMMARY
    % ----------------------------------------------------------------
    % Note:
    % This final section creates a summary for parameters used in this
    % workflow. This is very useful when others want to see the key
    % information of this workflow; they don't need to go through the whole
    % script.
    % ----------------------------------------------------------------
summary.removedRegionMin=data.ppmEndMin;
summary.removedRegionMax=data.ppmEndMax;
summary.removedWaterRegionMin=data.ppmWaterMin;
summary.removedWaterRegionMax=data.ppmWaterMax;
summary.alignmentMetric=data.distance_metric;
summary.alignmentMethod=data.method;
summary.normalization=data.normalization;
summary.scaling=analysis.scaling;

%%  Save Workspace
cd(paths.results);
save('matlab.mat')