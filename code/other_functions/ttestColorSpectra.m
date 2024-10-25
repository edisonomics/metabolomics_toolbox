function [H,p,sigPvals] = ttestColorSpectra(X,ppm,y,groups,save,pval)
%% 
% Function performs 2 sample t-test between two groups
% and subsequently colors indidividual
% datapoints along spectrum to correspond with their p-values
%
% INPUTS    X = spectral datset
%           ppm = chemical shift vector
%           y = Yvector indicating group membership within 'X'
%           groups = 1 x 2 vector indicating which two members of 'y' you want
%               to test
%           save = input string 'save' to automatically save the
%               output figure to the current directory. Otherwise, put empty
%               string,''
%           pval = number indicating the p-value threshold for recording to
%               ouput argument "sigPvals"
%
% OUTPUTS   H = figure object of colorized spectrum 
%           p = list of p-values for each datapoint along 'ppm'
%           sigPvals = object containing significant (<=.05) pvalues and corresponding ppms 

% Modified from STOCSY plotting script 
% May 2017 by JW, GJG, MTJ, MBC

%% Simulated dataset for debugging
% X=rand(10,1000);
% ppm=linspace(0,10,1000);
% y=[0,0,0,0,0,1,1,1,1,1];
% ppmLeft=[2.31, 4.32];
% ppmRight=[2.50, 4.54];
% groups=[0,1];

%% Perform T-test 
% Groups specified in "groups" variable
group1=min(groups);
group2=max(groups);
pos=find(y==group2);
neg=find(y==group1);
for k=1:length(ppm)
   [h(k),p(k)]=ttest2(X(pos,k),X(neg,k));
end
%%
corr =  p; % list of pvalues
covar = mean(X,1); % list of covariance values

%% Generate list of only signifcant p values and corresponding ppms
sigPvals = [];
ind2 = 1;
for j = 1:length(p)
      if p(j)<=pval
         sigPvals(1,ind2) = ppm(j);
         sigPvals(2,ind2) = p(j);
         ind2 = ind2+1;
      end
end

if isempty(sigPvals)
    fprintf(['\n\n\tNo p-values below ',num2str(pval),' were found. sigPvals is empty.\n\n'])
end
%% set up colormap for plot

cmap = flipud(jet(100)); % reversing colormap so low pvalues show up as red                                            
lines=NaN(size(cmap,1),size(corr,2)); 

%% Modify these values to define range of pvales and interval.
pvalrange=[0 .2]; % the range of pvalues we want to see
interval = .05;   % determines both the number of intervals for calculating 'lines'                     
%%
ind=1;
for k=pvalrange(1):interval/size(cmap,1):pvalrange(2)  % subdividing by number of colors k = range(1):interval:range(2);
lines(ind,find(corr>k))=covar(find(corr>k));
ind=ind+1;
end
        
%% Generate the Figure
        H = figure;
        
                plot(ppm,lines(1,:),'Color',cmap(1,:),'LineWidth',1.1);
                hold on
                for k=2:size(cmap,1)
                plot(ppm,lines(k,:),'Color',cmap(k,:),'LineWidth',1.1); 
                end
                set(gca,'XDir','rev');
                xlabel('Chemical Shift (ppm)')
                title('P-values of spectral features between groups')
                H.Colormap = jet(100);

                c = colorbar;
                    c.Label.String = 'P-value';
                    c.Label.FontSize = 12;
                    
                    c.TickLabelsMode = 'manual';
                    c.TickLabels = fliplr(pvalrange(1):(pvalrange(2)/10):pvalrange(2));


       %% Save the Figure (if opted and if figure is actually generated)
        
               if   strcmp(save,'save') % only run if saving figures is desired
                    filename = ['pValues_for_groups_',num2str(groups(1)),'_',num2str(groups(2)),'_ttest.fig'];
                    savefig(H,filename,'compact')
                    fprintf('Figure was saved...\n');
%             close all
               else
                   fprintf('Figures will not be saved...\n');
               end
               
end