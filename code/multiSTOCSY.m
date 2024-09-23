function [S_Table]= multiSTOCSY(XALN,ppm,peaklist,intensity_cutoff,corr_cutoff)
%{
   Chris Esselman
   5.18.23
   
   Function to do STOCSY on multiple driver peaks and automatically display
   correlated peaks to each driver peak in a struct

   Edit Log:
   Chris Esselman 5.31.23- Took off option to explicity say if peak picked
   or not. If the length of ppm is above a certain value indicating not peak
   picked, and the intensity_cutoff is zero will display a warning message
   that says the intensity_cutoff is higher

   Input    
              XALN = normalized data matrix

              ppm = chemical shift vector corresponding to XALN
 
              intensity_cutoff = the cutoff used to distinguish between
                                 background intensity and real peaks

              peaklist = the table of significant peaks to use when
              calculating STOCSY correlations for peaks

              corr_cutoff = correlation value cutoff when finding peaks
              that correlate with peak input used in STOCSY program
            
 Output      
               
               S_Table = Struct containing 1. initial peak list in first
               column. 2. ppm values of peaks highly correlated to each
               peak of peak list in second column. 3. correlation values
               of peaks correlated to inital driver peak.

               Next step is to export ppm values of peaks to format
               for easy database matching
%}

if (length(ppm) > 5000) && (intensity_cutoff == 0)
    warning(['If spectra is full resolution, the intensity cutoff should be ' ...
        'above 0']);
end
%Initialize struct and variable to index struct
S_Table = [];
place = 1;

%Loop through each important peak
for i=1:length(peaklist)

    %Create correlation vector
    [corr, ~]=STOCSY_NoFig_inside(peaklist(1,i),XALN,ppm);

    %Find ppm values and corresponding correlation values of peaks at or
    %above intensity cutoff and correlation cutoff
    [ppm_list,corr_list] = STOCSY_parser_inside(ppm,XALN,intensity_cutoff,corr,corr_cutoff);

    %Put driver peak, ppm values of highly correlated peaks, and
    %correlation values of the highly correlated peaks into struct
    S_Table(place).peak = peaklist(i);
    S_Table(place).ppm_correlated = ppm_list;
    S_Table(place).corr_values = corr_list;
    place = place + 1;
end
    function [corr,covar]=STOCSY_NoFig_inside(target,X,ppm)
        %
        %
        % This function is modified from STOCSY.
        % It does not generate a figure, just produces the outputs (GG-Mar-30-2020)
        % Code cleaned and more comments added (Chris Esselman 5.18.23)
        %
        % STOCSY(target,X,ppm)
        %
        %
        % Plots correlation/covariance projection of NMR spectrum to target
        % chemical shift or response vector

        % Arguments:
        %
        % target       Chemical shift of driver peak or response vector with length equal to size(X,1)
        % X            Data matrix of spectra
        % ppm          Chemical shift vector corresponding to X
        %close all

        %More comments located in STOCY with comments in ChrisE_zilla private
        %github directory but first portion just parses the inputs to ensure that
        %they are in the correct format.

        %Finds corresponding ppm index value for driver peak and calculates
        %correlation and covariance
        if length(target)==1
            [h,k]=min(abs(ppm-target));
            target_vect=X(:,k);
        else
            target_vect=target;
        end

        if size(target_vect,1)~=size(X,1) && size(target_vect,2)==size(X,1)
            target_vect=target_vect';
        end

        corr=(zscore(target_vect')*zscore(X))./(size(X,1)-1);
        covar=(target_vect-mean(target_vect))'*(X-repmat(mean(X),size(X,1),1))./(size(X,1)-1);
    end
    function [ppm_list,corr_list] = STOCSY_parser_inside(ppm,X_in,intensity_cutoff,corr,corr_cutoff)

        %{
    Edit Log
    Chirs Esselman 5.31.2023- Took off option to explicity say peak picked
    or not. Automatically does it. Change peak_pick_cutoff if not correct

    Main ideas for code comes from Sig_STOCSY, but much was incorrect.
    Ideas corrected and comments applied by Chris Esselman 5.18.23

 Inputs       
              ppm = chemical shift values corresponding to X matirx

              X_in = Data matrix

              intensity_cutoff = the cutoff used to distinguish between
                                 background intensity and real peaks

              corr = correlation values for driver peak created by STOCSY
              program
           
              corr_cutoff = correlation value cutoff when finding peaks
              that correlate with peak input used in STOCSY program

 Output       

            ppm_list = ppm values of peaks that are highly correlated to
            driver peak

            corr_list = correlation values to corresponding peaks of ppm
            list
                  
        %}

        %Create a modified average intensity table, which will help for calling peaks
        X_ave = mean(X_in,1);

        peak_picked_cutoff = 5000;
        %If peak picked, find peaks that are above intensity and correlation cutoff
        if length(ppm) < peak_picked_cutoff
            counter = 1;
            ppm_list = [];
            corr_list = [];
            for s = 1:size(X_in,2)
                if (abs(X_ave(1,s)) >= intensity_cutoff) && (corr(1,s) >= corr_cutoff)
                    ppm_list(1,counter)= ppm(1,s);
                    corr_list(1,counter)= corr(1,s);
                    counter = counter + 1;
                end
            end
            %If full resoltion, find local maxima that are above intensity and
            %correlation cutoff
        elseif length(ppm) > peak_picked_cutoff
            counter = 1;
            ppm_list = [];
            corr_list = [];
            for s = 2:(size(X_in,2)-1)
                if (X_ave(1,s) >= intensity_cutoff) && (corr(1,s) >= corr_cutoff)...
                        && (X_ave(1,s) > X_ave(1,s-1)) && (X_ave(1,s) > X_ave(s+1))
                    ppm_list(1,counter)= ppm(1,s);
                    corr_list(1,counter)= corr(1,s);
                    counter = counter + 1;
                end
            end
        end
    end
end