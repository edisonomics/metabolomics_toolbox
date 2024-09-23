function makeDatabaseInput(Multi_STOCSY_output,Filenametxt,which_database)
%{
makeDatabaseInput- put output of Multi_STOCSY into txt file for easy
                          copy and paste into 1D Query for Colmar. Txt file
                          is written into current working directory

Written by Chris Esselman 5.19.23
- Changed name for nmrbox push - 9.21.24

Input:
        Multi_Stocsy_output      variable in workspace after running
                                 Multi_Stocsy

        Filenametxt              String of desired file name
                                 ex/'STOCSY_file.txt'

Optional_Input:
    
        which_database          string- ex/ 'colmar'. Will be adding more
                                to this in the future

Output:

        Txt file saved to current working directory of multi_STOCSY output
        with each row being a different 1D trace

%}

arguments
    Multi_STOCSY_output
    Filenametxt string
    which_database string = 'colmar'
end

if strcmp(which_database,'colmar')
    %Want to print float with max of 2 digits before decimal and 4 digits after
    %decimal
    formatSpec = '%2.4f ';

    %Open the file to write to. Name of file is input name.
    Stocsy_file = fopen(Filenametxt,'w');

    % Parse through and combine peaks that are within 0.03ppm of each other
    for i = 1:length(Multi_STOCSY_output)
        vector_length = length(Multi_STOCSY_output(i).ppm_correlated);
        j = 1;
        subset_vector = [];
        while j <= vector_length
            start_index = j;
            while (j < vector_length) && ...
                    ((Multi_STOCSY_output(i).ppm_correlated(j + 1) - Multi_STOCSY_output(i).ppm_correlated(j)) <= 0.03)
                j = j + 1;
            end
            end_index = j;
            subset_median = median(Multi_STOCSY_output(i).ppm_correlated(start_index:end_index));
            subset_vector = [subset_vector, subset_median];
            j = j + 1;
        end
        Multi_STOCSY_output(i).ppm_correlated = subset_vector;
    end


    %Print each 1D trace as a row. Each ppm value seperated by a space
    for i= 1:length(Multi_STOCSY_output)
        fprintf(Stocsy_file,formatSpec,Multi_STOCSY_output(i).ppm_correlated);
        fprintf(Stocsy_file,'\n');
    end

    %Close the file and delete variable in workspace
    fclose(Stocsy_file);
    clear Stocsy_file
end
end