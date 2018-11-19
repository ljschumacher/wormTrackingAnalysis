% script finds all the .eps files in the directory, converts them to
% .pdf files, and removes the original .eps files. 
% script is made because some workstations don't have the epstopdf function
% so analysis figures are originally saved as .eps files.

% find all the .eps files
eps = dir('*.eps');
% loop through each file
for i = 1:length(eps)
    % get the file name
    epsname = eps(i).name;
    % convert file to .pdf
    system(['epstopdf ' epsname]);
    % remove .eps file
    system(['rm ' epsname]);
end
