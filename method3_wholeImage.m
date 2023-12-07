clearvars
clc

%folder = 'D:\Work\ALMC\mosaic-bio-2\data\Images for analysis - W131R';
folder = 'D:\Work\ALMC\mosaic-bio-2\data\Images for analysis - WT';

files = dir(fullfile(folder, '*.nd2'));

fid = fopen(fullfile(folder, 'correlation.csv'), 'w');

fprintf(fid, 'Filename, GFP-Cy5, TRITC-Cy5\n');

for iFile = 1:numel(files)

    %reader = BioformatsImage('../data/W131R_Img1.nd2');
    reader = BioformatsImage(fullfile(files(iFile).folder, files(iFile).name));

    try
        IGFP = double(getPlane(reader, 1, 'SDC-GFP', 1));
        ICy5 = double(getPlane(reader, 1, 'SDC-Cy5', 1));
        ITRITC = double(getPlane(reader, 1, 'SDC-TRITC', 1));

    catch ME
        
        fprintf(fid, '%s, %s, \n', files(iFile).name, ME.message);
        continue

    end
    %Calculate correlation coeffs

    rGFP_Cy5 = pearson(IGFP, ICy5);
    rTRITC_Cy5 = pearson(ITRITC, ICy5);

    fprintf(fid, '%s, %.3f, %.3f\n', files(iFile).name, rGFP_Cy5, rTRITC_Cy5);

end

fclose(fid);