clearvars
clc

%reader = BioformatsImage('../data/W131R_Img1.nd2');
%reader = BioformatsImage('../data/WT_Img.nd2');

folder = 'D:\Work\ALMC\mosaic-bio-2\data\Images for analysis - W131R';
%folder = 'D:\Work\ALMC\mosaic-bio-2\data\Images for analysis - WT';

files = dir(fullfile(folder, '*.nd2'));

fid = fopen(fullfile(folder, 'correlation_maskedCy5.csv'), 'w');

fprintf(fid, 'Filename, GFP-Cy5, TRITC-Cy5\n');

for iFile = 1:numel(files)

    reader = BioformatsImage(fullfile(files(iFile).folder, files(iFile).name));

    try
        IGFP = getPlane(reader, 1, 'SDC-GFP', 1);
        ICy5 = getPlane(reader, 1, 'SDC-Cy5', 1);
        ITRITC = getPlane(reader, 1, 'SDC-TRITC', 1);

    catch ME

        fprintf(fid, '%s, %s, \n', files(iFile).name, ME.message);
        continue

    end

    maskCy5 = makeMask(ICy5);

    %Make some output images
    IGFPnorm = double(IGFP);
    IGFPnorm = (IGFPnorm - min(IGFPnorm(:)))/(max(IGFPnorm(:)) - min(IGFPnorm(:)));

    ICy5norm = double(ICy5);
    ICy5norm = (ICy5norm - min(ICy5norm(:)))/(max(ICy5norm(:)) - min(ICy5norm(:)));

    ITRITCnorm = double(ITRITC);
    ITRITCnorm = (ITRITCnorm - min(ITRITCnorm(:)))/(max(ITRITCnorm(:)) - min(ITRITCnorm(:)));

    Iout1 = imfuse(IGFPnorm, ICy5norm);
    Iout1 = showoverlay(Iout1, bwperim(maskCy5), 'Color', [1 1 0]);
    Iout1 = insertText(Iout1, [11, 12], 'GFP-Cy5', 'TextColor', 'white', 'BoxOpacity', 0);

    Iout2 = imfuse(ITRITCnorm, ICy5norm);
    Iout2 = showoverlay(Iout2, bwperim(maskCy5), 'Color', [1 1 0]);
    Iout2 = insertText(Iout2, [11, 12], 'TRITC-Cy5', 'TextColor', 'white', 'BoxOpacity', 0);

    imwrite([Iout1 Iout2], ...
        fullfile(files(iFile).folder, [files(iFile).name(1:end - 4), '.png']))

    %Measure the Pearson correlation coeff within the mask
    rGFP_Cy5 = calculateCoeffs(maskCy5, IGFP, ICy5);
    rTRITC_Cy5 = calculateCoeffs(maskCy5, ITRITC, ICy5);

    fprintf(fid, '%s, %.3f, %.3f\n', files(iFile).name, rGFP_Cy5, rTRITC_Cy5);

end

fclose(fid);

function spotMask = makeMask(I)

d1 = imgaussfilt(I, 1/(sqrt(2)) * 3);
d2 = imgaussfilt(I, 1/(sqrt(2)) * 30);
diff = d1 - d2;

spotMask = diff > 1200;
spotMask = imclearborder(spotMask);

end

function rr = calculateCoeffs(mask, I1, I2)

dataI1 = double(I1(mask));
dataI2 = double(I2(mask));

rr = pearson(dataI1, dataI2);

end