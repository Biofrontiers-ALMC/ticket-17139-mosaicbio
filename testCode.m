clearvars
clc

dataDir = 'D:\Projects\ALMC Tickets\T17139-TobinBrown\data';

files = dir(fullfile(dataDir, '*.nd2'));

storeData = struct;

for iFile = 1:numel(files)

    reader = BioformatsImage(fullfile(files(iFile).folder, files(iFile).name));

    %reader = BioformatsImage(fullfile(dataDir, 'WT_highMOPC.nd2'));
    %reader = BioformatsImage(fullfile(dataDir, 'W131R_highMOPC.nd2'));

    %%
   Inucl = getPlane(reader, 1, 'DAPI', 1);
    
    IGFP = getPlane(reader, 1, 'EGFP', 1);

    maskNucl = imbinarize(Inucl, 'adaptive', 'sensitivity', 0.005);
    maskNucl = imopen(maskNucl, strel('disk', 3));
    maskNucl = bwareaopen(maskNucl, 50);

    dd = -bwdist(~maskNucl);
    dd(~maskNucl) = -Inf;
    dd = imhmin(dd, 2);

    L = watershed(dd);

    maskNucl(L == 0) = false;
    maskNucl = bwareaopen(maskNucl, 250);
%         imshowpair(Inucl, mask)


%     imshowpair(Inucl, maskNucl)
%%
    %Cell mask
    maskCell = imbinarize(IGFP, 'adaptive', 'sensitivity', 0.7);
    maskCell = imopen(maskCell, strel('diamond', 3));
    maskCell = bwareaopen(maskCell, 100);

%    maskCell = maskCell(2:(end - 1), 2:(end - 1));
    maskCell = imclearborder(maskCell);

    ddCell = -bwdist(~maskCell);
    ddCell(~maskCell) = -Inf;
    dd = imhmin(dd, 2);

    dd = imimposemin(dd, maskNucl);

    L = watershed(dd);
    maskCell(L == 0) = false;
    
    maskCell = bwareaopen(maskCell, 100);
    
    maskCell(maskNucl) = false;


    C = imfuse(IGFP, Inucl);
    showoverlay(C, bwperim(maskCell), 'Color', [1 1 0], 'Opacity', 100);


    %%
    %Measure cell properties
    Igfp = double(getPlane(reader, 1, 'EGFP', 1));
    Igfp_rand = Igfp;
    Igfp_rand = Igfp_rand(randperm(numel(Igfp_rand)));
    Igfp_rand = reshape(Igfp_rand, size(Igfp));

    Irfp = double(getPlane(reader, 1, 'TRITC', 1));

    Icy5 = double(getPlane(reader, 1, 'Cy5', 1));

    cellData = regionprops(maskCell, Inucl, 'MeanIntensity', 'PixelIdxList');

    %Exclude misidentified regions

    %The main outputs I would like to see are 1) colocalization between
    %(channel 3 and channel 2) and (channel 3 and channel 4), and the total
    %intensity per cell of channel 3.


    %Measure correlation between EGFP and TRITC channel
    for iCell = 1:numel(cellData)

        cellData(iCell).pccscore_rfpvgfp = pearson(Irfp(cellData(iCell).PixelIdxList), ...
             Igfp(cellData(iCell).PixelIdxList));

        cellData(iCell).pccscore_rfpvcy5 = pearson(Irfp(cellData(iCell).PixelIdxList), ...
             Icy5(cellData(iCell).PixelIdxList));


%         cellData(iCell).pccscore_rfpvgfp = corrcoef(Irfp(cellData(iCell).PixelIdxList), ...
%             Igfp(cellData(iCell).PixelIdxList));
% 
%         %     %Randomly permute the vector
%         %     rfpdata = Irfp(cellData(iCell).PixelIdxList);
%         %     rfpdata = rfpdata(randperm(numel(rfpdata)));
%         %
%         %     cellData(iCell).pccnull = corrcoef(Igfp(cellData(iCell).PixelIdxList), ...
%         %         rfpdata);
% 
%         cellData(iCell).pccscore_rfpvcy5 = corrcoef(Irfp(cellData(iCell).PixelIdxList), ...
%             Icy5(cellData(iCell).PixelIdxList));

        cellData(iCell).meanTritc = mean(Irfp(cellData(iCell).PixelIdxList));

    end

    storeData(iFile).filename = files(iFile).name;
    storeData(iFile).data = cellData;

end

return

%% Data analysis

for iF = 1:numel(storeData)

    histogram([storeData(iF).data.pccscore_rfpvcy5], ...
        'BinWidth', 0.05, ...
        'Normalization', 'probability')
    hold on

end
hold off

legend({storeData.filename}, 'Interpreter', 'none')

xlabel('Pearson''s correlation coefficient')
ylabel('Normalized counts')
title('Ch3 v Ch4')


%%







