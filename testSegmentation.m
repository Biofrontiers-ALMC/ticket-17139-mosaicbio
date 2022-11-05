clearvars
clc

dataDir = 'D:\CU-Projects\mosaic-bio\data';

files = dir(fullfile(dataDir, '*.nd2'));

storeData = struct;

cropPx = 15;

for iFile = 1

    reader = BioformatsImage(fullfile(files(iFile).folder, files(iFile).name));

    %%
    Inucl = getPlane(reader, 1, 'DAPI', 1);
    IGFP = getPlane(reader, 1, 'EGFP', 1);
    
    %Crop the images
    Inucl = Inucl(cropPx:(size(Inucl, 1) - cropPx), ...
        cropPx:(size(Inucl, 2) - cropPx));
    IGFP = IGFP(cropPx:(size(IGFP, 1) - cropPx), ...
        cropPx:(size(IGFP, 2) - cropPx));
   
    maskNucl = imbinarize(Inucl, 'adaptive', 'sensitivity', 0.005);
    maskNucl = imopen(maskNucl, strel('disk', 3));
    maskNucl = bwareaopen(maskNucl, 50);

    dd = -bwdist(~maskNucl);
    dd(~maskNucl) = -Inf;
    dd = imhmin(dd, 2);

    L = watershed(dd);

    maskNucl(L == 0) = false;
    maskNucl = bwareaopen(maskNucl, 250);

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
    
    C = imfuse(IGFP, Inucl);
    showoverlay(C, maskCell, 'Color', [0 1 1], 'Opacity', 20);
    
%%

    


end







