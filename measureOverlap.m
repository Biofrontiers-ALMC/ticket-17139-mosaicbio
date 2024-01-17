clearvars
clc

%fileDir = 'D:\Work\ALMC\mosaic-bio-2\data\Images for analysis - W131R';
fileDir = 'D:\Work\ALMC\mosaic-bio-2\data\Images for analysis - WT';
file = 'WT_Img038.nd2';

reader = BioformatsImage(fullfile(fileDir, file));

ICy5 = getPlane(reader, 1, 'SDC-Cy5', 1);

IGFP = getPlane(reader, 1, 'SDC-GFP', 1);
ITRITC = getPlane(reader, 1, 'SDC-TRITC', 1);
%%
maskCy5 = makeMask(ICy5);
maskTRITC = makeMask(ITRITC);
maskGFP = makeMask(IGFP);

%Look at overlapping masks?
areaOverlap_Cy5_TRITC = maskCy5 & maskTRITC;
areaOverlap_Cy5_GFP = maskCy5 & maskGFP;

imshow(areaOverlap_Cy5_TRITC)

%Blob analysis - let's do TRITC intensity within each Cy5 region
objData = regionprops(maskCy5, ITRITC, 'MaxIntensity');
TRITCints = cat(1, objData.MaxIntensity);

numLocalizedTRITC = nnz(TRITCints > prctile(ITRITC(:), 99));

objData = regionprops(maskCy5, IGFP, 'MaxIntensity');
GFPints = cat(1, objData.MaxIntensity);

numLocalizedGFP = nnz(GFPints > prctile(IGFP(:), 99));

%Return stats
pcOverlap_Cy5_TRITC = (nnz(areaOverlap_Cy5_TRITC)/nnz(maskCy5)) * 100
pcOverlap_Cy5_GFP = (nnz(areaOverlap_Cy5_GFP)/nnz(maskCy5)) * 100

pcTRITClocalized = (numLocalizedTRITC / numel(objData)) * 100
pcGFPocalized = (numLocalizedGFP / numel(objData)) * 100

%Number of blobs
%%
figure(1)
imshowpair(maskCy5, ICy5)

figure(2)
imshowpair(maskTRITC, ITRITC)

%%
Iout = showoverlay(maskCy5, maskTRITC, 'Opacity', 80);
showoverlay(Iout, maskGFP, 'Color', [1 0 1], 'Opacity', 80)

function spotMask = makeMask(I)

d1 = imgaussfilt(I, 1/(sqrt(2)) * 3);
d2 = imgaussfilt(I, 1/(sqrt(2)) * 30);
diff = d1 - d2;

spotMask = diff > 1200;
spotMask = imclearborder(spotMask);

end