clearvars
clc

%reader = BioformatsImage('../data/W131R_Img1.nd2');
reader = BioformatsImage('../data/WT_Img.nd2');

I = getPlane(reader, 1, 'SDC-GFP', 1);

%%

d1 = imgaussfilt(I, 1/(sqrt(2)) * 3);
d2 = imgaussfilt(I, 1/(sqrt(2)) * 50);
diff = d1 - d2;

spotMask = diff > 400;
spotMask = imclearborder(spotMask);

imshowpair(I, spotMask)

%%

dataGFP = double(I(spotMask));

ICy5 = getPlane(reader, 1, 'SDC-Cy5', 1);
ITRITC = getPlane(reader, 1, 'SDC-TRITC', 1);

dataCy5 = double(ICy5(spotMask));
dataTRITC = double(ITRITC(spotMask));

rGFP_Cy5 = pearson(dataGFP, dataCy5)
rTRITC_Cy5 = pearson(dataTRITC, dataCy5)


% %Make a mask
% mask = imbinarize(I, 'adaptive', 'sensitivity', 0.8);
% mask = imopen(mask, strel('diamond', 3));
% mask = bwareaopen(mask, 100);
% 
% mask = imclearborder(mask);
% 
% imshowpair(I, mask)