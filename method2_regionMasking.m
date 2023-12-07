clearvars
clc

reader = BioformatsImage('../data/W131R_Img1.nd2');
%reader = BioformatsImage('../data/WT_Img.nd2');

IGFP = getPlane(reader, 1, 'SDC-GFP', 1);
ICy5 = getPlane(reader, 1, 'SDC-Cy5', 1);
ITRITC = getPlane(reader, 1, 'SDC-TRITC', 1);


%% Make masks

maskGFP = makeMask(IGFP);
maskCy5 = makeMask(ICy5);
maskTRITC = makeMask(ITRITC);

% figure;
% imshow([imfuse(IGFP, maskGFP), imfuse(ICy5, maskCy5), imfuse(ITRITC, maskTRITC)])

%% Calculate correlation coeffs

%GFP-Cy5
rGFP_Cy5 = calculateCoeffs(maskCy5, ICy5, maskGFP, IGFP)
rTRITC_Cy5 = calculateCoeffs(maskTRITC, ITRITC, maskCy5, ICy5)


% 
% % %Make a mask
% % mask = imbinarize(I, 'adaptive', 'sensitivity', 0.8);
% % mask = imopen(mask, strel('diamond', 3));
% % mask = bwareaopen(mask, 100);
% % 
% % mask = imclearborder(mask);
% % 
% % imshowpair(I, mask)

function spotMask = makeMask(I)

d1 = imgaussfilt(I, 1/(sqrt(2)) * 3);
d2 = imgaussfilt(I, 1/(sqrt(2)) * 30);
diff = d1 - d2;

spotMask = diff > 1200;
spotMask = imclearborder(spotMask);

end

function rr = calculateCoeffs(mask1, I1, mask2, I2)

mask = mask1 | mask2;

dataI1 = double(I1(mask));
dataI2 = double(I2(mask));

rr = pearson(dataI1, dataI2);

end
