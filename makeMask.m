function spotMask = makeMask(I)
%MAKEMASK  Create a mask using the difference of Gaussians approach
%
%  M = MAKEMASK(I) will return a binary mask M after filtering the image
%  using the difference of Gaussians approach.

d1 = imgaussfilt(I, 1/(sqrt(2)) * 3);
d2 = imgaussfilt(I, 1/(sqrt(2)) * 30);
diff = d1 - d2;

spotMask = diff > 1200;
spotMask = imclearborder(spotMask);

end