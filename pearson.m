function r = pearson(A, B)
%PEARSON  Compute the Pearson correlation coefficient
%
%  R = PEARSON(A, B) returns the Pearson correlation coefficient.

meanA = mean(A, 'all');
meanB = mean(B, 'all');

Asub = A(:) - meanA;
Bsub = B(:) - meanB;

r = sum(Asub .* Bsub)/(sqrt(sum(Asub.^2) .* sum(Bsub.^2)));
