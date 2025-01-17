function [coeffs,crossOverPts,R2] = piecewiselm(xData,yData,n)

% PIECEWISELM performs piecewise linear regression with least squares
%   [coef,breakPt,R2] = piecewiselm(x,y,n)
%   performs n-segmented linear regression. This is an extension of
%   two-segmented linear regression described in Bogartz (1968).
%
%   R. S. Bogartz (1968): "A least squares method for fitting intercepting
%   line segments to a set of data points." Psychol Bull 70(6), 749-755.
%
%   x      : independent variable
%   y      :   dependent variable
%   n      : number of segments
%   coef   : [a1 a0 b1 b0 ...] y=a1x+a0; y=b1x+b0; ...
%   breakPt: [x0 y0 x1 y1 ...]
%   R2     : R-squared
%
%(c) S. Okazaki 2023.02
%% Create segment table
% Segment:
% a0-a1 b0-b1 c0-c1 ...
%
% Table:
% a0 a1 b0 b1 c0 c1 ...
% a0 a1 b0 b1 c0 c1 ...
% ...


%% Organise the data
xData_ordered = unique(xData,'sorted');
nX = length(xData_ordered);
endTable = nchoosek(2:nX-2,n-1);
[nRows,nCols] = size(endTable);
flags = ones(nRows,1);
if nCols > 1
    for ii = 1:nRows
        for jj = 1:nCols-1
            a = endTable(ii,jj);
            b = endTable(ii,jj+1);
            if b-a < 2
                flags(ii) = 0;
                break;
            end
        end
    end
end
ii = flags == 1;
endTable = endTable(ii,:);
startTable = endTable + 1;
nRows = size(endTable,1);
startTable = [ones(nRows,1) startTable];
endTable = [endTable ones(nRows,1)*nX];
nCols = size(endTable,2);
indexTable = zeros(nRows,2*nCols);

for ii = 1:nCols
    indexTable(:,2*ii-1) = startTable(:,ii);
    indexTable(:,2*ii)   =   endTable(:,ii);
end


%% Regression
nRows = size(indexTable,1);
sumSquares = zeros(nRows,1);
bTable = zeros(nRows,n*2);

for ii = 1:nRows
    for jj = 1:n
        i0 = indexTable(ii,2*jj-1);
        i1 = indexTable(ii,2*jj);
        x0 = xData_ordered(i0);
        x1 = xData_ordered(i1);
        xi = find(xData >= x0 & xData <= x1);
        xx = xData(xi);
        yy = yData(xi);
        lm      = polyfit(xx,yy,1);
        yyFit   = polyval(lm,xx);
        yyResidual = yy - yyFit;
        sumSquaresResidual = sum(yyResidual.^2);

        sumSquares(ii) = sumSquares(ii) + sumSquaresResidual;
        b1 = lm(1);
        b0 = lm(2);
        bTable(ii,2*jj-1) = b1;
        bTable(ii,2*jj)   = b0;
    end
    progress(ii,nRows);
end
[M,I] = min(sumSquares);
Ms = find(sumSquares == M);
if length(Ms) > 1
    warning('Two or more minimum sums of squares confirmed.')
end
coeffs = bTable(I,:);
crossOverPts = zeros(1,2*(n-1));
for ii = 1:n-1
    a1 = coeffs(2*ii-1);
    a0 = coeffs(2*ii);
    b1 = coeffs(2*(ii+1)-1);
    b0 = coeffs(2*(ii+1));
    crossOverPtX = (b0-a0) / (a1 - b1);
    crossOverPtY = a1 * crossOverPtX + a0;
    crossOverPts(2*ii-1) = crossOverPtX;
    crossOverPts(2*ii)   = crossOverPtY;
end
sumSquaresTotal = sum((yData - mean(yData)).^2);
R2 = 1 - M / sumSquaresTotal;

end