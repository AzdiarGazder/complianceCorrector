function [coeffs, crossOverPts, R2] = piecewiselm(xData, yData, n)
%% Function description:
% This function performs least squares-based piecewise linear regression
% modelling for n-segments. This is an extension of the 2-segment linear 
% regression model described in Bogartz (1968).
%
%% Author:
% Dr. Azdiar Gazder, 2025, azdiaratuowdotedudotau
%
%% Acknowledgements:
% Satoshi Okazaki
% For developing the original "piecewiselm.m" function, rev. Feb 2023.
%
% For details, please refer to the following reference:
% https://au.mathworks.com/matlabcentral/fileexchange/124190-piecewise-linear-model
%
% For the original reference, please refer to:
% R. S. Bogartz (1968): "A least squares method for fitting intercepting
% line segments to a set of data points." Psychol Bull 70(6), 749-755.
%
%% Syntax:
%  piecewiselm3(x, y, n)
%
%% Input:
%  xData    - @double, an [n,1] list of unique x-values
%  yData    - @double, an [n,1] list of corresponding y-values
%  n        - @double, a number corresponding to the number of segments
%
%% Output:
%  coeffs       - @double, a 2-column vector such that each row is a slope
%                 and intercept pair representing the linear regression
%                 coefficients of a segment
%  crossOverPts - @double, a 2-column vector such that each row is an x
%                 and y coordinate indicating where adjacent segments
%                 intersect
%  R2           - @double, a column vector of the coefficient of
%                 determination quantifying the overall goodness-of-fit of
%                 the piecewise linear model
%
%% Options:
%
%%


% Validate input
if n < 1
    error('The number of segments n must be at least 1.');
end

%% Pre-process: Aggregate data based on unique x-values
[xu, ~, ic] = unique(xData, 'sorted');
nX = numel(xu);

% Ensure there are enough unique points for n segments (each segment
% needs >=2 unique points)
if nX < 2*n
    error('Not enough unique data points for %d segments (need at least %d unique points).', n, 2*n);
end

% Aggregate statistics for each unique x-value
w     = accumarray(ic, 1);
sumY  = accumarray(ic, yData);
sumY2 = accumarray(ic, yData.^2);

% Cumulative sums for weighted regression error computation
Wcum   = cumsum(w);
WXcum  = cumsum(w .* xu);
WYcum  = cumsum(sumY);
WXXcum = cumsum(w .* xu.^2);
WXYcum = cumsum(xu .* sumY);
WYYcum = cumsum(sumY2);

% Pre-compute the cost matrix: costMatrix(ii,jj) is the weighted sum of
% squares error (SSE) for a segment spanning xu(i:j)
costMatrix = inf(nX, nX);
for ii = 1:nX-1
    for jj = ii+1:nX
        if ii == 1
            s_w  = Wcum(jj);
            s_x  = WXcum(jj);
            s_y  = WYcum(jj);
            s_xx = WXXcum(jj);
            s_xy = WXYcum(jj);
            s_yy = WYYcum(jj);
        else
            s_w  = Wcum(jj) - Wcum(ii-1);
            s_x  = WXcum(jj) - WXcum(ii-1);
            s_y  = WYcum(jj) - WYcum(ii-1);
            s_xx = WXXcum(jj) - WXXcum(ii-1);
            s_xy = WXYcum(jj) - WXYcum(ii-1);
            s_yy = WYYcum(jj) - WYYcum(ii-1);
        end

        % At least two unique points are required in the segment
        if (jj - ii + 1) < 2
            continue;
        end

        denom = s_w * s_xx - s_x^2;
        if abs(denom) < eps
            slopes = 0;
        else
            slopes = (s_w * s_xy - s_x * s_y) / denom;
        end
        intercepts = (s_y - slopes * s_x) / s_w;

        % Compute weighted SSE error for this segment
        err = s_yy - slopes * s_xy - intercepts * s_y;
        costMatrix(ii, jj) = err;
    end
end


%% Dynamic programming for optimal segmentation
% Let K = number of segments = n
K = n;
dp = inf(K, nX);
segPtr = zeros(K, nX);  % pointers for backtracking segmentation boundaries

% Base case: one segment from 1 to jj.
for jj = 2:nX
    dp(1, jj) = costMatrix(1, jj);
end

% Recurrence: for s segments and endpoint jj, try all valid previous
% breakpoints
for s = 2:K
    for jj = 2*s:nX
        for ii = 2*(s-1):jj-1
            candidate = dp(s-1, ii) + costMatrix(ii+1, jj);
            if candidate < dp(s, jj)
                dp(s, jj) = candidate;
                segPtr(s, jj) = ii;
            end
        end
    end
end

% Optimal total error is dp(K, nX)
M = dp(K, nX);

% Backtrack to recover segmentation boundaries (indices into xu)
boundaries = zeros(1, K+1);
boundaries(K+1) = nX;
for s = K:-1:2
    boundaries(s) = segPtr(s, boundaries(s+1));
end
boundaries(1) = 1;


%% Regression on the original Data for each segment & save segment bounds
% Use the segmentation (in terms of xu) to extract segments from the
% original data
coeffs = zeros(K, 2);
segBounds = zeros(K,2);  % store [min(x) max(x)] for each segment
for s = 1:K
    if s == 1
        segStart = 1;
        segEnd   = boundaries(2);
    else
        segStart = boundaries(s) + 1;
        segEnd   = boundaries(s+1);
    end
    xLower = xu(segStart);
    xUpper = xu(segEnd);
    % Select original data points that lie within [xLower, xUpper]
    xi = find(xData >= xLower & xData <= xUpper);
    xx = xData(xi);
    yy = yData(xi);
    lm = polyfit(xx, yy, 1);
    coeffs(s, :) = [lm(1), lm(2)]; % [slope, intercept]

    % Save segment boundaries from original data
    segBounds(s, :) = [min(xx), max(xx)];
    %     segBounds(s,1) = min(xx);
    %     segBounds(s,2) = max(xx);
end


%% Calculate the crossover points between consecutive segments
% For each adjacent pair of segments, compute the intersection of their
% regression lines
if K > 1
    slopes     = coeffs(:,1);
    intercepts = coeffs(:,2);
    crossOverPts = zeros(K-1, 2);
    for s = 1:K-1
        % Compute the intersection of line s and line s+1
        % (a1*x + b1 = a2*x + b2) => x = (b2 - b1)/(a1 - a2)
        denom = slopes(s) - slopes(s+1);
        if abs(denom) < eps
            % If slopes are nearly equal, use the mid-point of the segment
            % boundaries
            x_int = (segBounds(s,2) + segBounds(s+1,1)) / 2;
        else
            x_int = (intercepts(s+1) - intercepts(s)) / denom;
        end
        % Clamp the computed intersection so that it lies between the
        % maximum x of the left segment and the minimum x of the right
        % segment
        x_left  = segBounds(s,2);
        x_right = segBounds(s+1,1);
        if x_left > x_right
            % If segments overlap (which should rarely happen), fallback
            % to the average
            x_int = (x_left + x_right) / 2;
        else
            x_int = max(x_int, x_left);
            x_int = min(x_int, x_right);
        end
        y_int = slopes(s)*x_int + intercepts(s);
        % Place x and y values into a row vector
        crossOverPts(s, :) = [x_int, y_int];
    end
else
    crossOverPts = [];
end


%% Calculate the coefficient of determination (R²)
sumSquaresTotal = sum((yData - mean(yData)).^2);
R2 = 1 - M / sumSquaresTotal;

end
