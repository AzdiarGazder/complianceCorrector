function [xNew, yNew, varargout] = calcPLRM(xData, yData, n, varargin)

%% Pre-define options
% polyOrder = get_option(varargin, 'order', ones(n,1));
fontSize = get_option(varargin, 'fontSize', 14);
markerSize = get_option(varargin, 'markerSize', 10);
lineWidth = get_option(varargin, 'lineWidth', 1.5);
xAxisLabel = get_option(varargin, 'xLabel', 'xData');
yAxisLabel = get_option(varargin, 'yLabel', 'yData');

% polyOrder = polyOrder(:);
% if numel(polyOrder) ~= n
%     error('The number of polynominal orders should be equal to n');
%     return;
% end

% %% Validate polyOrder size
% if length(polyOrder) ~= n
%     error('The size of polyOrder must be [n 1] where n is the number of segments.');
% end

%% Perform piecewise linear regression modelling
% [~, crossOverPts, R2] = piecewiselm(xData, yData, n);
% slopes = coeffs(:,1); intercepts = coeffs(:,2); GOF = R2;
[coeffs, crossOverPts, ~] = piecewiselm(xData, yData, n);

% Plot the input data
figure;
plot(xData, yData, '.k', 'MarkerSize', markerSize);
xlabel(xAxisLabel, 'FontSize', fontSize);
ylabel(yAxisLabel, 'FontSize', fontSize);
grid on;
hold all;

% Find the closest crossover points
xc = crossOverPts(:,1);
deltax = abs(xData - xc(:)');
[~, closestIdx] = min(deltax, [], 1);
xClosest = xData(closestIdx);
yClosest = yData(closestIdx);

% Plot a green circle around the closest crossover points
plot(xClosest, yClosest, 'ko', 'LineWidth', lineWidth, 'MarkerSize', markerSize);

% Initialise variables to store fitted data
xNew = [];
yNew = [];
varargout = cell(1, n);

% Loop through the segments based on the number of crossover points (n)
for ii = 1:n
    if ii == 1
        % First segment
        xSegment = xData(1:closestIdx(ii)-1);
        ySegment = yData(1:closestIdx(ii)-1);
    elseif ii == n
        % Last segment
        xSegment = xData(closestIdx(ii-1):end);
        ySegment = yData(closestIdx(ii-1):end);
    else
        % Middle segments
        xSegment = xData(closestIdx(ii-1):closestIdx(ii)-1);
        ySegment = yData(closestIdx(ii-1):closestIdx(ii)-1);
    end

    % % Fit the polynomial based on the degree in polyOrder(ii)
    % p = polyfit(xSegment, ySegment, polyOrder(ii));
    % yFitSegment = polyval(p, xSegment);

    % Fit the polynomial based on the best-fit coefficients from PLRM
    yFitSegment = polyval(coeffs(ii,:), xSegment);

    % If this is not the first segment, adjust the fitted line to maintain continuity
    if ii > 1
        yFitSegment = yFitSegment - yFitSegment(1);
        yFitSegment = yNew(end) + yFitSegment;
    end

    % Concatenate the fitted data for each segment
    xNew = [xNew; xSegment];
    yNew = [yNew; yFitSegment];

    % Plot the fitted segment
    plot(xSegment, yFitSegment, 'LineWidth', lineWidth);

    % Store each segment's data in varargout
    segmentData.x = xSegment;
    segmentData.y = yFitSegment;
    % segmentData.P = p;  % Store the best-fit polynomial coefficients
    segmentData.P = coeffs;  % Store the PLRM best-fit polynomial coefficients
    varargout{ii} = segmentData;
end

end
