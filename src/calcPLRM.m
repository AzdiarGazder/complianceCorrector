function [xNew,yNew] = calcPLRM(xData, yData, n, varargin)

%% Pre-define options
fontSize = get_option(varargin,'fontSize',14);
markerSize = get_option(varargin,'markerSize',10);
lineWidth = get_option(varargin,'lineWidth',0.5);


%% Perform piecewise linear regression modelling
[coeffs,crossOverPts,~] = piecewiselm(xData,yData,n);


slope = coeffs(1:2:end);
intercept = coeffs(2:2:end);

% Plot the input data
figure;
plot(xData, yData, '.k', 'MarkerSize', markerSize);
xlabel('Time (s)', 'FontSize', fontSize);
ylabel('Displacement (mm)', 'FontSize', fontSize);
grid on;
hold all;


% Find the closest crossover points
xc = crossOverPts(1:2:end);
for ii = 1:length(xc)
    deltax = abs(xData - xc(ii));
    [~, closestIdx(ii)] = min(deltax);
    xClosest(ii) = xData(closestIdx(ii));
    yClosest(ii) = yData(closestIdx(ii));
end
% Plot a green circle around the closest crossover points
plot(xClosest, yClosest, 'ko', 'LineWidth', lineWidth, 'MarkerSize', markerSize);

% Overlay the best fit lines on the data
% disp('...')
% for ii = 1:3
%     % Use the equation of line to get fitted/regressed y values
%     yFit(:,ii) = slope(ii) * xData + intercept(ii);
%     if sign(intercept(ii)) >=0
%         message = sprintf('Equation (%1.0f): y = %.16f * x + %.16f', ii, slope(ii), abs(intercept(ii)));
%     else
%         message = sprintf('Equation (%1.0f): y = %.16f * x - %.16f', ii, slope(ii), abs(intercept(ii)));
%     end
%     fprintf('%s\n', message);
% end
% disp('...')

% Calculate a new best fit polynominal for the first segment
x1 = xData(1:closestIdx(1)-1);
y1 = yData(1:closestIdx(1)-1);
p1 = polyfit(x1,y1,2);
yF1 = p1(1).*x1.^2 + p1(2).*x1 + p1(3);
plot(x1,yF1,'-r')

% Calculate a new best fit line for the second segment
x2 = xData(closestIdx(1):closestIdx(2)-1);
y2 = yData(closestIdx(1):closestIdx(2)-1);
p2 = polyfit(x2,y2,1);
yF2 = p2(1).*x2 + p2(2);
plot(x2,yF2,'-g')

% Calculate a new best fit line for the third segment
x3 = xData(closestIdx(2):end);
y3 = yData(closestIdx(2):end);
p3 = polyfit(x3,y3,1);
yF3 = p3(1).*x3 + p3(2);
plot(x3,yF3,'-b')

% Concantenate the fitted data
xNew = [x1; x2; x3];
yNew = [yF1; yF2; yF3];

end
