
%% Fit: 'untitled fit 1'.
[xData, yData] = prepareCurveData( T, C );

% Set up fittype and options.
ft = fittype( 'smoothingspline' );

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft );

% Plot fit with data.
figure( 'Name', 'untitled fit 1' );
h = plot( fitresult, xData, yData );
legend( h, 'C vs. T', 'untitled fit 1', 'Location', 'NorthEast' );
% Label axes
xlabel T
ylabel C
grid on


