function status = odeProgressBar(t,~,flag,varargin)
    
    persistent progressFigure;
    persistent progressGraph;
    persistent integrationTimeStr;
    persistent progressBar;
    persistent initialClock;
    
    switch(flag)
        case 'firstInit'
            firstTimeStep = varargin{1};
            finalTime = varargin{2};
            screenSize = get(groot,'ScreenSize');
            progressFigure = figure('Name', 'Pulse temporal integration progress bar', 'NumberTitle', 'off', ...
                'MenuBar', 'none', ...
                'Position', [floor(screenSize(3)/4) floor(screenSize(4)*5/12) floor(screenSize(3)/2) floor(screenSize(4)/6)]);
            integrationTimeStr = uicontrol('Parent', progressFigure, 'Style', 'text', 'Units', 'normalized', ...
                'Position', [0.11 0.65 0.89 0.15], 'HorizontalAlignment', 'left', 'String', 'Computational time: 0 s');
            progressGraph = axes('Parent', progressFigure, 'Units', 'normalized', 'OuterPosition', [0 0 1 0.6], ...
                'Box', 'on', 'XScale', 'log', 'Ytick', [], 'Xlim', [firstTimeStep finalTime]);
            xlabel('Integration time (s)');
            hold on;
            progressBar = area(progressGraph, [firstTimeStep firstTimeStep], [1 1]);
            initialClock = clock;
        case 'lastDone'
            close(progressFigure)
            clear progressBar integrationTimeStr progressGraph progressFigure initialClock
        case 'init'
    
        case 'done'
    
        otherwise
            progressBar.XData = [1e-18 t(end)];
            integrationTimeStr.String = sprintf('Computational time: %.1f s', etime(clock, initialClock));
    end
    status = 0;
    drawnow;

end