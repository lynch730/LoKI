function solve(bz)
    % solve solves the electron Boltzmann equation according to the setup
    % specified by the user
    
    % logging start of the boltzmann calculations
    start = tic;
    notify(bz, 'genericStatusMessage', ...
           StatusEventData('\t- Solving Boltzmann ...\n', 'status'));
    
    % select the time-dependent or time-independent solutions
    if bz.isTimeDependent
        bz.obtainTimeDependentSolution();
    else
        bz.obtainTimeIndependentSolution();
    end
    
    % logging end of the boltzmann calculations
    str = sprintf('\\t    Finished (%f seconds).\\n', toc(start));
    notify(bz, 'genericStatusMessage', StatusEventData(str, 'status'));

end