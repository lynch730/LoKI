function eedf = matrixInversion(matrix, energyGrid)
    % matrixInversion introduce the normalization condition intro de matrix in the arguments and obtain the solution by
    % direct inverion (mldivide matlab routine) and returns the corresponding eedf
    
    % local copies of energy grid variables
    energyCell = energyGrid.cell;
    energyStep = energyGrid.step;
    cellNumber = energyGrid.cellNumber;
    
    % include normalization condition for the EEDF in the Boltzmann matrix
    matrix(1,:) = matrix(1,:) + energyStep*sqrt(energyCell);
    
    % invert the Boltzmann matrix
    eedf = (matrix\([1 zeros(1,cellNumber-1)]'))';
    
    % renormalize of the EEDF
    eedf = eedf/dot(eedf,sqrt(energyCell)*energyStep);

end