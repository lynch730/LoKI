function eedf = linearSolver(bz)
    % linearSolver invert the matrix of the discretized Boltzmann equation
    % without considering non-linear operators (non-conservative ionization
    % or attachment, growth models, e-e collisions) in order to obtain an
    % eedf
    
    % sum all contributions to the boltzmann matrix and reescale it to
    % avoid very small numbers
    matrix =    bz.workCond.reducedFieldSI^2.0 * bz.fieldMatrix ...
              + bz.elasticMatrix ...
              + bz.CARMatrix ...
              + bz.inelasticMatrix ...
              + bz.ionizationConservativeMatrix ...
              + bz.attachmentConservativeMatrix;
    matrix = 1.0e20 .* matrix;

    % invert the matrix
    eedf = matrixInversion(matrix, bz.energyGrid);
    
    % store eedf in the properties of the boltzmann object
    bz.eedf = eedf;

end