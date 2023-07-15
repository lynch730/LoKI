

clear; clc;

setupFile = 'default_lokib_pulse_setup.in';

clear Boltzmann Collision EedfGas EedfState Output
setup = Setup(setupFile);
electronKinetics = setup.initializeSimulation();

electronKinetics.solve();
setup.finishSimulation();

