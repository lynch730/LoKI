

clear; clc;


setupFile = 'default_lokib_pulse_setup.in';
% save setupFile 
% clear all
% load setupfile

setup = Setup(setupFile);
electronKinetics = setup.initializeSimulation();

electronKinetics.solve();
setup.finishSimulation();
