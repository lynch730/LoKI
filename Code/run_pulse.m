

clear all; clc;
setup = Setup('pulseN200.in');
electronKinetics = setup.initializeSimulation();
electronKinetics.solve();
setup.finishSimulation();

clear all; clc;
setup = Setup('pulseN400.in');
electronKinetics = setup.initializeSimulation();
electronKinetics.solve();
setup.finishSimulation();

clear all; clc;
setup = Setup('pulseN600.in');
electronKinetics = setup.initializeSimulation();
electronKinetics.solve();
setup.finishSimulation();

clear all; clc;
setup = Setup('pulseN800.in');
electronKinetics = setup.initializeSimulation();
electronKinetics.solve();
setup.finishSimulation();

clear all; clc;
setup = Setup('pulseN1000.in');
electronKinetics = setup.initializeSimulation();
electronKinetics.solve();
setup.finishSimulation();

clear all; clc;
setup = Setup('pulseN1200.in');
electronKinetics = setup.initializeSimulation();
electronKinetics.solve();
setup.finishSimulation();

clear all; clc;
setup = Setup('pulseN1400.in');
electronKinetics = setup.initializeSimulation();
electronKinetics.solve();
setup.finishSimulation();

clear all; clc;
setup = Setup('pulseN1600.in');
electronKinetics = setup.initializeSimulation();
electronKinetics.solve();
setup.finishSimulation();

clear all; clc;
setup = Setup('pulseN1800.in');
electronKinetics = setup.initializeSimulation();
electronKinetics.solve();
setup.finishSimulation();

clear all; clc;
setup = Setup('pulseN2000.in');
electronKinetics = setup.initializeSimulation();
electronKinetics.solve();
setup.finishSimulation();

