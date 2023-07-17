
function run_pulse(N)

	setup = Setup(['pulseN',sprintf('%i', N),'.in']);
	electronKinetics = setup.initializeSimulation();
	electronKinetics.solve();
	setup.finishSimulation();

end
