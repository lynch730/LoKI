
function run_pulse(N)

    if ~exist('N', 'var')
        N = 600;
    end

	setup = Setup(['pulse_N',sprintf('%i', N),'.in']);
	electronKinetics = setup.initializeSimulation();
	electronKinetics.solve();
	setup.finishSimulation();

end
