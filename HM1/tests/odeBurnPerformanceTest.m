classdef odeBurnPerformanceTest < matlab.perftest.TestCase
    %odeBurnPerformanceTest Benchmarks for the HM1 hot loop.
    %  Unit: one ode_burn RHS evaluation (called ~1e6+ times by the fsolve
    %  sweeps of main_task1..4). System: one ode45 burn-arc integration at
    %  loose vs shooting-grade tolerances.
    %  Run with: results = runperf('odeBurnPerformanceTest')

    properties
        p    % costate/thrust parameter struct
        z0   % initial state [x; y; vx; vy; m; lam_m]
    end

    properties (TestParameter)
        RelTol = struct('Loose', 1e-6, 'Tight', 1e-10)
    end

    methods (TestClassSetup)
        function addHm1ToPath(testCase)
            hm1 = fileparts(fileparts(mfilename('fullpath')));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(hm1));
        end
    end

    methods (TestMethodSetup)
        function setupProblem(testCase)
            % Representative Task 1 solution neighbourhood (Q = 2.5, yf = 0.04)
            testCase.p  = struct('T', 1.5, 'Q', 2.5, 'c', 0.6, ...
                                 'lam_vx0', 0.6, 'lam_vy0', 3.8, 'lam_y', 14);
            testCase.z0 = [0; 0; 0; 0; 1; 1];
        end
    end

    methods (Test)
        function testRhsEvaluation(testCase)
            pp = testCase.p;
            z  = [0.05; 0.01; 0.4; 0.2; 0.7; 1.2];
            while testCase.keepMeasuring
                dz = ode_burn(0.15, z, pp);
            end
            testCase.verifySize(dz, [6 1]);
        end

        function testBurnArcIntegration(testCase, RelTol)
            pp   = testCase.p;
            ic   = testCase.z0;
            opts = odeset('RelTol', RelTol, 'AbsTol', RelTol*1e-2);
            while testCase.keepMeasuring
                [~, Z] = ode45(@(t,z) ode_burn(t, z, pp), [0 0.3], ic, opts);
            end
            % Mass equation is linear: exact propellant bookkeeping
            testCase.verifyEqual(Z(end,5), 1 - pp.Q*0.3, 'AbsTol', 1e-8);
        end
    end
end
