classdef rk4ZohTest < matlab.unittest.TestCase
    %rk4ZohTest Unit tests for rk4_zoh.m (fixed-step RK4 ZOH propagator).
    %  Reference solutions come from ode45 at tight tolerance on the same
    %  ode_descent dynamics with the control held constant.

    properties (Constant)
        Vc = 0.0777                    % effective Tsiolkovsky number
        x0 = [0.3; 1; -0.2; -0.6; 1]   % non-dim initial state
        u  = [0.5; 1.0]                % non-dim ZOH control
        dt = 0.4                       % non-dim ZOH interval
    end

    methods (TestClassSetup)
        function addHm2ToPath(testCase)
            hm2 = fileparts(fileparts(mfilename('fullpath')));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(hm2));
        end
    end

    methods (Test)
        function testMatchesOde45Reference(testCase)
            xRef = rk4ZohTest.ode45Reference();
            xRk4 = rk4_zoh(testCase.x0, testCase.u, testCase.dt, testCase.Vc, 8);
            testCase.verifyEqual(xRk4, xRef, 'AbsTol', 1e-8);
        end

        function testMassRowIsExact(testCase)
            % dm/dt = -Vc*|u| is constant under ZOH, so RK4 integrates the
            % mass row exactly for any number of substeps
            xRk4 = rk4_zoh(testCase.x0, testCase.u, testCase.dt, testCase.Vc, 1);
            mExpected = testCase.x0(5) - testCase.Vc*norm(testCase.u)*testCase.dt;
            testCase.verifyEqual(xRk4(5), mExpected, 'AbsTol', 1e-14);
        end

        function testFourthOrderConvergence(testCase)
            % Halving the substep size must cut the error by ~2^4
            xRef = rk4ZohTest.ode45Reference();
            err = arrayfun(@(n) norm(rk4_zoh(testCase.x0, testCase.u, ...
                testCase.dt, testCase.Vc, n) - xRef), [1 2 4 8]);
            order = log2(err(1:end-1) ./ err(2:end));
            testCase.verifyGreaterThan(min(order), 3.5);
        end
    end

    methods (Static)
        function xRef = ode45Reference()
            opts = odeset('RelTol', 1e-12, 'AbsTol', 1e-14);
            [~, Y] = ode45(@(t, x) ode_descent(x, rk4ZohTest.u, rk4ZohTest.Vc), ...
                           [0 rk4ZohTest.dt], rk4ZohTest.x0, opts);
            xRef = Y(end, :).';
        end
    end
end
