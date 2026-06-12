classdef descentDynamicsPerformanceTest < matlab.perftest.TestCase
    %descentDynamicsPerformanceTest Performance benchmarks for the HM2 hot loop.
    %  Unit level: a single ode_descent evaluation and one rk4_zoh
    %  ZOH-interval propagation -- the building blocks of the trapezoidal
    %  and ZOH defect constraints that fmincon evaluates at every iteration.
    %  System level: the ode45 replay of one ZOH interval at the
    %  fidelity-check tolerances used by fwd_integrate.
    %  Run with: results = runperf('descentDynamicsPerformanceTest')

    properties (Constant)
        Vc = 0.0777    % effective Tsiolkovsky number (Table 1 data)
        dt = 0.0444    % one ZOH interval, non-dim (tf_nd / (N-1), N = 50)
    end

    properties
        x0
        u0
    end

    properties (TestParameter)
        nSub = struct('one', 1, 'two', 2, 'eight', 8)
    end

    methods (TestClassSetup)
        function addHm2ToPath(testCase)
            hm2 = fileparts(fileparts(mfilename('fullpath')));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(hm2));
        end
    end

    methods (TestMethodSetup)
        function setupState(testCase)
            % Representative mid-descent state and near-hover control
            testCase.x0 = [0.1; 0.5; -0.2; -0.4; 0.8];
            testCase.u0 = [0.1; 0.9];
        end
    end

    methods (Test)
        function testOdeDescentEvaluation(testCase)
            % A single call is ~60 ns -- below the framework precision --
            % so each measured sample is a batch of 1000 evaluations
            x = testCase.x0;  u = testCase.u0;  vc = testCase.Vc;
            while testCase.keepMeasuring
                for k = 1:1000
                    dx = ode_descent(x, u, vc);
                end
            end
            testCase.verifySize(dx, [5 1]);
        end

        function testRk4ZohPropagation(testCase, nSub)
            % Batch of 100 propagations per sample: the 1-substep case is
            % sub-microsecond and noise-dominated when measured singly
            x  = testCase.x0;  u = testCase.u0;
            vc = testCase.Vc;  h = testCase.dt;
            while testCase.keepMeasuring
                for k = 1:100
                    xn = rk4_zoh(x, u, h, vc, nSub);
                end
            end
            testCase.verifySize(xn, [5 1]);
        end

        function testOde45ZohReplayInterval(testCase)
            x  = testCase.x0;  u = testCase.u0;
            vc = testCase.Vc;  h = testCase.dt;
            opts = odeset('RelTol', 1e-10, 'AbsTol', 1e-12);
            while testCase.keepMeasuring
                [~, Y] = ode45(@(t, xx) ode_descent(xx, u, vc), [0 h], x, opts);
            end
            testCase.verifySize(Y(end,:), [1 5]);
        end
    end
end
