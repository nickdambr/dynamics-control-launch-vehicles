classdef odeDescentTest < matlab.unittest.TestCase
    %odeDescentTest Unit tests for ode_descent.m.
    %  Non-dim throughout (g = 1, state [x; y; vx; vy; m]).

    methods (TestClassSetup)
        function addHm2ToPath(testCase)
            hm2 = fileparts(fileparts(mfilename('fullpath')));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(hm2));
        end
    end

    methods (Test)
        function testDerivativeDefinition(testCase)
            % Hand-computed derivative, |u| = 1
            x  = [1; 2; 0.3; -0.4; 0.8];
            u  = [0.6; 0.8];                 % |u| = 1
            Vc = 0.0777;
            dx = ode_descent(x, u, Vc);
            expected = [0.3; -0.4; 0.6/0.8; 0.8/0.8 - 1; -Vc];
            testCase.verifyEqual(dx, expected, 'AbsTol', 1e-15);
        end

        function testBallisticCoast(testCase)
            % Zero thrust: free fall, no mass flow
            x  = [0.5; 1; 0.2; -0.1; 0.9];
            dx = ode_descent(x, [0; 0], 0.0777);
            testCase.verifyEqual(dx, [0.2; -0.1; 0; -1; 0], 'AbsTol', 1e-15);
        end

        function testHoverEquilibrium(testCase)
            % Ty = m cancels gravity: vx, vy derivatives are zero
            m  = 0.73;
            x  = [0; 1; 0; 0; m];
            dx = ode_descent(x, [0; m], 0.0777);
            testCase.verifyEqual(dx(3), 0, 'AbsTol', 1e-15);
            testCase.verifyEqual(dx(4), 0, 'AbsTol', 1e-15);
        end

        function testMassFlowDependsOnlyOnThrustMagnitude(testCase)
            % dm = -Vc*|u|: independent of thrust direction
            x  = [0; 1; 0; 0; 1];
            Vc = 0.0777;
            dx1 = ode_descent(x, [1; 0],  Vc);
            dx2 = ode_descent(x, [0; -1], Vc);
            testCase.verifyEqual(dx1(5), -Vc,    'AbsTol', 1e-15);
            testCase.verifyEqual(dx2(5), dx1(5), 'AbsTol', 1e-15);
        end
    end
end
