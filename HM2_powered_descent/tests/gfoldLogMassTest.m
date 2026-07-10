classdef gfoldLogMassTest < matlab.unittest.TestCase
    %gfoldLogMassTest Unit tests for the GFOLD log-mass kernels.
    %  Covers ode_descent_uacc.m (acceleration-ZOH RHS) and lti_zoh.m (exact
    %  LTI ZOH discretisation). Non-dim throughout (g = 1).

    properties (Constant)
        Vc = 0.0777;       % V_ref/c
        dt = 0.0444;       % tf_nd/(N-1) with N = 50
    end

    methods (TestClassSetup)
        function addHm2ToPath(testCase)
            hm2 = fileparts(fileparts(mfilename('fullpath')));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(hm2));
        end
    end

    methods (Test)
        % ---- ode_descent_uacc.m -------------------------------------------
        function testUaccDerivativeDefinition(testCase)
            % Hand-computed derivative with the acceleration held constant.
            x    = [1; 2; 0.3; -0.4; 0.8];
            uacc = [0.6; 0.8];                 % ||u|| = 1
            dx   = ode_descent_uacc(x, uacc, testCase.Vc);
            expected = [0.3; -0.4; 0.6; 0.8 - 1; -testCase.Vc * 0.8 * 1];
            testCase.verifyEqual(dx, expected, 'AbsTol', 1e-15);
        end

        function testAccelerationIsDirect(testCase)
            % vx_dot = ux and vy_dot = uy - 1 are independent of the mass,
            % unlike ode_descent.m where the acceleration is T/m.
            uacc = [0.5; 1.2];
            dx1  = ode_descent_uacc([0;1;0;0;0.9], uacc, testCase.Vc);
            dx2  = ode_descent_uacc([0;1;0;0;0.3], uacc, testCase.Vc);
            testCase.verifyEqual(dx1(3:4), [0.5; 0.2], 'AbsTol', 1e-15);
            testCase.verifyEqual(dx2(3:4), [0.5; 0.2], 'AbsTol', 1e-15);
        end

        function testMassFlowScalesWithMass(testCase)
            % m_dot = -Vc*m*||u|| (i.e. d(ln m)/dt = -Vc*||u|| is constant).
            uacc = [0; 1];                     % ||u|| = 1
            dxA  = ode_descent_uacc([0;1;0;0;1.0], uacc, testCase.Vc);
            dxB  = ode_descent_uacc([0;1;0;0;0.5], uacc, testCase.Vc);
            testCase.verifyEqual(dxA(5), -testCase.Vc * 1.0, 'AbsTol', 1e-15);
            testCase.verifyEqual(dxB(5), -testCase.Vc * 0.5, 'AbsTol', 1e-15);
        end

        function testBallisticCoast(testCase)
            % Zero acceleration: free fall, no mass flow.
            dx = ode_descent_uacc([0.5;1;0.2;-0.1;0.9], [0;0], testCase.Vc);
            testCase.verifyEqual(dx, [0.2; -0.1; 0; -1; 0], 'AbsTol', 1e-15);
        end

        % ---- lti_zoh.m ----------------------------------------------------
        function testZohClosedForm(testCase)
            % Discrete matrices match the analytic double-integrator + integrator.
            h = testCase.dt;   vc = testCase.Vc;
            [Abar, Bbar, cbar] = lti_zoh(h, vc);
            A_exp = eye(5);  A_exp(1,3) = h;  A_exp(2,4) = h;
            B_exp = [h^2/2, 0,     0;
                     0,     h^2/2, 0;
                     h,     0,     0;
                     0,     h,     0;
                     0,     0,    -vc*h];
            c_exp = [0; -h^2/2; 0; -h; 0];
            testCase.verifyEqual(Abar, A_exp, 'AbsTol', 1e-12);
            testCase.verifyEqual(Bbar, B_exp, 'AbsTol', 1e-12);
            testCase.verifyEqual(cbar, c_exp, 'AbsTol', 1e-12);
        end

        function testZohMatchesOde45(testCase)
            % One ZOH step equals an ode45 integration of the LTI dynamics
            % with the control held constant over the interval.
            h = testCase.dt;   vc = testCase.Vc;
            [Abar, Bbar, cbar] = lti_zoh(h, vc);
            xi0 = [0.3; 1; -0.2; -0.6; 0];        % [x;y;vx;vy;z]
            w   = [0.4; 1.1; 1.3];                % [ux;uy;sigma]
            rhs = @(~, xi) [xi(3); xi(4); w(1); w(2)-1; -vc*w(3)];
            opts = odeset('RelTol', 1e-12, 'AbsTol', 1e-14);
            [~, XI] = ode45(rhs, [0, h], xi0, opts);
            testCase.verifyEqual(Abar*xi0 + Bbar*w + cbar, XI(end,:).', ...
                'AbsTol', 1e-9);
        end

        function testMassRowConsistency(testCase)
            % With the lossless cone active (sigma = ||u||), the LTI z-update
            % matches the exact nonlinear u-ZOH mass depletion.
            h = testCase.dt;   vc = testCase.Vc;
            uacc = [0.5; 1.0];   umag = norm(uacc);
            % Nonlinear replay of the mass channel (z = ln m is exactly linear).
            opts = odeset('RelTol', 1e-12, 'AbsTol', 1e-14);
            [~, X] = ode45(@(~,x) ode_descent_uacc(x, uacc, vc), [0 h], ...
                           [0;1;0;0;1], opts);
            z_nl  = log(X(end,5));               % nonlinear log-mass change
            z_lti = -vc * umag * h;              % LTI prediction with sigma=||u||
            testCase.verifyEqual(z_nl, z_lti, 'AbsTol', 1e-9);
        end
    end
end
