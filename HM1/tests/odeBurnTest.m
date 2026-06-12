classdef odeBurnTest < matlab.unittest.TestCase
    %odeBurnTest Unit tests for ode_burn.m (HM1 powered-flight dynamics).
    %  Checks the linear-tangent steering law, the costate equations and
    %  two analytic limits (ballistic flight, vertical Tsiolkovsky burn).
    %  All quantities are non-dimensional (g = 1).

    methods (TestClassSetup)
        function addHm1ToPath(testCase)
            hm1 = fileparts(fileparts(mfilename('fullpath')));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(hm1));
        end
    end

    methods (Test)
        function testKinematicsAndMassFlow(testCase)
            % dx = vx, dy = vy, dm = -Q regardless of the costates
            p = struct('T', 1.2, 'Q', 2, 'c', 0.6, ...
                       'lam_vx0', 1, 'lam_vy0', 0.5, 'lam_y', 0.2);
            z = [0.1; 0.2; 0.3; 0.4; 0.9; 1.1];
            dz = ode_burn(0.5, z, p);
            testCase.verifySize(dz, [6 1]);
            testCase.verifyEqual(dz(1), z(3), 'AbsTol', 1e-15);
            testCase.verifyEqual(dz(2), z(4), 'AbsTol', 1e-15);
            testCase.verifyEqual(dz(5), -p.Q, 'AbsTol', 1e-15);
        end

        function testThrustAlongConstantCostate(testCase)
            % lam = (1, 0) constant -> phi = 0: horizontal thrust, dvy = -g
            p = struct('T', 1.2, 'Q', 2, 'c', 0.6, ...
                       'lam_vx0', 1, 'lam_vy0', 0, 'lam_y', 0);
            m = 0.8;
            z = [0; 0; 0.5; 0.1; m; 1];
            dz = ode_burn(0.3, z, p);
            testCase.verifyEqual(dz(3), p.T/m, 'AbsTol', 1e-14);
            testCase.verifyEqual(dz(4), -1, 'AbsTol', 1e-14);
        end

        function testLinearTangentSwitch(testCase)
            % lam_vy(t) = lam_vy0 - lam_y*t vanishes at t* = lam_vy0/lam_y,
            % where the thrust must be exactly horizontal (phi = 0)
            p = struct('T', 1.2, 'Q', 2, 'c', 0.6, ...
                       'lam_vx0', 1, 'lam_vy0', 2, 'lam_y', 4);
            tStar = p.lam_vy0 / p.lam_y;
            m = 1;
            z = [0; 0; 0; 0; m; 1];
            dz = ode_burn(tStar, z, p);
            testCase.verifyEqual(dz(3), p.T/m, 'AbsTol', 1e-14);
            testCase.verifyEqual(dz(4), -1, 'AbsTol', 1e-14);
        end

        function testCostateMassEquation(testCase)
            % dlam_m/dt = (T/m^2) * |lam_v|, with |lam_v| = 5 for lam = (3,4)
            p = struct('T', 1.2, 'Q', 2, 'c', 0.6, ...
                       'lam_vx0', 3, 'lam_vy0', 4, 'lam_y', 0);
            m = 0.5;
            z = [0; 0; 0; 0; m; 1];
            dz = ode_burn(0, z, p);
            testCase.verifyEqual(dz(6), p.T/m^2 * 5, 'AbsTol', 1e-12);
        end

        function testBallisticLimitMatchesAnalytic(testCase)
            % T = Q = 0: pure ballistic flight under unit gravity
            p = struct('T', 0, 'Q', 0, 'c', 0.6, ...
                       'lam_vx0', 1, 'lam_vy0', 1, 'lam_y', 0);
            z0 = [0; 0; 0.3; 0.5; 1; 1];
            tf = 0.8;
            opts = odeset('RelTol', 1e-10, 'AbsTol', 1e-12);
            [~, Z] = ode45(@(t,z) ode_burn(t, z, p), [0 tf], z0, opts);
            testCase.verifyEqual(Z(end,1), z0(3)*tf,              'AbsTol', 1e-9);
            testCase.verifyEqual(Z(end,2), z0(4)*tf - 0.5*tf^2,   'AbsTol', 1e-9);
            testCase.verifyEqual(Z(end,3), z0(3),                 'AbsTol', 1e-10);
            testCase.verifyEqual(Z(end,4), z0(4) - tf,            'AbsTol', 1e-9);
            testCase.verifyEqual(Z(end,5), z0(5),                 'AbsTol', 1e-12);
        end

        function testVerticalBurnTsiolkovskyWithGravity(testCase)
            % lam = (0, 1) constant -> phi = 90 deg: vertical burn.
            % Analytic: vy(t) = (T/Q)*ln(m0/(m0 - Q*t)) - t  with m0 = 1
            p = struct('T', 1.2, 'Q', 2, 'c', 0.6, ...
                       'lam_vx0', 0, 'lam_vy0', 1, 'lam_y', 0);
            z0 = [0; 0; 0; 0; 1; 1];
            tf = 0.3;
            opts = odeset('RelTol', 1e-12, 'AbsTol', 1e-14);
            [~, Z] = ode45(@(t,z) ode_burn(t, z, p), [0 tf], z0, opts);
            vyAnalytic = (p.T/p.Q) * log(1/(1 - p.Q*tf)) - tf;
            testCase.verifyEqual(Z(end,4), vyAnalytic, 'AbsTol', 1e-9);
            testCase.verifyEqual(Z(end,3), 0,          'AbsTol', 1e-12);
            testCase.verifyEqual(Z(end,5), 1 - p.Q*tf, 'AbsTol', 1e-12);
        end
    end
end
