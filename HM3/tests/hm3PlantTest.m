classdef hm3PlantTest < matlab.unittest.TestCase
    % Unit tests for load_hw3_params / build_plant_rigid / build_plant_full.
    %  Pins the Greensite pitch-plane physics at max-qbar (t = 72 s): Table-1
    %  coefficients, the unstable airframe pole at +sqrt(A6), and the INS
    %  bending contamination of Eq. (2).

    properties
        p    % nominal parameter struct, loaded once per class
    end

    methods (TestClassSetup)
        function addHm3ToPath(testCase)
            hm3 = fileparts(fileparts(mfilename('fullpath')));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(hm3));
        end

        function loadNominalParams(testCase)
            testCase.p = load_hw3_params();
        end
    end

    methods (Test)
        function testParamsMatchTable1(testCase)
            % LPV data interpolated at t = 72 s must agree with Table 1
            pp = testCase.p;
            testCase.verifyEqual(pp.A6,  3.3818, 'AbsTol', 5e-3);
            testCase.verifyEqual(pp.K1,  4.5647, 'AbsTol', 5e-3);
            testCase.verifyEqual(pp.V,   937.70, 'AbsTol', 0.5);
            testCase.verifyEqual(pp.wBM, 18.9,   'AbsTol', 0.05);
            testCase.verifyEqual(pp.a4,  -27.2710, 'AbsTol', 5e-2);
        end

        function testDynamicPressureSelfConsistent(testCase)
            % qbar = 0.5 * rho0 * exp(-Alt/Hscale) * V^2 (exponential atmosphere)
            pp = testCase.p;
            qbarExpected = 0.5 * 1.225 * exp(-pp.Alt/8000) * pp.V^2;
            testCase.verifyEqual(pp.qbar, qbarExpected, 'AbsTol', 1e-6);
        end

        function testUncertaintyScalingAppliesToA6K1(testCase)
            % Task-3 corner scaling multiplies mu_alpha = A6 and mu_c = K1 only
            pp = testCase.p;
            ps = load_hw3_params('mu_alpha_scale', 1.3, 'mu_c_scale', 0.7);
            testCase.verifyEqual(ps.A6, 1.3*pp.A6, 'AbsTol', 1e-12);
            testCase.verifyEqual(ps.K1, 0.7*pp.K1, 'AbsTol', 1e-12);
            testCase.verifyEqual(ps.a3, pp.a3,     'AbsTol', 1e-12);
        end

        function testRigidPlantDimensionsAndNames(testCase)
            G = build_plant_rigid(testCase.p);
            testCase.verifySize(G.A, [4 4]);
            testCase.verifySize(G.B, [4 2]);
            testCase.verifyEqual(G.InputName,  {'delta'; 'alpha_w'});
            testCase.verifyEqual(G.OutputName(1), {'theta_m'});
        end

        function testRigidAirframeUnstablePole(testCase)
            % Aerodynamically unstable airframe: dominant pole at ~ +sqrt(A6)
            % (the a1/a4 drift coupling shifts it ~1% off the pitch-only value)
            G = build_plant_rigid(testCase.p);
            testCase.verifyEqual(max(real(pole(G))), sqrt(testCase.p.A6), ...
                'RelTol', 0.02);
        end

        function testRigidMeasurementsEqualTrueStates(testCase)
            % No bending mode -> INS rows coincide with the true-state rows
            G = build_plant_rigid(testCase.p);
            testCase.verifyEqual(G.C(1,:), G.C(5,:), 'AbsTol', 1e-15);  % theta
            testCase.verifyEqual(G.C(3,:), G.C(6,:), 'AbsTol', 1e-15);  % z
        end

        function testFullPlantBendingMode(testCase)
            % States 5-6 carry the bending oscillator (wBM, zBM) exactly
            pp = testCase.p;
            G  = build_plant_full(pp);
            [wn, zeta] = damp(G);
            [~, k] = min(abs(wn - pp.wBM));
            testCase.verifyEqual(wn(k),   pp.wBM, 'AbsTol', 1e-9);
            testCase.verifyEqual(zeta(k), pp.zBM, 'AbsTol', 1e-9);
        end

        function testInsMeasurementsContaminatedByBending(testCase)
            % Eq. (2): eta leaks into theta_m (+sigma_ins) and z_m (-phi_ins)
            pp = testCase.p;
            G  = build_plant_full(pp, 'ins');
            testCase.verifyEqual(G.C(1,5),  pp.sigma_ins, 'AbsTol', 1e-15);
            testCase.verifyEqual(G.C(3,5), -pp.phi_ins,   'AbsTol', 1e-15);
        end

        function testTrueMeasurementsBypassBending(testCase)
            % 'true' feedback: zero bending columns in the measurement block
            G = build_plant_full(testCase.p, 'true');
            testCase.verifyEqual(G.C(1:4, 5:6), zeros(4, 2), 'AbsTol', 1e-15);
        end

        function testFullPlantRejectsUnknownMeas(testCase)
            testCase.verifyError(@() build_plant_full(testCase.p, 'bogus'), ...
                'build_plant_full:meas');
        end
    end
end
