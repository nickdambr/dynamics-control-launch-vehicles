classdef hm3LoopTest < matlab.unittest.TestCase
    %hm3LoopTest Unit tests for load_wind_profile, assemble_loop,
    %  design_controller and simulate_gust_response.
    %  Pins the closed-loop conclusions of the homework: the Task-1 PD design
    %  meets the |GM|/|PM| targets, the bare full model is bending-unstable,
    %  and the deep notch gain-stabilises it (Task 2).

    properties (Constant)
        % Task-1 PD design (auto-tuner output, pinned for regression)
        Kref = struct('Kp_th', 1.9800, 'Kd_th', 1.3997, ...
                      'Kp_z', -1e-3,   'Kd_z', -1e-3)
    end

    properties
        p
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

    methods (TestMethodSetup)
        function muteConditionallyStableMarginWarning(testCase)
            ws = warning('off', 'Control:analysis:MarginUnstable');
            testCase.addTeardown(@() warning(ws));
        end
    end

    methods (Test)
        function testGustProfileShape(testCase)
            % 1-cosine gust: zero before onset, peak Vg at mid-gust, alphaw = vw/V
            pp = testCase.p;
            w  = load_wind_profile(pp, 'Vg', 8.0, 'Tg', 3.0, 't0', 1.0);
            testCase.verifyEqual(max(abs(w.vw(w.t < 1.0))), 0, 'AbsTol', 1e-15);
            [vwPeak, iPeak] = max(w.vw);
            testCase.verifyEqual(vwPeak, 8.0, 'AbsTol', 1e-6);
            testCase.verifyEqual(w.t(iPeak), 1.0 + 1.5, 'AbsTol', w.t(2)-w.t(1));
            testCase.verifyEqual(w.alphaw, w.vw/pp.V, 'AbsTol', 1e-15);
        end

        function testStepProfileShape(testCase)
            w = load_wind_profile(testCase.p, 'profile', 'step', 'Vg', 5.0);
            testCase.verifyEqual(max(abs(w.vw(w.t < 1.0))), 0, 'AbsTol', 1e-15);
            testCase.verifyEqual(w.vw(end), 5.0, 'AbsTol', 1e-12);
        end

        function testDefaultGustAmplitudeFromDrywind(testCase)
            % Severe dry-wind dispersion at 15.1 km is the documented default
            w = load_wind_profile(testCase.p);
            testCase.verifyEqual(w.Vg, 6.38, 'AbsTol', 0.05);
        end

        function testWindProfileRejectsUnknownProfile(testCase)
            testCase.verifyError( ...
                @() load_wind_profile(testCase.p, 'profile', 'sinusoid'), ...
                'load_wind_profile:profile');
        end

        function testRigidLoopMeetsMarginTargets(testCase)
            % Pinned Task-1 gains reproduce |GM| ~ 6 dB, |PM| ~ 30 deg
            G = build_plant_rigid(testCase.p);
            [L, T] = assemble_loop(G, testCase.Kref);
            [Gm, Pm] = margin(L);
            testCase.verifyTrue(isstable(T));
            testCase.verifyEqual(abs(20*log10(Gm)), 6.0, 'AbsTol', 0.2);
            testCase.verifyEqual(abs(Pm), 30.0, 'AbsTol', 0.5);
        end

        function testDesignControllerMeetsTargets(testCase)
            G = build_plant_rigid(testCase.p);
            [K, m] = design_controller(G, [], 'verbose', false);
            testCase.verifyEqual(abs(m.GM_dB),  6.0,  'AbsTol', 0.1);
            testCase.verifyEqual(abs(m.PM_deg), 30.0, 'AbsTol', 0.5);
            testCase.verifyTrue(m.stable);
            testCase.verifyEqual(K.Kp_th, testCase.Kref.Kp_th, 'AbsTol', 5e-3);
            testCase.verifyEqual(K.Kd_th, testCase.Kref.Kd_th, 'AbsTol', 5e-3);
        end

        function testDesignControllerRestoresWarningState(testCase)
            % Regression: the tuner mutes Control:analysis:MarginUnstable
            % internally and must restore the caller's state on exit
            warning('on', 'Control:analysis:MarginUnstable');
            G = build_plant_rigid(testCase.p);
            design_controller(G, [], 'verbose', false);
            st = warning('query', 'Control:analysis:MarginUnstable');
            testCase.verifyEqual(st.state, 'on');
        end

        function testBareFullModelIsBendingUnstable(testCase)
            % Task 2, Step B: TVC + delay with no bending filter -> unstable
            pp = testCase.p;
            Gf = build_plant_full(pp, 'ins');
            [~, T] = assemble_loop(Gf, testCase.Kref, build_tvc(pp));
            testCase.verifyFalse(isstable(T));
        end

        function testDeepNotchStabilisesFullModel(testCase)
            % Task 2 retained design: deep notch gain-stabilises the resonance
            pp = testCase.p;
            Gf = build_plant_full(pp, 'ins');
            Hn = build_notch_filter(pp.wBM, 0.002, 0.7, +1);
            [L, T] = assemble_loop(Gf, testCase.Kref, build_tvc(pp)*Hn);
            testCase.verifyTrue(isstable(T));
            testCase.verifyLessThan(20*log10(abs(freqresp(L, pp.wBM))), -10);
        end

        function testGustResponseAngleOfAttackBudget(testCase)
            % alpha = theta + zdot/V + alpha_w, and peaks match the histories
            pp = testCase.p;
            G  = build_plant_rigid(pp);
            [~, T] = assemble_loop(G, testCase.Kref);
            w  = load_wind_profile(pp);
            r  = simulate_gust_response(T, w);
            testCase.verifyEqual(r.alpha, r.theta + r.zdot/pp.V + r.alphaw, ...
                'AbsTol', 1e-15);
            testCase.verifyEqual(r.peak_theta, max(abs(r.theta)), 'AbsTol', 1e-15);
            testCase.verifyEqual(r.theta(1), 0, 'AbsTol', 1e-15);
        end
    end
end
