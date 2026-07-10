classdef hm3LoopPerformanceTest < matlab.perftest.TestCase
    % Performance benchmarks for the HM3 hot path.
    %  Unit: one assemble_loop call (fminsearch hits it ~400x, the Task-2
    %  sweeps ~150x). System: full tuner cost (assemble_loop + margin +
    %  isstable) and the lsim gust replay. Workflow: one design_controller
    %  search. Run: runperf('hm3LoopPerformanceTest')

    properties (Constant)
        % Task-1 PD design (pinned), same as hm3LoopTest
        Kref = struct('Kp_th', 1.9800, 'Kd_th', 1.3997, ...
                      'Kp_z', -1e-3,   'Kd_z', -1e-3)
    end

    properties
        p
        Grigid
        Gfull
        Wchain      % TVC + delay + deep notch (Task-2 retained chain)
        Trigid      % closed loop for the gust replay
        wind
    end

    methods (TestClassSetup)
        function addHm3ToPath(testCase)
            hm3 = fileparts(fileparts(mfilename('fullpath')));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(hm3));
        end
    end

    methods (TestMethodSetup)
        function buildModels(testCase)
            % Model construction outside the measurement boundary; the
            % conditionally stable loop makes margin() warn on every call
            ws = warning('off', 'Control:analysis:MarginUnstable');
            testCase.addTeardown(@() warning(ws));
            testCase.p      = load_hw3_params();
            testCase.Grigid = build_plant_rigid(testCase.p);
            testCase.Gfull  = build_plant_full(testCase.p, 'ins');
            testCase.Wchain = build_tvc(testCase.p) * ...
                build_notch_filter(testCase.p.wBM, 0.002, 0.7, +1);
            [~, testCase.Trigid] = assemble_loop(testCase.Grigid, testCase.Kref);
            testCase.wind = load_wind_profile(testCase.p);
        end
    end

    methods (Test)
        function testAssembleLoopRigid(testCase)
            % fminsearch cost kernel: close the rigid loop (connect +
            % getLoopTransfer + minreal)
            G = testCase.Grigid;  K = testCase.Kref;
            while testCase.keepMeasuring
                [L, T] = assemble_loop(G, K);
            end
            testCase.verifyEqual(order(T), 4);
            testCase.verifyNotEmpty(L);
        end

        function testAssembleLoopFullChain(testCase)
            % Task-2 sweep kernel: 6-state plant + TVC + delay + notch
            G = testCase.Gfull;  K = testCase.Kref;  Wa = testCase.Wchain;
            while testCase.keepMeasuring
                [L, T] = assemble_loop(G, K, Wa);
            end
            testCase.verifyEqual(order(T), 13);   % 6 + 5 (TVC+Pade) + 2 (notch)
            testCase.verifyNotEmpty(L);
        end

        function testTunerCostEvaluation(testCase)
            % One full cost evaluation as design_controller performs it:
            % assemble_loop + margin + isstable on the rigid loop
            G = testCase.Grigid;  K = testCase.Kref;
            while testCase.keepMeasuring
                [L, T] = assemble_loop(G, K);
                [Gm, Pm] = margin(L);
                stable = isstable(T);
            end
            testCase.verifyTrue(stable);
            testCase.verifyTrue(isfinite(Gm) && isfinite(Pm));
        end

        function testGustResponseLsim(testCase)
            % Time-domain replay: lsim over the 2401-point severe gust
            T = testCase.Trigid;  w = testCase.wind;
            while testCase.keepMeasuring
                r = simulate_gust_response(T, w);
            end
            testCase.verifySize(r.theta, [numel(w.t) 1]);
        end

        function testDesignControllerSearch(testCase)
            % Workflow level: one complete PD margin-matching search
            G = testCase.Grigid;
            testCase.startMeasuring();
            K = design_controller(G, [], 'verbose', false);
            testCase.stopMeasuring();
            testCase.verifyEqual(K.Kp_th, testCase.Kref.Kp_th, 'AbsTol', 5e-3);
        end
    end
end
