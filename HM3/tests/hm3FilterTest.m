classdef hm3FilterTest < matlab.unittest.TestCase
    % Unit tests for build_tvc and build_notch_filter.
    %  Pins Eq.-3 actuator (unity DC gain, 2nd-order + Pade order) and Eq.-4
    %  section (depth zN/zD at centre, unity gain far off, RHP zeros for NMP).

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

    methods (Test)
        function testTvcUnityDcGain(testCase)
            % Actuator and Pade sections both have unit static gain
            Wtvc = build_tvc(testCase.p);
            testCase.verifyEqual(dcgain(Wtvc), 1, 'AbsTol', 1e-9);
        end

        function testTvcOrderIsActuatorPlusPade(testCase)
            % 2nd-order actuator + nth-order Pade delay
            testCase.verifyEqual(order(build_tvc(testCase.p)),    5);  % default n = 3
            testCase.verifyEqual(order(build_tvc(testCase.p, 2)), 4);
        end

        function testTvcDelayPhaseAtLowFrequency(testCase)
            % Pade reproduces exp(-j*w*tau) phase lag well below its bandwidth;
            % total phase = actuator + delay
            pp = testCase.p;
            wTest = 2.0;                                   % rad/s, low frequency
            h  = freqresp(build_tvc(pp), wTest);
            hA = freqresp(tf(pp.wTVC^2, [1, 2*pp.zTVC*pp.wTVC, pp.wTVC^2]), wTest);
            phaseDelay = angle(h) - angle(hA);
            testCase.verifyEqual(phaseDelay, -wTest*pp.tau, 'AbsTol', 1e-6);
        end

        function testNotchUnityGainFarFromCentre(testCase)
            wx = testCase.p.wBM;
            Hx = build_notch_filter(wx, 0.002, 0.7, +1);
            testCase.verifyEqual(dcgain(Hx), 1, 'AbsTol', 1e-9);
            testCase.verifyEqual(abs(freqresp(Hx, 1e4)), 1, 'AbsTol', 1e-3);
        end

        function testNotchDepthIsZetaRatio(testCase)
            % |Hx(j*wx)| = zN/zD for both numerator-sign variants
            wx = testCase.p.wBM;  zN = 0.002;  zD = 0.7;
            Hn  = build_notch_filter(wx, zN, zD, +1);
            Hll = build_notch_filter(wx, zN, zD, -1);
            testCase.verifyEqual(abs(freqresp(Hn,  wx)), zN/zD, 'AbsTol', 1e-9);
            testCase.verifyEqual(abs(freqresp(Hll, wx)), zN/zD, 'AbsTol', 1e-9);
        end

        function testDefaultVariantIsNonMinimumPhase(testCase)
            % sgn = -1 (assignment Eq. 4 as printed): both zeros in the RHP
            Hx = build_notch_filter(testCase.p.wBM, 0.2, 0.5);
            testCase.verifyTrue(all(real(zero(Hx)) > 0));
        end

        function testMinimumPhaseVariantHasLhpZeros(testCase)
            Hx = build_notch_filter(testCase.p.wBM, 0.2, 0.5, +1);
            testCase.verifyTrue(all(real(zero(Hx)) < 0));
        end

        function testNotchRejectsNegativeFrequency(testCase)
            testCase.verifyError(@() build_notch_filter(-1, 0.2, 0.5), ...
                'MATLAB:validators:mustBePositive');
        end

        function testNotchRejectsInvalidSign(testCase)
            testCase.verifyError(@() build_notch_filter(18.9, 0.2, 0.5, 0), ...
                'MATLAB:validators:mustBeMember');
        end

        function testTvcRejectsNonIntegerPadeOrder(testCase)
            testCase.verifyError(@() build_tvc(testCase.p, 2.5), ...
                'MATLAB:validators:mustBeInteger');
        end
    end
end
