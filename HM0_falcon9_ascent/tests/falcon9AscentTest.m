classdef falcon9AscentTest < matlab.unittest.TestCase
    %falcon9AscentTest Physics and regression checks for the HM0 scripts.
    %  Runs main.m (dimensional) and main2.m (non-dimensional) once in the
    %  class setup, harvests their workspace results, and verifies the
    %  propellant bookkeeping, the event ordering and the agreement between
    %  the two independent implementations.
    %  Note: running the scripts regenerates the PNGs in figures/ (the
    %  repository pipeline owns those files, so this is intended).

    properties
        dim   % results harvested from main.m
        nd    % results harvested from main2.m
    end

    properties (TestParameter)
        impl = {'dim', 'nd'}
    end

    methods (TestClassSetup)
        function runScripts(testCase)
            hm0 = fileparts(fileparts(mfilename('fullpath')));
            testCase.dim = runAscentScript(fullfile(hm0, 'main.m'));
            testCase.nd  = runAscentScript(fullfile(hm0, 'main2.m'));
            testCase.addTeardown(@() close('all'));
        end
    end

    methods (Test)
        function testPropellantBookkeeping(testCase, impl)
            % dm/dt = -Qdot is constant: m(tb) = m0 - Qdot*tb exactly
            S = testCase.(impl);
            testCase.verifyEqual(S.mass(end), S.m0 - S.Qdot1*S.tb1, ...
                'RelTol', 1e-8);
        end

        function testMachOneBeforeMaxQ(testCase, impl)
            S = testCase.(impl);
            testCase.verifyNotEmpty(S.im1);
            testCase.verifyLessThan(S.t(S.im1), S.t(S.imQ));
        end

        function testAltitudeStaysAboveGround(testCase, impl)
            S = testCase.(impl);
            testCase.verifyGreaterThanOrEqual(min(S.h), -1e-6);
        end

        function testFinalAltitudeInPlausibleBand(testCase, impl)
            % Loose sanity band for the Falcon 9 first-stage MECO altitude
            S = testCase.(impl);
            testCase.verifyGreaterThan(S.h(end), 50e3);
            testCase.verifyLessThan(S.h(end), 120e3);
        end

        function testDimensionalVsNonDimensionalAgreement(testCase)
            % The two scripts integrate the same physics with different
            % scalings: results must agree to well below 1%
            testCase.verifyEqual(testCase.nd.h(end),    testCase.dim.h(end),    'RelTol', 5e-3);
            testCase.verifyEqual(testCase.nd.Vmag(end), testCase.dim.Vmag(end), 'RelTol', 5e-3);
            testCase.verifyEqual(testCase.nd.mass(end), testCase.dim.mass(end), 'RelTol', 1e-6);
            testCase.verifyEqual(testCase.nd.qmax,      testCase.dim.qmax,      'RelTol', 1e-2);
        end
    end
end

function S = runAscentScript(scriptPath)
%runAscentScript Run an HM0 main script and harvest the variables needed by
%  the tests. The script's leading `clear` wipes this workspace first, then
%  the script repopulates it; results are collected after the run.
    run(scriptPath);
    S = struct('t', t, 'h', h, 'mass', mass, 'Vmag', Vmag, ...
               'qdyn', qdyn, 'Mach', Mach, 'qmax', qmax, ...
               'm0', m0, 'Qdot1', Qdot1, 'tb1', tb1, ...
               'im1', im1, 'imQ', imQ);
end
