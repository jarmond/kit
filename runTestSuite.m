% run a suite of tests to check basic functionality of KiT
% for tracking of Kinetochores in spinning disk or LLSM data
%
% Jonathan U Harrison 2019-02-12
%%%%%%%%%%%%%%%%%%%%%%%%%%

%see here for details on unit tests and tags
%actually run via results = runtests({'runTestSuite'},'Tag','Diagnose');

classdef runTestSuite < matlab.unittest.TestCase
    methods (Test, TestTags = {'Data'})
        function loadExampleData (testCase)
            jonathanDataTest();
        end
    end
    methods (Test, TestTags = {'Detection'})
        function detection (testCase)
            jonathanExampleDetectionTest();
        end
        function refinement (testCase)
            [spots, movie, job] = jonathanExampleDetectionTest();
            jonathanExampleRefinementTest(spots,movie,job);
        end
    end
    methods (Test, TestTags = {'Plane'})
        function planeFit (testCase)
            jonathanExamplePlaneFitTest();
        end
    end
    methods (Test, TestTags = {'Tracking'})
        function tracking (testCase)
            jonathanExampleTrackTest(); 
        end
    end
    methods (Test, TestTags = {'Sisters'})
        function groupSisters (testCase)
            jonathanExampleGroupSistersTest();
        end
    end
    methods (Test, TestTags = {'Diagnose'})
        function diagnose (testCase)
            jonathanDiagnoseAndTest();
        end
    end
end