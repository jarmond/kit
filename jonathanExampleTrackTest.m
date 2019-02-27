function job = jonathanExampleTrackTest(job)
% Perform and test tracking
%
% take output from jonathanExampleDetectionTest.m and
% jonathanExampleRefinementTest.m and jonathanExampleFitTest
% rather than redoing
% Jonathan U Harrison 2019-02-12
%%%%%%%%%%

%%%%% can get output from metaphase example by rerunning previous tests
if nargin <1
    [spots, movie, job] = jonathanExampleDetectionTest();
    job = jonathanExampleRefinementTest(spots,movie,job);
    job = jonathanExamplePlaneFitTest(job);
end

channel = 1;
verbose = 1;
job.dataStruct{channel} = kitGenerateTracks(job.dataStruct{channel},verbose);
    
end


    