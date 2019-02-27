function job = jonathanExampleGroupSistersTest(job,verbose,movie,theta)
% Perform and test tracking
%
% take output from jonathanExampleDetectionTest.m and
% jonathanExampleRefinementTest.m and jonathanExampleFitTest.m and
% jonathanExampleGroupSistersTest.m
% rather than redoing
%
% movie input optional for visualisation
% Jonathan U Harrison 2019-02-12
%%%%%%%%%%

channel = 1;

%%%%% can get output from metaphase example by rerunning previous tests
if nargin <1
    [spots, movie, job] = jonathanExampleDetectionTest();
    job = jonathanExampleRefinementTest(spots,movie,job);
    job = jonathanExamplePlaneFitTest(job);
    job = jonathanExampleTrackTest(job);
end

if nargin<2 || isempty(verbose)
verbose = 1;
end

if nargin<4
    %set default options for grouping sisters
    opts.useAlignment = 1;
    opts.maxAngle = 30;
    opts.maxDist = 1.5;
    opts.minOverlap = 10;
    opts.useAnaphase=0;
    opts.robust=0;
else %use input, can set up as optimization problem
    opts.useAlignment = 1; 
    opts.maxAngle = theta(1);
    opts.maxDist = theta(2);
    opts.minOverlap = theta(3);
    opts.useAnaphase=1;
    opts.robust=0;
end


job.dataStruct{channel}.dataProperties.groupSisters=opts;

job.dataStruct{channel} = kitGroupSisters(job.dataStruct{channel},...
    verbose);

%how many sisters does this give us on average through the movie?
    nFrames = job.metadata.nFrames;
    nSisters = size(job.dataStruct{channel}.sisterList,1);
    doesSisterPairExist = zeros(nFrames,nSisters);
    for j = 1:nSisters
        %put ones where a paired sister exists
        doesSisterPairExist(:,j)= ~isnan(...
            job.dataStruct{channel}.sisterList(j).coords1(:,1));
    end
    %multiply by 2 since we count the pair twice
    nSistersMovieAverage = 2*mean(sum(doesSisterPairExist,2));
    fprintf('On average through the movie we find %f sisters\n',...
        nSistersMovieAverage);

%can we plot which sisters we have matched up
if verbose && ~isempty(movie)
    frameToShow = min(50,nFrames); %piack a frame to plot sisters on
    jonathanVisualiseTrackedAndGroupedSpots(job,movie,channel,frameToShow);
elseif verbose
    %load movie if input movie is empty. eg:
    [~, movie, ~] = jonathanExampleDetectionTest(job);
    frameToShow = min(50,nFrames); %piack a frame to plot sisters on
    jonathanVisualiseTrackedAndGroupedSpots(job,movie,channel,frameToShow);
end