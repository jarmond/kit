function lazyKTs = kitIdentifyLazyKTs(job,channel)
%provide indices of tracks that are suggested to be lazy or lagging
%kinetochores
%
%Jonathan U. Harrison 2019-06-19
%%%%%%%%%%%%%%%%
if nargin<2
    channel=1;
end
if isfield(job,'dataStruct')
    dataStruct = job.dataStruct{channel};
else
    error('No dataStruct present in job');
end
lazyKTs=[];
nFrames = job.metadata.nFrames;
for iFrame=1:nFrames
    laggingIdx = dataStruct.planeFit(iFrame).laggingIdx; %indexes into initcoord
    for j=laggingIdx
        for k=1:length(dataStruct.trackList)
            if sum((dataStruct.trackList(k).coords(iFrame,1:3) - dataStruct.planeFit(iFrame).planeCoord(j,1:3)).^2)<eps             
                lazyKTs = [lazyKTs, k];
            end
        end
    end
end
end
