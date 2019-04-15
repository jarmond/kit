function T = kitMakeAnalysisTable(job,channel)
%%
%%Take a processed, paired jobset file and convert analysis to a convenient
%%table format
%% Will account for filtered spots and ignore any that have been filtered out
%%Jonathan U Harrison 2019-04-12
%%%%%%%%%%%%%

if nargin <2
    channel=1;
end

if ~isfield(job,'dataStruct')
    error('Expect a jobset file that has already been processed and paired');
end

dataStruct = job.dataStruct{channel};

nFrames = job.metadata.nFrames;
if nFrames > 1
    warning('Categories cannot yet be used for movies. Making a basic table')
    if isfield(dataStruct,'sisterList')         
        nSisters = length(dataStruct.sisterList);
        %create columns for table
        position = zeros(nFrames*nSisters*2,3);
        amplitude = zeros(nFrames*nSisters*2,3);
        frame = repmat(repmat((1:nFrames)',nSisters,1),2,1);
        sisterPairID = repmat(1:nSisters,nFrames,1);
        sisterPairID = repmat(sisterPairID(:),2,1);
        sisterID = [ones(nFrames*nSisters,1);2*ones(nFrames*nSisters,1)];
        for k = [1,2]
            for j = 1:nSisters
                for i = 1:nFrames
                    if k==1 %use coords1 or coords2
                        position((k-1)*nFrames*nSisters + (j-1)*nFrames + i,:) = ...
                            dataStruct.sisterList(j).coords1(i,1:3);
                        amplitude((k-1)*nFrames*nSisters + (j-1)*nFrames + i,:) = ...
                            dataStruct.sisterList(j).coords1(i,4:6);
                    elseif k==2
                        position((k-1)*nFrames*nSisters + (j-1)*nFrames + i,:) = ...
                            dataStruct.sisterList(j).coords2(i,1:3);
                        amplitude((k-1)*nFrames*nSisters + (j-1)*nFrames + i,:) = ...
                            dataStruct.sisterList(j).coords2(i,4:6); 
                    end
                end
            end
        end

        varnames = {'Position','Amplitude','Frame','SisterPairID','SisterID'};
        T = table(position,amplitude,frame,sisterPairID,sisterID,...
            'VariableNames', varnames);
    else
        error('No sisterList found but need paired joblist to make table from');
    end
else
    %put all spots into a table structure
    T = table(dataStruct.initCoord.allCoord(:,1:3), ...
        dataStruct.initCoord.allCoord(:,4:6), ...
        'VariableNames',{'Position','Amplitude'});
    % add category information from manual labelling via gui
    if isfield(job,'categories') %have added some kind of categories
        categoryNames = fieldnames(job.categories);
        for j=1:numel(fieldnames(job.categories))
            T.(categoryNames{1}) = zeros(size(T,1),1);
            T.(categoryNames{1})(job.categories.(categoryNames{1})) = 1;
        end
    end
end

