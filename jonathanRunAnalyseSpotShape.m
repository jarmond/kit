function T = jonathanRunAnalyseSpotShape(jobset,saveToCSV,channel,verbose)
if nargin<2
    saveToCSV=0;
end
if nargin < 3
    channel = 1;
end
if nargin<4
    verbose=0;
end

%get the data
job = kitLoadAllJobs(jobset);
nMovs = length(job);

for jobInd = 1:nMovs
        fprintf('Processing the %dth job...\n',jobInd);
        if ~isfield(job{jobInd},'dataStruct')
            error('Expect a jobset file that has already been processed and paired');
        end
        dataStruct = job{jobInd}.dataStruct{channel};
        if ~isfield(dataStruct,'failed') || dataStruct.failed
            warning('Movie %d has failed so skipping \n',jobInd);
            continue
        end
        stt = jonathanAnalyseSpotShape(job{jobInd},'verbose',verbose);
        
        % if verbose
%             dt = 4.7;
%             pairID = 35;
%         tbl = kitMakeAnalysisTable(job201202,0);
%         figure;
%         title('Spot scale (pxls) over time')
%         subplot(4,1,1);
%         ylabel('x (pxls)')
%         subplot(4,1,2);
%         ylabel('y (pxls)')
%         subplot(4,1,3);
%         ylabel('z (pxls)')
%         time_vec = ((1:job201202.metadata.nFrames)-1)*dt;
%         for k=1:2
%             for j=1:3
%                 subplot(4,1,j);
%                 hold all;
%                 plot(time_vec,stt(:,j,k,pairID),'linewidth',2)
%                 axis([0,550,-Inf,Inf]);
%                 xlabel('Time (s)')
%                 set(gca,'fontsize',20)
%             end
%             subplot(4,1,4)
%             hold all;
%             plot(time_vec,tbl((tbl.SisterPairID==pairID)&(tbl.SisterID==k),:).Position(:,1),'c','linewidth',2);
%             axis([0,550,-Inf,Inf]);
%             xlabel('Time (s)')
%             ylabel('Position (um)')
%             set(gca,'fontsize',20)
%         end
        % end
        
        
        %%%%%%%%%%%%%%%%%%%%
        nFrames = size(stt,1);
        nSisters = size(stt,4);
        frame = repmat(repmat((1:nFrames)',nSisters,1),2,1);
        sisterPairID = repmat(1:nSisters,nFrames,1);
        sisterPairID = repmat(sisterPairID(:),2,1);
        sisterID = [ones(nFrames*nSisters,1);2*ones(nFrames*nSisters,1)];
        stretch = zeros(nFrames*nSisters*2,3);
        for k = [1,2]
            for j = 1:nSisters
                for i = 1:nFrames
                    stretch((k-1)*nFrames*nSisters + (j-1)*nFrames + i,:) = ...
                        stt(i,:,k,j);
                end
            end
        end
        varnames = {'Stretch','Frame','SisterPairID','SisterID'};
        T = table(stretch,frame,sisterPairID,sisterID,...
            'VariableNames', varnames);
        %tbl = kitMakeAnalysisTable(job{jobInd},0);
        outname = kitGenerateOutputFilename(job{jobInd});
        outname = strcat(outname(1:(end-3)),'csv');
        if exist(outname)
            tbl = readtable(outname);
            T = join(tbl,T);
        else
            fprintf('Could not find csv file of tracks to join spot shape info to. Has the csv file been generated?\n');
        end
        
        if saveToCSV
            outname = kitGenerateOutputFilename(job{jobInd});
            outname = strcat(outname(1:(end-4)),'_spot_stretch.csv'); %replace file ending of mat with csv
            writetable(T, outname);
        end
%     catch
%         fprintf('Unsuccessful in converting to table: job %d run correctly?\n',jobInd);
%         continue
%     end
end
