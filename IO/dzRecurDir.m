

function [FileList__,Options]=dzRecurDir(Dir,Options)

% default in spm way, i.e. each file in a cell
% if specified in consice way, each directory is a cell


try, tmp=Options.fileFilter;       catch, Options.fileFilter={'*'};  end
try, tmp=Options.dirFilter;        catch, Options.dirFilter={'*'};   end
try, tmp=Options.method;           catch, Options.method='spm';      end
try, tmp=Options.fun.names;        catch, Options.fun.names={''};    end
% everydirec/everyfile/firstfile
if isdir(Dir)
    % filter files__;
    files__=[]; firstfileDoneFlag=0;
    for fF=1:length(Options.fileFilter)
        tmpfiles=dir([Dir,filesep,Options.fileFilter{fF}]);
        tmpfiles=tmpfiles(~logical([tmpfiles(:).isdir]));
        files__=[files__;tmpfiles];
    end
    nFiles__=length(files__);
    % filter directorys
    dirs__=[];
    for dF=1:length(Options.dirFilter)
        tmpdirs =dir([Dir,filesep,Options.dirFilter{dF}]); tmpdirs=tmpdirs(3:end);
        tmpdirs=tmpdirs(logical([tmpdirs(:).isdir]));
        dirs__=[dirs__;tmpdirs];
    end
    nDirs__=length(dirs__);
    %     nDirs__ =length(dirs__);
    %     if length(dirs__)==1&&strcmpi(dirs__(1).name,'.')||strcmpi(dirs__(1).name,'..'), nDirs__=0; end
    %     if length(dirs__)==2&&strcmpi(dirs__(1).name,'.')||strcmpi(dirs__(2).name,'..'), nDirs__=0; end
    newDir=[files__;dirs__]; % files__ on the top
    fileFlag__=0;
    idxCell__=1; % file
    FileList__={}; % in case of empty output
    for ii__=1:length(newDir)
        idx=ii__;
        if newDir(ii__).isdir
            if idxCell__>1&&ischar(FileList__(1)), FileList__={FileList__}; end % in case FileList__ is char and FileList__{2}=* expression is not valid
            FileList__{idxCell__,1}=dzRecurDir([Dir,filesep,newDir(ii__).name],Options);FileList__{idxCell__,1}
            if ~isempty(FileList__{idxCell__,1})
                % ------------------------------  removeble  ------------------------------ %
                if strcmpi(Options.method,'concise')&&ischar(FileList__{idxCell__,1})...
                        %                         &&length(FileList__{idxCell__,1})==1 % 'concise' & 'char' & single directory
                    if ~isempty([Options.fun.names{:}])
                        inFiles=FileList__(idxCell__,1);
                        for k__=1:length(Options.fun.names)
                            if strcmpi(Options.fun.methods{k__},'everyfile'), continue; end % fix me
                            if strcmpi(Options.fun.methods{k__},'firstfile'), continue; end
                            %                             if strcmpi(Options.fun.methods{k__},'inherit'),   continue; end
                            eval(Options.cmd{k__});
                        end
                        FileList__{idxCell__,1}=inFiles; % remember for 'concise' mode, u need to align with 'spm' mode
                    end
                end
                % ------------------------------     end     ------------------------------ %
                idxCell__=idxCell__+1;
            else
                FileList__(idxCell__)=[]; % delete empty ones
            end
        else
            if firstfileDoneFlag, continue; end % if 'firsfile' operation already done, then go on
            File__=[Dir,filesep,newDir(ii__).name];
            % ------------------------------  removeble  ------------------------------ %
            if ~isempty([Options.fun.names{:}])
                inFile=File__;
                for k__=1:length(Options.fun.names)
                    % specificity 'Dcm2Nii'_____________________________
                    if strcmpi(Options.fun.methods{k__},'firstfile')&&ii__==1
                        eval(Options.cmd{k__});
                        firstfileDoneFlag=1;
                        continue;
                    end
                    %                      _____________________________
                    if ~strcmpi(Options.fun.methods{k__},'everyfile'), continue; end
                    eval(Options.cmd{k__});
                end
                File__=inFile; % remember for 'concise' mode, u need to align with 'spm' mode
            end
            % ------------------------------     end    ------------------------------ %
            if strcmpi(Options.method,'concise')   % in a concise way
                if ~fileFlag__ % first file flag
                    ii_Start=idx;
                    idxCell__=ii_Start+1;
                    if ~nDirs__ % no dorectory
                        FileList__=File__;
                    else
                        FileList__{ii_Start,1}=File__;
                    end
                else
                    if ~nDirs__ % no dorectory
                        FileList__=strvcat(FileList__,File__);
                    else
                        FileList__{ii_Start,1}=strvcat(FileList__{ii_Start},File__);
                    end
                end
                fileFlag__=fileFlag__+1;
            elseif strcmpi(Options.method,'spm')     % in spm way, each file in a cell
                FileList__{idxCell__,1}=[File__];
                idxCell__=idxCell__+1;
            else
                error('Unknown method type: %f',Options.method);
            end
        end
    end
else
    error('%s not a directory',Dir);
end

% ------------------------------  removeble  ------------------------------ %
% subjmatch & post-preprocess
everydirectDoneFlag__=0;
if ~nDirs__&&strcmpi(Options.method,'concise')
    for k__=1:length(Options.fun.names) % fix me
        if strcmpi(Options.fun.methods{k__},'everydirect')
            fprintf('!!!warning: Since only files exist, methods will be converted from ''everydirect'' ''everyfile''');
            everydirectDoneFlag__=1;
        elseif  everydirectDoneFlag__&&strcmpi(Options.fun.methods{k__},'inherit')
            
        else
            continue;
        end
        inFiles={FileList__};
        eval(Options.cmd{k__});
        FileList__=inFiles; % return as cell
    end
end
% ------------------------------  removeble  ------------------------------ %

end