

function Options=dzRecurDir_options(oper)


switch lower(oper)
    case 'subjmatch'
        Options.fun.names={'subjmatch'};
        Options.fun.methods={'everydirect'};
        Options.cmd{1}=['fname=inFiles{1}(1,:);',...
                        'subjidx=cellfun(@cat,repmat({2},size(Options.SubjID)),repmat({filesep},size(Options.SubjID)),Options.SubjID,''UniformOutput'',false);',... % /subj...
                        'subjidx=cellfun(@cat,repmat({2},size(Options.SubjID)),subjidx,repmat({filesep},size(Options.SubjID)),''UniformOutput'',false);',... % /subj../
                        'subjidx=cellfun(@strfind,repmat({fname},size(subjidx)),subjidx,''UniformOutput'',false);',...
                        'subjidx=[subjidx{:}];',...
                        'if length(subjidx)~=1, error(''file name error in subj match in post-preprocess''); end;',...
                        'subjidx1=strfind(fname(subjidx+1:end),filesep);',...
                        'Options.subj=fname(subjidx+1:subjidx+subjidx1-1);',...  % 1st subj
                        'Options.subjidx=strmatch(Options.subj,Options.SubjID);' % 2nd subjidx
                       ]; % specific arg
    case 'moveeveryfile'
        Options.fun.names={'movefile'}; 
        Options.fun.methods={'everyfile'};
        Options.cmd{1}=['tmpfile=inFile;',...
                        '[tmpp,tmpf,tmpe]=fileparts(tmpfile);',...
                        'tmpidx=strfind(tmpp,Options.RawDir);',... % Options.RawDir
                        'tmpoutputdir=[Options.DestDir,filesep,tmpp(tmpidx(1)+length(Options.RawDir)+1:end)];',... % Options.DestDir
                        'if exist(tmpoutputdir,''dir'')~=7, mkdir(tmpoutputdir); end; ',... % mkdir, if exist, then ignore
                        'movefile(tmpfile,tmpoutputdir,''f'');',...
                        'inFile=[tmpoutputdir,filesep,tmpf,tmpe];'
                       ]; % specific arg: 'inFiles', see below & only the first file will be converted
    otherwise
        error('No Such Operation');
end


end