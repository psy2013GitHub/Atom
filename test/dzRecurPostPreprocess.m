

function outcell=dzRecurPostPreprocess(incell,Options)

for ii=1:length(incell)
    tmp=incell{ii};
    if ~ischar(tmp{1})
        outcell{ii}=dzRecurPostPreprocess(tmp);
    else
        inFiles=tmp;
        for jj=1:length(Options.cmd)
            eval(Options.cmd);
        end
        outcell=Options.fname;
    end
end
end


