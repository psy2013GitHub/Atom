
function emptycell_flag=delblankcell(incell)
if iscell(incell)
    if isempty(incell), emptycell_flag=0; return; end
    emptycell_flag=zeros(length(incell),1);
    for ii=1:length(incell)
        incell2=incell{ii};
        emptycell_flag(ii)=delblankcell(incell2);
    end
else
    if isempty(incell), 
        emptycell_flag=1;
    else
        emptycell_flag=0;
    end
end
end