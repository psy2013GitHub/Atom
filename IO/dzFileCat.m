function outcell=dzFileCat(incell)

outcell={};
if iscell(incell)
    for ii=1:length(incell), outcell=[outcell;dzFileCat(incell{ii})]; end
elseif ischar(incell)
    outcell=[outcell,{incell}];           
end
return
end