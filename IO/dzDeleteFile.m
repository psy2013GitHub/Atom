



function dzDeleteFile(incell)

if iscell(incell)
   for ii=1:length(incell), dzDeleteFile(incell{ii}); end
elseif ischar(incell)
    delete(incell);
end


return
end