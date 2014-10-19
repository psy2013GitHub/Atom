
function outcell=dzchar2cell(inchar,direc)

hflag=0;
if strcmpi(direc,'h') % horizontal
    hflag=1;
    len=size(inchar,2);
else                   % vertical
    len=size(inchar,1);
end

outcell=cell(len,1); 

for ii=1:len
    if(hflag)
        outcell{ii}=inchar(:,ii);
    else
        outcell{ii}=inchar(ii,:);
    end
end


return
end