

function Atom=Atom_init(Atom)

% This Function Is Used For Parameter Initialization.

Atom_defaults;
global DEFAULTS;

%- TR
try, TR=Atom.TR;
catch,
    if isfield(Atom,'SliceTiming')&&isfield(Atom.SliceTiming,'TR')&&~isnumeric(Atom.SliceTiming.TR)
        Atom=Atom.SliceTiming.TR;
    elseif 1 % fix me
    end
end



Atom=recurCpField(Atom,DEFAULTS);

return;
end

function Atom=recurCpField(Atom,Defaults)
if isstruct(Atom)&&isstruct(Defaults)
    fdnames=fieldnames(Defaults);
    for ii=1:length(fdnames)
        if isfield(Atom,fdnames{ii})
            Atom.(fdnames{ii})=recurCpField(Atom.(fdnames{ii}),Defaults.(fdnames{ii}));
        else
            Atom.(fdnames{ii})=Defaults.(fdnames{ii});
        end
    end
end

return
end