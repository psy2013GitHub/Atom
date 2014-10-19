

function Atom_movefile(Incell,Dir,Options)

if nargin<3, Options=struct(); end
Style='Raw'; if isfield(Options,'Style')&&~isempty(Options.Style), Style=Options.Style; end

for ii=1:length(Incell)
    fprintf('.');
    if ischar(Incell{ii})
        copyfile(Incell{ii},Dir);
    else
        if iscell(Incell{ii})
            for ses=1:length(Incell{ii})
                if ~strcmpi(Style,'Raw'); tmpDir=[Dir,filesep,Style,'_',num2str(ses)]; else, tmpDir=[Dir,filesep,]; end
                mkdir(tmpDir);
                Atom_movefile(Incell{ii}(ses),tmpDir);
            end
        end
    end
end
end