

function regressor = merge_rp_para(root, dir_lst)

regressor = struct('name',[],'val',[]);

nValidDir = 0;

rp_data_lst = {}; nr_lst = []; 
for d = 1:length(dir_lst)

    dir_path = [root, filesep, dir_lst(d).name];
    rp_txt = spm_select('List', dir_path, '^rp.*\.txt$');
    
    
    if isempty(rp_txt)
        continue
    end
    
    tmp_rp_data = load([dir_path, filesep, rp_txt]);
    rp_data_lst = [rp_data_lst, {tmp_rp_data}];
    nr_lst = [nr_lst, size(tmp_rp_data, 1)];
    
    nValidDir = nValidDir + 1;

end

Img_range = [0, cumsum(nr_lst)];
TotalRows = sum(nr_lst);
for d = 1:nValidDir
    
    left = Img_range(d)+1;
    right = Img_range(d+1);
    
    for c = 1:6
        reg.name = ['run_',num2str(d),'_motion_',num2str(c)];
        val = zeros(TotalRows,1);
        val(left:right) = rp_data_lst{d}(:,c);  
        reg.val = val;
        regressor(end+1) = reg;
    end
try    
    regressor = regressor(2:end);
catch
    regressor = [];
end

return
end






























return
end