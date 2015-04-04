

function reg = design_specific_regressor(total_nImg_lst, nImgDiscard_pre, nImgDiscard_pro, RT, regressor_cell)

reg = struct('name',[], 'val',[]);

nRun = length(total_nImg_lst);

remained_Img_lst = total_nImg_lst - nImgDiscard_pre - nImgDiscard_pro;

Img_range = [0, cumsum(remained_Img_lst)];

TotalImgNumber = Img_range(end);

for reg_idx = 1:length(regressor_cell)
    
    curr_reg = regressor_cell{reg_idx};
    
    if strcmpi(curr_reg, 'linear_trend')
        for run_idx = 1:nRun
            nImg_run = remained_Img_lst(run_idx);
            
            reg(end+1).name = ['linear_trend_', num2str(run_idx)];
            fill_val = 1:nImg_run;
            reg(end).val = get_reg_val(fill_val, run_idx,TotalImgNumber, Img_range);
        end
    elseif strcmpi(curr_reg, 'nonlinear_fourier_trend')
        for run_idx = 1:nRun
            nImg_run = remained_Img_lst(run_idx);
            
            % sin
            reg(end+1).name = ['sin1_trend_', num2str(run_idx)];
            fill_val = sin(2*pi/(nImg_run*RT) * (1:nImg_run));
            reg(end).val = get_reg_val(fill_val, run_idx, TotalImgNumber, Img_range);
            
            reg(end+1).name = ['sin2_trend_', num2str(run_idx)];
            fill_val = sin(4*pi/(nImg_run*RT) * (1:nImg_run));
            reg(end).val = get_reg_val(fill_val, run_idx, TotalImgNumber, Img_range);
            
            reg(end+1).name = ['sin3_trend_', num2str(run_idx)];
            fill_val = sin(6*pi/(nImg_run*RT) * (1:nImg_run));
            reg(end).val = get_reg_val(fill_val, run_idx, TotalImgNumber, Img_range);
            % cos
            reg(end+1).name = ['cos1_trend_', num2str(run_idx)];
            fill_val = cos(2*pi/(nImg_run*RT) * (1:nImg_run));
            reg(end).val = get_reg_val(fill_val, run_idx, TotalImgNumber, Img_range);
            
            reg(end+1).name = ['cos2_trend_', num2str(run_idx)];
            fill_val = cos(4*pi/(nImg_run*RT) * (1:nImg_run));
            reg(end).val = get_reg_val(fill_val, run_idx, TotalImgNumber, Img_range);
            
            reg(end+1).name = ['cos3_trend_', num2str(run_idx)];
            fill_val = cos(6*pi/(nImg_run*RT) * (1:nImg_run));
            reg(end).val = get_reg_val(fill_val, run_idx, TotalImgNumber, Img_range);
            
        end
        
    elseif strcmpi(curr_reg, 'cofound_mean')
        for run_idx = 1:nRun
            nImg_run = remained_Img_lst(run_idx);
            
            reg(end+1).name = ['cofound_mean_', num2str(run_idx)];
            fill_val = ones(nImg_run, 1);
            reg(end).val = get_reg_val(fill_val, run_idx, TotalImgNumber, Img_range);
        end
        
    else
        
        error('Unsuported Yet!');
    end
    
end

try
    reg = reg(2:end);
catch
    reg = [];
end

return
end

function [left, right] = get_range_idx(img_range, run_idx)

left = img_range(run_idx) + 1;
right = img_range(run_idx+1);

end


function reg_val = get_reg_val(fill_val, run_idx, total_img_number, img_range)

reg_val = zeros(total_img_number, 1);

[left, right] = get_range_idx(img_range, run_idx);

reg_val(left:right) = fill_val;

end
