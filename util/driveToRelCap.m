function [cellSimParams] = driveToRelCap(cellSimParams,des_cellRelativeCapacity)
cell_pow_set = cellSimParams.cell_pow_set;
pow_num = length(cell_pow_set);

model = cellSimParams.model;
cur_cellRelativeCapacity = cellSimParams.currentCellRelativeCapacity;

if(cur_cellRelativeCapacity <= des_cellRelativeCapacity)
    desired_state_reached = 1;
else
    desired_state_reached = 0;
    model.param.set('period',cellSimParams.slotIntervalInSeconds*cellSimParams.driveToSOH_timeAccelerationFactor);
    model.param.set('t_factor',cellSimParams.driveToSOH_timeAccelerationFactor);
    
    soc_k = cellSimParams.cur_soc_3c;
    soc_grid_boundaries = cellSimParams.soc_grid_boundaries;
    soc_num = length(soc_grid_boundaries)-1;
    soc_k_bin_idx = find(soc_grid_boundaries(2:end-1)>soc_k,1);
    if(isempty(soc_k_bin_idx))
        soc_k_bin_idx = soc_num;
    end
end

while(~desired_state_reached)
    disp(strcat({'driving SOH from '},num2str(cur_cellRelativeCapacity),{' to '},num2str(des_cellRelativeCapacity),'...'));
    
    if(soc_k_bin_idx==1)
        ref_pow_idx = (randi(pow_num));
        while(cell_pow_set(ref_pow_idx)<=0)
            ref_pow_idx = (randi(pow_num));
        end
    elseif(soc_k_bin_idx==soc_num)
        ref_pow_idx = (randi(pow_num));
        while(cell_pow_set(ref_pow_idx)>=0)
            ref_pow_idx = (randi(pow_num));
        end
    else
        ref_pow_idx = (randi(pow_num));
    end    
    
    out = simulateBatteryCell(cellSimParams,ref_pow_idx,1);
    cellSimParams = out.cellSimParams;
    cur_cellRelativeCapacity = cellSimParams.currentCellRelativeCapacity;
    
    if(cur_cellRelativeCapacity <= des_cellRelativeCapacity)
        desired_state_reached = 1;
    else
        desired_state_reached = 0;
        soc_k = cellSimParams.cur_soc_3c;
        soc_grid_boundaries = cellSimParams.soc_grid_boundaries;
        soc_num = length(soc_grid_boundaries)-1;
        soc_k_bin_idx = find(soc_grid_boundaries(2:end-1)>soc_k,1);
        if(isempty(soc_k_bin_idx))
            soc_k_bin_idx = soc_num;
        end
    end
end
model.param.set('period',cellSimParams.slotIntervalInSeconds);
model.param.set('t_factor',1);

cell_pow_set = cellSimParams.cell_pow_set;
zero_pow_idx = find(cell_pow_set == 0,1);
out = simulateBatteryCell(cellSimParams,zero_pow_idx,1); % keep idle for one time slot
cellSimParams = out.cellSimParams;
end