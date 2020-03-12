function [cellSimParams,desired_state_reached,attempts,cell_reset] = driveToSoc(cellSimParams,des_soc_bin_idx)
model = cellSimParams.model;
soc_grid_boundaries = cellSimParams.soc_grid_boundaries;
soc_num = length(soc_grid_boundaries)-1;
cell_pow_set = cellSimParams.cell_pow_set;
pow_num = length(cell_pow_set);

cur_soc = cellSimParams.cur_soc_coulombic; % Cheating a bit here
cur_soc_bin_idx = find(soc_grid_boundaries(2:end-1)>cur_soc,1);
if(isempty(cur_soc_bin_idx))
    cur_soc_bin_idx = soc_num;
end

if(cur_soc_bin_idx == des_soc_bin_idx)
    desired_state_reached = 1;
    cellSimParams.cur_soc_3c = cur_soc;        
else
    desired_state_reached = 0;
end

des_SOC_low = soc_grid_boundaries(des_soc_bin_idx);
des_SOC_high = soc_grid_boundaries(des_soc_bin_idx+1);
des_SOC = (des_SOC_low + des_SOC_high)/2;
model.param.set('des_SOC',des_SOC);

cell_reset = 0;
charging_set = 0;
discharging_set = 0;
attempts = 0;
attempts_max = cellSimParams.driveToSOC_attempts_max;
cycles = 0;
timePeriodScale = cellSimParams.driveToSOC_timePeriodScaleFactor;
timePeriod = (timePeriodScale^cycles)*cellSimParams.slotIntervalInSeconds*cellSimParams.driveToSOC_timeAccelerationFactor;
model.param.set('period',timePeriod);
model.param.set('t_factor',cellSimParams.driveToSOC_timeAccelerationFactor);


while(~desired_state_reached)
    disp(strcat({'driving SOC from '},num2str(cur_soc),{' to ['},num2str(des_SOC_low),',',num2str(des_SOC_high),'] with attempt #',num2str(attempts),'...'));
    if(cur_soc_bin_idx<des_soc_bin_idx)
        ref_pow_idx = (randi(pow_num));
        while(cell_pow_set(ref_pow_idx)<=0)
            ref_pow_idx = (randi(pow_num));
        end
                
        disp(strcat({'charging in progress with power: '},num2str(cell_pow_set(ref_pow_idx))));
        if(~charging_set)
            model.sol('sol2').feature('t1').feature('st1').setIndex('stopcondActive', true, 2); % charging
            model.sol('sol2').feature('t1').feature('st1').setIndex('stopcondActive', false, 3); % discharging
            charging_set = 1;
        end
        if(discharging_set)
            discharging_set = 0;
            cycles = cycles + 1;
            timePeriod = (timePeriodScale^cycles)*cellSimParams.slotIntervalInSeconds*cellSimParams.driveToSOC_timeAccelerationFactor;
            model.param.set('period',timePeriod);
        end
    else
        ref_pow_idx = (randi(pow_num));
        while(cell_pow_set(ref_pow_idx)>=0)
            ref_pow_idx = (randi(pow_num));
        end
        
        disp(strcat({'discharging in progress with power: '},num2str(cell_pow_set(ref_pow_idx))));
        if(~discharging_set)
            model.sol('sol2').feature('t1').feature('st1').setIndex('stopcondActive', false, 2); % charging
            model.sol('sol2').feature('t1').feature('st1').setIndex('stopcondActive', true, 3); % discharging
            discharging_set = 1;
        end        
        if(charging_set)
            charging_set = 0;
            cycles = cycles + 1;
            timePeriod = (timePeriodScale^cycles)*cellSimParams.slotIntervalInSeconds*cellSimParams.driveToSOC_timeAccelerationFactor;
            model.param.set('period',timePeriod);
        end
    end
    
    [out,cell_reset] = simulateBatteryCell(cellSimParams,ref_pow_idx,0);
    if(cell_reset)
        attempts = 1;
    end
    cellSimParams = out.cellSimParams;
    cur_soc = cellSimParams.cur_soc_coulombic; % Cheating a bit here
    cur_soc_bin_idx = find(soc_grid_boundaries(2:end-1)>cur_soc,1);
    if(isempty(cur_soc_bin_idx))
        cur_soc_bin_idx = soc_num;
    end
    
    if(cur_soc_bin_idx == des_soc_bin_idx)
        desired_state_reached = 1;
        cellSimParams.cur_soc_3c = cur_soc;
    elseif(cell_pow_set(ref_pow_idx)~=0)
        attempts = attempts + 1;
        if(attempts>attempts_max)
            cellSimParams.cur_soc_3c = cur_soc;
            break;
        end
    end
end

if(charging_set)
    model.sol('sol2').feature('t1').feature('st1').setIndex('stopcondActive', false, 2); % charging
end
if(discharging_set)
    model.sol('sol2').feature('t1').feature('st1').setIndex('stopcondActive', false, 3); % discharging
end
model.param.set('des_SOC',1);
model.param.set('period',cellSimParams.slotIntervalInSeconds);
model.param.set('t_factor',1);
if(cur_soc_bin_idx~=des_soc_bin_idx)
    warning(strcat('want: ',num2str(des_soc_bin_idx),' acheived:',num2str(cur_soc_bin_idx)));
end
end