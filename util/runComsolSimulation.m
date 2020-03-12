function [cellSimData] = runComsolSimulation(cellSimParams,comsol_model_filename)
soc_grid_boundaries = cellSimParams.soc_grid_boundaries;
cell_pow_set = cellSimParams.cell_pow_set;

sample_num = cellSimParams.sample_num;

soc_num = length(soc_grid_boundaries)-1;
pow_num = length(cell_pow_set);

energy_loss_samples = Inf([soc_num,pow_num,sample_num]);
energy_applied_samples = Inf([soc_num,pow_num,sample_num]);
capacity_loss_factor_samples = Inf([soc_num,pow_num,sample_num]);
soc_k_samples_at_load = Inf([soc_num,pow_num,sample_num]);
soc_kp1_samples_at_load = Inf([soc_num,pow_num,sample_num]);
soc_k_samples_coulombic= Inf([soc_num,pow_num,sample_num]);
soc_kp1_samples_coulombic = Inf([soc_num,pow_num,sample_num]);
soc_k_samples_3c= Inf([soc_num,pow_num,sample_num]);
soc_kp1_samples_3c = Inf([soc_num,pow_num,sample_num]);
simTimeRatio_samples = Inf([soc_num,pow_num,sample_num]);
terminal_voltage_samples = Inf([soc_num,pow_num,sample_num]);
internal_resistance_samples = Inf([soc_num,pow_num,sample_num]);
relative_capacity_samples = Inf([soc_num,pow_num,sample_num]);
attempts_to_driveToSOC_samples = Inf([soc_num,pow_num,sample_num]);

map_fill_count = 0;
map_fill_count_max = soc_num*pow_num*sample_num;
validDataCount = 0;
invalidDataCount = 0;

try
    mphopen -clear;
catch ME1
    message0 = 'Cannot find COMSOL server';
    message1 = 'A connection to COMSOL could not be established';
    if(contains(ME1.message,message0))
        try
            mphstart
        catch ME2
            message2 = 'Already connected to a server';
            if(contains(ME2.message,message2))
                warning(message2);
            else
                warning(message1);
            end
        end
    else
        warning(message1);
    end    
end

try
    model = mphopen(comsol_model_filename);
catch ME1
    message0 = 'Cannot find COMSOL server';
    if(contains(ME1.message,message0))
        try
            mphstart;
            import com.comsol.model.util.*;
            model = mphopen(comsol_model_filename);
        catch ME2
            message1 = 'A connection to COMSOL could not be established';
            message2 = 'Already connected to a server';
            if(contains(ME2.message,message2))
                warning(message2);
            else
                error(message1);
            end
        end
    else
        rethrow(ME1);
    end    
end

import com.comsol.model.util.*;

cellSimParams.cellsIterated = 0;
cellSimParams.model = model;

[cellSimParams] = initializeNewCell(cellSimParams);
des_soc_bin_idx = findUnfilledSoC(energy_loss_samples);
desired_state_reached = 0;
attempts_to_driveToSOC = 0;
while(~desired_state_reached)
    [cellSimParams,desired_state_reached,attempts,cell_reset] = driveToSoc(cellSimParams,des_soc_bin_idx);
    if(~desired_state_reached)
        des_soc_bin_idx = findUnfilledSoC(energy_loss_samples);        
    end
    if(cell_reset)
        attempts_to_driveToSOC = attempts;
    else
        attempts_to_driveToSOC = attempts_to_driveToSOC + attempts;
    end
end
soc_k_bin_idx = des_soc_bin_idx;
mapFilled = 0;
while(~mapFilled)
    ref_pow_idx = findUnfilledPow(energy_loss_samples,soc_k_bin_idx);
    changed_flag = 0;
    if(isempty(ref_pow_idx))
        des_soc_bin_idx = findUnfilledSoC(energy_loss_samples);
        if(isempty(des_soc_bin_idx))
            mapFilled = 1;
            break % map filled
        end
        desired_state_reached = 0;
        attempts_to_driveToSOC = 0;
        while(~desired_state_reached)
            [cellSimParams,desired_state_reached,attempts,cell_reset] = driveToSoc(cellSimParams,des_soc_bin_idx);
            if(~desired_state_reached)
                des_soc_bin_idx = findUnfilledSoC(energy_loss_samples);
            end
            if(cell_reset)
                attempts_to_driveToSOC = attempts;
            else
                attempts_to_driveToSOC = attempts_to_driveToSOC + attempts;
            end
%             if(attempts>cellSimParams.driveToSOC_attempts_max)
%                 keyboard
%             end
        end
        soc_k_bin_idx = des_soc_bin_idx;
        changed_flag = 1;
        ref_pow_idx = findUnfilledPow(energy_loss_samples,soc_k_bin_idx);
    end
    sample_idx = find(isinf(energy_loss_samples(soc_k_bin_idx,ref_pow_idx,:)), 1);
    attempts_to_driveToSOC_samples(soc_k_bin_idx,ref_pow_idx,sample_idx) = attempts_to_driveToSOC;
    soc_k_samples_at_load(soc_k_bin_idx,ref_pow_idx,sample_idx) = cellSimParams.cur_soc_at_load;
    soc_k_samples_coulombic(soc_k_bin_idx,ref_pow_idx,sample_idx) = cellSimParams.cur_soc_coulombic;
    soc_k_samples_3c(soc_k_bin_idx,ref_pow_idx,sample_idx) = cellSimParams.cur_soc_3c;
    
    out = simulateBatteryCell(cellSimParams,ref_pow_idx,0);    
    
    relative_capacity_samples(soc_k_bin_idx,ref_pow_idx,sample_idx) = cellSimParams.currentCellRelativeCapacity;
    validData = 0;    
    simTimeRatio = out.simTimeRatio;
    cellSimParams = out.cellSimParams;
    if(out.simTimeRatio>0)    
        capacity_loss_factor_samples(soc_k_bin_idx,ref_pow_idx,sample_idx) = out.capacity_loss_factor; %Ah
        energy_loss_samples(soc_k_bin_idx,ref_pow_idx,sample_idx) = out.energy_loss; %Wh   
        simTimeRatio_samples(soc_k_bin_idx,ref_pow_idx,sample_idx) = simTimeRatio;
        internal_resistance_samples(soc_k_bin_idx,ref_pow_idx,sample_idx) = out.mean_internal_resistance;
        terminal_voltage_samples(soc_k_bin_idx,ref_pow_idx,sample_idx) = out.mean_terminal_voltage;
        energy_applied_samples(soc_k_bin_idx,ref_pow_idx,sample_idx) = out.energy_applied;        
        
        soc_kp1_samples_at_load(soc_k_bin_idx,ref_pow_idx,sample_idx) = cellSimParams.cur_soc_at_load;
        soc_kp1_samples_coulombic(soc_k_bin_idx,ref_pow_idx,sample_idx) = cellSimParams.cur_soc_coulombic;
        soc_kp1_samples_3c(soc_k_bin_idx,ref_pow_idx,sample_idx) = cellSimParams.cur_soc_3c;
        
        soc_kp1_bin_idx = find(soc_grid_boundaries(2:end-1)>cellSimParams.cur_soc_3c,1);
        if(isempty(soc_kp1_bin_idx))
            soc_kp1_bin_idx = soc_num;
        end
        soc_k_bin_idx = soc_kp1_bin_idx;
        
        validData = 1;
        validDataCount = validDataCount + 1;
        map_fill_count = map_fill_count + 1;
    else
        energy_loss_samples(soc_k_bin_idx,ref_pow_idx,sample_idx) = nan;
        energy_applied_samples(soc_k_bin_idx,ref_pow_idx,sample_idx) = nan;
        capacity_loss_factor_samples(soc_k_bin_idx,ref_pow_idx,sample_idx) = nan;
        soc_kp1_samples_at_load(soc_k_bin_idx,ref_pow_idx,sample_idx) = nan;
        soc_kp1_samples_coulombic(soc_k_bin_idx,ref_pow_idx,sample_idx) = nan;
        soc_kp1_samples_3c(soc_k_bin_idx,ref_pow_idx,sample_idx) = nan;
        internal_resistance_samples(soc_k_bin_idx,ref_pow_idx,sample_idx) = nan;
        terminal_voltage_samples(soc_k_bin_idx,ref_pow_idx,sample_idx) = nan;
        map_fill_count = map_fill_count + 1;
        invalidDataCount = invalidDataCount + 1;
    end
    disp(strcat({' simTimeRatio: '},num2str(simTimeRatio),...
        {', Valid: '},num2str(validData),...
        {', ValidDataCount: '},num2str(validDataCount),...
        {', InvalidDataCount: '},num2str(invalidDataCount),...
        {', Filled map ratio: '},num2str(map_fill_count),{'/'},num2str(map_fill_count_max)));
end

cellSimData = struct;
cellSimData.simTimeRatio_samples = simTimeRatio_samples;
cellSimData.relative_capacity_samples = relative_capacity_samples;
cellSimData.energy_loss_samples = energy_loss_samples;
cellSimData.energy_applied_samples = energy_applied_samples;
cellSimData.capacity_loss_factor_samples = capacity_loss_factor_samples;
cellSimData.soc_k_samples_at_load = soc_k_samples_at_load;
cellSimData.soc_k_samples_coulombic = soc_k_samples_coulombic;
cellSimData.soc_k_samples_3c = soc_k_samples_3c;
cellSimData.soc_kp1_samples_at_load = soc_kp1_samples_at_load;
cellSimData.soc_kp1_samples_coulombic = soc_kp1_samples_coulombic;
cellSimData.soc_kp1_samples_3c = soc_kp1_samples_3c;
cellSimData.internal_resistance_samples = internal_resistance_samples;
cellSimData.terminal_voltage_samples = terminal_voltage_samples;
cellSimData.attempts_to_driveToSOC_samples = attempts_to_driveToSOC_samples;
cellSimData.cellsIterated = cellSimParams.cellsIterated;
end