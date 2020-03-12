function [cellSimData] = addCalenderAgingAndAaggregateSamples(cellSimData_all_sims,cellSimParams,soc_grid_bin_mean)
simNum = length(cellSimData_all_sims);
soc_num = size(cellSimData_all_sims{1}.energy_loss_samples,1);
pow_num = size(cellSimData_all_sims{1}.energy_loss_samples,2);
sample_num_per_sim = size(cellSimData_all_sims{1}.energy_loss_samples,3);
sample_num = sample_num_per_sim*simNum;
slotIntervalInSeconds = cellSimParams.slotIntervalInSeconds;
slotIntervalInHrs = slotIntervalInSeconds/3600;
experimental_socs = [0;0.5;1];
experimental_rel_cap_at_900_calenderDays =[0.965;0.908;0.865];
experimental_calender_capacity_loss_factor_per_slot =  (1-experimental_rel_cap_at_900_calenderDays)/(900*24/slotIntervalInHrs);
interpolated_calender_capacity_loss_factor_per_slot = interp1(experimental_socs,experimental_calender_capacity_loss_factor_per_slot,soc_grid_bin_mean);

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
cellsIterated = 0;
sample_num_offset = 0;
for simIdx = 1:simNum    
    simTimeRatio_samples(:,:,sample_num_offset+1:sample_num_offset+sample_num_per_sim) = cellSimData_all_sims{simIdx}.simTimeRatio_samples;
    relative_capacity_samples(:,:,sample_num_offset+1:sample_num_offset+sample_num_per_sim) = cellSimData_all_sims{simIdx}.relative_capacity_samples;
    energy_loss_samples(:,:,sample_num_offset+1:sample_num_offset+sample_num_per_sim) = cellSimData_all_sims{simIdx}.energy_loss_samples;
    energy_applied_samples(:,:,sample_num_offset+1:sample_num_offset+sample_num_per_sim) = cellSimData_all_sims{simIdx}.energy_applied_samples;    
    capacity_loss_factor_samples(:,:,sample_num_offset+1:sample_num_offset+sample_num_per_sim) = cellSimData_all_sims{simIdx}.capacity_loss_factor_samples;
    soc_k_samples_at_load(:,:,sample_num_offset+1:sample_num_offset+sample_num_per_sim) = cellSimData_all_sims{simIdx}.soc_k_samples_at_load;
    soc_kp1_samples_at_load(:,:,sample_num_offset+1:sample_num_offset+sample_num_per_sim) = cellSimData_all_sims{simIdx}.soc_kp1_samples_at_load;
    soc_k_samples_coulombic(:,:,sample_num_offset+1:sample_num_offset+sample_num_per_sim) = cellSimData_all_sims{simIdx}.soc_k_samples_coulombic;
    soc_kp1_samples_coulombic(:,:,sample_num_offset+1:sample_num_offset+sample_num_per_sim) = cellSimData_all_sims{simIdx}.soc_kp1_samples_coulombic;
    soc_k_samples_3c(:,:,sample_num_offset+1:sample_num_offset+sample_num_per_sim) = cellSimData_all_sims{simIdx}.soc_k_samples_3c;
    soc_kp1_samples_3c(:,:,sample_num_offset+1:sample_num_offset+sample_num_per_sim) = cellSimData_all_sims{simIdx}.soc_kp1_samples_3c;
    terminal_voltage_samples(:,:,sample_num_offset+1:sample_num_offset+sample_num_per_sim) = cellSimData_all_sims{simIdx}.terminal_voltage_samples;
    internal_resistance_samples(:,:,sample_num_offset+1:sample_num_offset+sample_num_per_sim) = cellSimData_all_sims{simIdx}.internal_resistance_samples;
    attempts_to_driveToSOC_samples(:,:,sample_num_offset+1:sample_num_offset+sample_num_per_sim) = cellSimData_all_sims{simIdx}.attempts_to_driveToSOC_samples;
    cellsIterated = cellsIterated + cellSimData_all_sims{simIdx}.cellsIterated;
    
    sample_num_offset = sample_num_offset + sample_num_per_sim;
end

capacity_loss_factor_incl_calender_samples = Inf([soc_num,pow_num,sample_num]);
for soc_idx = 1:soc_num
    capacity_loss_factor_incl_calender_samples(soc_idx,:,:) = capacity_loss_factor_samples(soc_idx,:,:);
    for pow_idx = 1:pow_num
        for sample_idx = 1: sample_num
            if(~isnan(capacity_loss_factor_incl_calender_samples(soc_idx,pow_idx,sample_idx)))
                capacity_loss_factor_incl_calender_samples(soc_idx,pow_idx,sample_idx) = max(capacity_loss_factor_incl_calender_samples(soc_idx,pow_idx,sample_idx),...
                    interpolated_calender_capacity_loss_factor_per_slot(soc_idx));
            end
        end
    end    
end


cellSimData = struct;
cellSimData.simTimeRatio_samples = simTimeRatio_samples;
cellSimData.relative_capacity_samples = relative_capacity_samples;
cellSimData.energy_loss_samples = energy_loss_samples;
cellSimData.energy_applied_samples = energy_applied_samples;
cellSimData.capacity_loss_factor_samples = capacity_loss_factor_samples;
cellSimData.capacity_loss_factor_incl_calender_samples = capacity_loss_factor_incl_calender_samples;
cellSimData.soc_k_samples_at_load = soc_k_samples_at_load;
cellSimData.soc_k_samples_coulombic = soc_k_samples_coulombic;
cellSimData.soc_k_samples_3c = soc_k_samples_3c;
cellSimData.soc_kp1_samples_at_load = soc_kp1_samples_at_load;
cellSimData.soc_kp1_samples_coulombic = soc_kp1_samples_coulombic;
cellSimData.soc_kp1_samples_3c = soc_kp1_samples_3c;
cellSimData.internal_resistance_samples = internal_resistance_samples;
cellSimData.terminal_voltage_samples = terminal_voltage_samples;
cellSimData.attempts_to_driveToSOC_samples = attempts_to_driveToSOC_samples;
cellSimData.cellsIterated = cellsIterated;
end

