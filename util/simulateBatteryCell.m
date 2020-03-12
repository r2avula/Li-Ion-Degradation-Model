function [out,cell_reset] = simulateBatteryCell(cellSimParams,ref_pow_idx,drivingToDesiredCapacity)
cell_reset = 0;
if(cellSimParams.initialRelCap-cellSimParams.currentCellRelativeCapacity > cellSimParams.allowedRelativeCapacityChange)
    [cellSimParams] = initializeNewCell(cellSimParams);
    cell_reset = 1;
end
cell_1C_capacityInAh = cellSimParams.cell_1C_capacityInAh;
cell_nominalVoltage = cellSimParams.cell_nominalVoltage;
cur_soc_3c = cellSimParams.cur_soc_3c;
cur_soc_coulombic = cellSimParams.cur_soc_coulombic;

cur_voltage = cellSimParams.cur_voltage;
slotIntervalInSeconds = cellSimParams.slotIntervalInSeconds;
model = cellSimParams.model;
cell_pow_set = cellSimParams.cell_pow_set;
ref_pow = cell_pow_set(ref_pow_idx);

soc_limit_buffer = 0.05;
volt_limit_buffer = 0.2;

model.param.set('P_ref',ref_pow);

run_sim = 1;
if(ref_pow>0 && cur_soc_coulombic >= cellSimParams.SOC_high) || (ref_pow<0 && cur_soc_coulombic <= cellSimParams.SOC_low)
    run_sim = 0;
end

SOC_low_changed = 0;
E_min_changed = 0;
SOC_high_changed = 0;
E_max_changed = 0;
something_changed = 0;
run_sim_executed = 0;
capacity_loss_factor = 0;
energy_loss = 0;
energy_applied = 0;
mean_terminal_voltage = 0;
mean_internal_resistance = 0;
simTimeRatio = 0;

error_message = 'Simulation failed to start!';

if(run_sim)  
    fprintf('Running COMSOL simulation ...');
    try
        model.study('std2').run;% Apply requested power
        run_sim_executed = 1;
    catch ME
        error_message = 'Stop condition fulfilled for initial values';
        if(contains(ME.message,error_message))
            if(ref_pow>0)
                if(cur_soc_coulombic <= cellSimParams.SOC_low+soc_limit_buffer)
                    model.param.set('SOC_low',max(cur_soc_coulombic-soc_limit_buffer,0));
                    SOC_low_changed = 1;
                    something_changed = 1;
                end
                if(cur_voltage <= cellSimParams.cell_voltage_low+volt_limit_buffer)
                    model.param.set('E_min',cur_voltage-volt_limit_buffer);
                    E_min_changed = 1;
                    something_changed = 1;
                end
            end
            if(ref_pow<=0)
                if(cur_soc_coulombic >= cellSimParams.SOC_high-soc_limit_buffer)
                    model.param.set('SOC_high',min(cur_soc_coulombic+soc_limit_buffer,1));
                    SOC_high_changed = 1;
                    something_changed = 1;
                end
                if(cur_voltage >= cellSimParams.cell_voltage_high-volt_limit_buffer)
                    model.param.set('E_max',cur_voltage+volt_limit_buffer);
                    E_max_changed = 1;
                    something_changed = 1;
                end
            end
            if(something_changed)
                try
                    model.study('std2').run;% Apply requested power
                    run_sim_executed = 1;
                catch ME2
                    if(~contains(ME2.message,error_message))
                        error_message = strcat(error_message,' + ',ME2.message);
                    end
                end
            end            
        else
            rethrow(ME);
        end
    end
end

if(run_sim_executed)
    fprintf('SUCCESS');
    simTimeInSeconds = floor(mphglobal(model,'t','dataset','dset2','solnum','end'));
    simTimeRatio = (simTimeInSeconds)/slotIntervalInSeconds;    
    
    cellSimParams.cur_soc_at_load = mphglobal(model,'SOCcell_load','dataset','dset2','solnum','end');
    cellSimParams.cur_soc_coulombic = mphglobal(model,'SOCcell','dataset','dset2','solnum','end');
    cellSimParams.cur_voltage = mphglobal(model,'Ecell','dataset','dset2','solnum','end');  
    
    soc_grid_boundaries = cellSimParams.soc_grid_boundaries;
    soc_num = length(soc_grid_boundaries)-1;
    soc_k_coulombic_bin_idx = find(soc_grid_boundaries(2:end-1)>cellSimParams.cur_soc_coulombic,1);
    if(isempty(soc_k_coulombic_bin_idx))
        soc_k_coulombic_bin_idx = soc_num;
    end    
    
    if(drivingToDesiredCapacity)
        cellSimParams.currentCellRelativeCapacity = mphglobal(model,'Rel_cap','dataset','dset2','solnum','end')*100;
        cellSimParams.cur_soc_3c = cellSimParams.cur_soc_coulombic;         % Cheating a bit here
    else
        old_rel_cap = cellSimParams.currentCellRelativeCapacity;
        new_rel_cap = mphglobal(model,'Rel_cap','dataset','dset2','solnum','end')*100;
        capacity_loss_factor = (old_rel_cap-new_rel_cap)/100/simTimeRatio; %in A
        try
            energy_loss = mphglobal(model,strcat('timeint(0,',num2str(simTimeInSeconds),',Power_loss)'),'dataset','dset2','solnum','end')/3600; %in Wh
            energy_applied = mphglobal(model,strcat('timeint(0,',num2str(simTimeInSeconds),',Power_act)'),'dataset','dset2','solnum','end')/3600; %in Wh
            mean_terminal_voltage = mphglobal(model,strcat('timeint(0,',num2str(simTimeInSeconds),',Ecell)'),'dataset','dset2','solnum','end')/simTimeInSeconds; %in V
            mean_internal_resistance = mphglobal(model,strcat('timeint(0,',num2str(simTimeInSeconds),',R_cell_int)'),'dataset','dset2','solnum','end')/simTimeInSeconds; %in ohm
            
            if(soc_k_coulombic_bin_idx == 1)
                soc_kp1_estimate_3c = cellSimParams.SOC_low;    %% SOC reset at Low Voltage
            elseif(soc_k_coulombic_bin_idx == soc_num)
                soc_kp1_estimate_3c = cellSimParams.SOC_high;    %% SOC reset at high Voltage
            else
                cell_cur_en_cap = old_rel_cap*cell_1C_capacityInAh*cell_nominalVoltage/100;
                gamma_tau = simTimeInSeconds/3600; % in h
                soc_kp1_estimate_3c = cur_soc_3c +...
                    ((gamma_tau*mean_terminal_voltage/2/mean_internal_resistance)*(sqrt(max((mean_terminal_voltage*mean_terminal_voltage)+4*mean_internal_resistance*ref_pow,0))-mean_terminal_voltage))/(cell_cur_en_cap);
                soc_kp1_estimate_3c = min(max(soc_kp1_estimate_3c,cellSimParams.SOC_low),cellSimParams.SOC_high);
            end
            cellSimParams.cur_soc_3c = soc_kp1_estimate_3c;            
        catch ME
            if(contains(ME.message,'timeint'))
                disp('Warning: Timeint not performed.');
            else
                rethrow(ME);
            end
        end
        
        cellSimParams.currentCellRelativeCapacity = new_rel_cap;
    end
    
    cellSimParams.model = model;
else    
    fprintf(strcat('FAIL: ',error_message));
end

fprintf('\n');
if(SOC_low_changed)
    model.param.set('SOC_low',cellSimParams.SOC_low);
end
if(E_min_changed)
    model.param.set('E_min',cellSimParams.cell_voltage_low);
end
if(SOC_high_changed)
    model.param.set('SOC_high',cellSimParams.SOC_high);
end
if(E_max_changed)
    model.param.set('E_max',cellSimParams.cell_voltage_high);
end

out = struct;
out.cellSimParams = cellSimParams;
out.capacity_loss_factor = capacity_loss_factor;
out.energy_loss = energy_loss;
out.energy_applied = energy_applied;
out.mean_terminal_voltage = mean_terminal_voltage;
out.mean_internal_resistance = mean_internal_resistance;
out.simTimeRatio = simTimeRatio;
end