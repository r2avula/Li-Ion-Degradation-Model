function [out,cell_reset] = simulateBatteryCell(cellSimParams,ref_pow_idx,drivingToDesiredState)
cell_reset = 0;
if(cellSimParams.initialRelCap-cellSimParams.currentCellRelativeCapacity > cellSimParams.allowedRelativeCapacityChange)
    [cellSimParams] = initializeNewCell(cellSimParams);
    cell_reset = 1;
end

cur_soc = cellSimParams.cur_soc;
cur_voltage = cellSimParams.cur_voltage;
slotIntervalInSeconds = cellSimParams.slotIntervalInSeconds;
model = cellSimParams.model;
cell_pow_set = cellSimParams.cell_pow_set;
ref_pow = cell_pow_set(ref_pow_idx);

model.param.set('P_ref',ref_pow);

run_sim = 1;
if(ref_pow>0 && cur_soc >= cellSimParams.SOC_high) || (ref_pow<0 && cur_soc <= cellSimParams.SOC_low)
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
                if(cur_soc <= cellSimParams.SOC_low+0.05)
                    model.param.set('SOC_low',max(cur_soc-0.05,0.05));
                    SOC_low_changed = 1;
                    something_changed = 1;
                end
                if(cur_voltage <= cellSimParams.cell_voltage_low+0.2)
                    model.param.set('E_min',max(cur_voltage-0.2,2.2));
                    E_min_changed = 1;
                    something_changed = 1;
                end
            end
            if(ref_pow<=0)
                if(cur_soc >= cellSimParams.SOC_high-0.05)
                    model.param.set('SOC_high',min(cur_soc+0.05,0.95));
                    SOC_high_changed = 1;
                    something_changed = 1;
                end
                if(cur_voltage >= cellSimParams.cell_voltage_high-0.2)
                    model.param.set('E_max',min(cur_voltage+0.2,4.5));
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
    if(drivingToDesiredState)
        cellSimParams.currentCellRelativeCapacity = mphglobal(model,'Rel_cap','dataset','dset2','solnum','end')*100;
    else
        old_rel_cap = cellSimParams.currentCellRelativeCapacity;
        new_rel_cap = mphglobal(model,'Rel_cap','dataset','dset2','solnum','end')*100;
        capacity_loss_factor = (old_rel_cap-new_rel_cap)/100/simTimeRatio; %in A
        try
            energy_loss = mphglobal(model,strcat('timeint(0,',num2str(simTimeInSeconds),',Power_loss)'),'dataset','dset2','solnum','end')/3600; %in Wh
            energy_applied = mphglobal(model,strcat('timeint(0,',num2str(simTimeInSeconds),',Power_act)'),'dataset','dset2','solnum','end')/3600; %in Wh
            mean_terminal_voltage = mphglobal(model,strcat('timeint(0,',num2str(simTimeInSeconds),',Ecell)'),'dataset','dset2','solnum','end')/simTimeInSeconds; %in V
            mean_internal_resistance = mphglobal(model,strcat('timeint(0,',num2str(simTimeInSeconds),',R_cell_int)'),'dataset','dset2','solnum','end')/simTimeInSeconds; %in ohm
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
    cellSimParams.cur_soc_at_load = mphglobal(model,'SOCcell_load','dataset','dset2','solnum','end');
    cellSimParams.cur_soc_coulombic = mphglobal(model,'SOCcell','dataset','dset2','solnum','end');
    cellSimParams.cur_soc = cellSimParams.cur_soc_coulombic;
    cellSimParams.cur_voltage = mphglobal(model,'Ecell','dataset','dset2','solnum','end');
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