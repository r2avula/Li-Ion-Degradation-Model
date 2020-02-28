function [cellSimData_out] = getDegradationData(config)
batteryNominalVoltage = config.batteryNominalVoltage;
cell_nominalVoltage = config.cell_nominalVoltage; %in V
cellsInSeries = ceil(batteryNominalVoltage/cell_nominalVoltage);
if(cellsInSeries~=floor(cellsInSeries))
    warning('Battery voltage is modified!');
end
p_pu = config.powerQuantPU; % in W
slotIntervalInSeconds = config.slotIntervalInSeconds;
slotIntervalInHours = slotIntervalInSeconds/3600; %in h
batteryNominalVoltage = cellsInSeries*cell_nominalVoltage;% in V
converterEfficiency = (config.converterEfficiency)/100;     
batteryRatedCapacityInAh = config.batteryRatedCapacityInAh; %in Ah
cell_SOC_high = config.cell_SOC_high;
cell_SOC_low = config.cell_SOC_low;
z_cap = (batteryRatedCapacityInAh*batteryNominalVoltage); % in Wh

cell_1C_capacityInAh = config.cell_1C_capacityInAh; %in Ah
cell_1C_power = floor(cell_1C_capacityInAh*cell_nominalVoltage); %in W
legsInParallel = round(batteryRatedCapacityInAh/cell_1C_capacityInAh);
d_max_ch_ess = cell_1C_power*legsInParallel*cellsInSeries/converterEfficiency;
d_max_disch_ess = -cell_1C_power*legsInParallel*cellsInSeries*converterEfficiency;

d_rated = config.converterRatedPower;
if(ischar(d_rated))
    d_rated = str2double(d_rated);
end
d_max_ch = min(d_rated,d_max_ch_ess);
d_max_disch = max(-d_rated,d_max_disch_ess);

d_max_ch_pu = floor(d_max_ch/p_pu);
d_max_disch_pu = ceil(d_max_disch/p_pu);
out_pow_set = (d_max_disch_pu:d_max_ch_pu)*p_pu;
pow_num = length(out_pow_set);

bat_pow_set = zeros(1,pow_num);
for pow_idx = 1:pow_num
    if(out_pow_set(pow_idx)<0)
        bat_pow_set(pow_idx) = out_pow_set(pow_idx)/converterEfficiency;
    else
        bat_pow_set(pow_idx) = out_pow_set(pow_idx)*converterEfficiency;
    end
end

e_pu = p_pu*slotIntervalInHours; % in Wh
z_min_pu = floor(cell_SOC_low*z_cap/e_pu);
z_max_pu = floor(cell_SOC_high*z_cap/e_pu);
soc_num = z_max_pu-z_min_pu+1;

z_grid = (z_min_pu:z_max_pu)*e_pu;
soc_grid_boundaries = z_grid/z_cap;
if(soc_grid_boundaries(1)<cell_SOC_low && soc_grid_boundaries(2) >cell_SOC_low)
    soc_grid_boundaries(1) = cell_SOC_low;
end
if(soc_grid_boundaries(end)>=cell_SOC_high)
    soc_grid_boundaries(end)= cell_SOC_high;
else
    soc_grid_boundaries = [soc_grid_boundaries cell_SOC_high];
end

soc_grid_bin_mean = zeros(soc_num,1);
for bin_idx = 1:soc_num
    soc_grid_bin_mean(bin_idx) = (soc_grid_boundaries(bin_idx) + soc_grid_boundaries(bin_idx+1))/2;
end

cell_pow_set = round(bat_pow_set/(cellsInSeries*legsInParallel),config.paramsPrecisionDigits);

deglifePartitions_num =  config.deglifePartitions_num;
deglifePartitions = linspace(1,0.8,deglifePartitions_num+1);
allowedRelativeCapacityChange = (deglifePartitions(1)-deglifePartitions(2))*100;

cellSimData_out = cell(deglifePartitions_num,1);

totalDegSampleNum = config.degSampleNum; 
degSampleNum_per_sim = 20; % to avoid COMSOL process abruptly killed by Linux kernel for taking long memeory
simRandRingersNum = floor(totalDegSampleNum/degSampleNum_per_sim);
for partition_idx = 1:deglifePartitions_num
    cellSimParams = struct;
    cellSimParams.soc_grid_boundaries = soc_grid_boundaries;
    cellSimParams.cell_pow_set = cell_pow_set;
    cellSimParams.initialRelCap = deglifePartitions(partition_idx)*100;
    cellSimParams.allowedRelativeCapacityChange = allowedRelativeCapacityChange; % 20
    cellSimParams.sample_num = degSampleNum_per_sim; 
    cellSimParams.slotIntervalInSeconds = slotIntervalInSeconds;
    cellSimParams.SOC_low = config.cell_SOC_low;
    cellSimParams.SOC_high = config.cell_SOC_high;
    cellSimParams.SOC_init = config.cell_SOC_init;
    cellSimParams.SOC_init = config.cell_SOC_init;
    cellSimParams.cell_voltage_high = config.cell_voltage_high;
    cellSimParams.cell_voltage_low = config.cell_voltage_low;
    cellSimParams.driveToSOH_timeAccelerationFactor = config.driveToSOH_timeAccelerationFactor;
    cellSimParams.driveToSOC_timeAccelerationFactor = config.driveToSOC_timeAccelerationFactor;
    cellSimParams.driveToSOC_timePeriodScaleFactor = config.driveToSOC_timePeriodScaleFactor;
    driveToSOC_attempts_max = config.driveToSOC_attempts_max;
    if(ischar(driveToSOC_attempts_max))
        driveToSOC_attempts_max = str2double(driveToSOC_attempts_max);
    end
    cellSimParams.driveToSOC_attempts_max = driveToSOC_attempts_max;
    fileNamePrefix = ['cache' filesep 'cellSimData_' num2str(batteryRatedCapacityInAh) '_'];   
    cellSimData_cur_partition = cell(simRandRingersNum,1);
    for simRandRinger = 1:simRandRingersNum
        cellSimParams.simRandRinger = simRandRinger;
        [filename,fileExists] = findFileName(cellSimParams,fileNamePrefix,'cellSimParams');
        if(fileExists)
            load(filename,'cellSimData');
            disp(strcat({'cellSimData loaded from '},filename,' .'));
        else
            rng(simRandRinger,'twister');
            [cellSimData] = runComsolSimulation(cellSimParams,config.comsol_model_filename);
            save(filename,'cellSimData','cellSimParams')
            disp(strcat({'cellSimData saved to '},filename,' .'));
        end
        cellSimData_cur_partition{simRandRinger} = cellSimData;
    end    
    cellSimData_out{partition_idx} = addCalenderAgingAndAaggregateSamples(cellSimData_cur_partition,cellSimParams,soc_grid_bin_mean);
end
end

