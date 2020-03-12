function [mean_nrmse_energyLoss_percentage,mean_nrmse_soc_kp1_percentage] = plotMeanDegradationMaps_1hour(config,cellSimData_all_partitions)
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
d_num = length(out_pow_set);

bat_pow_set = zeros(1,d_num);
for pow_idx = 1:d_num
    if(out_pow_set(pow_idx)<0)
        bat_pow_set(pow_idx) = out_pow_set(pow_idx)/converterEfficiency;
    else
        bat_pow_set(pow_idx) = out_pow_set(pow_idx)*converterEfficiency;
    end
end

e_pu = p_pu*slotIntervalInHours; % in Wh
z_min_pu = floor(cell_SOC_low*z_cap/e_pu);
z_max_pu = floor(cell_SOC_high*z_cap/e_pu);
z_num = z_max_pu-z_min_pu+1;

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

soc_grid_bin_mean = zeros(z_num,1);
for bin_idx = 1:z_num
    soc_grid_bin_mean(bin_idx) = (soc_grid_boundaries(bin_idx) + soc_grid_boundaries(bin_idx+1))/2;
end

cell_pow_set = round(bat_pow_set/(cellsInSeries*legsInParallel),config.paramsPrecisionDigits);

deglifePartitions_num =  config.deglifePartitions_num;

[soc_num_sim,pow_num_sim,~] = size(cellSimData_all_partitions{2}.simTimeRatio_samples);
if(d_num ~= pow_num_sim || z_num ~= soc_num_sim)
    error('Something is wrong!');
end

cell_cur_set_sim = cell_pow_set/cell_nominalVoltage; % in A

batterySelfDischargeRatePerMonth = 0; % factor in [0,1]
tau = 30*24/-log(1-batterySelfDischargeRatePerMonth); %h
gamma = 1 - exp(-slotIntervalInHours/tau); % factor in [0,1]
if(gamma >0)
    gamma_tau = gamma*tau; % in h
else
    gamma_tau = slotIntervalInHours; % in h
end

mean_capacityLoss_factor_comsol = nan*zeros(soc_num_sim,pow_num_sim,deglifePartitions_num);
mean_capacityLoss_factor_incl_calender_ageing = nan*zeros(soc_num_sim,pow_num_sim,deglifePartitions_num);
mean_resEnergyLoss_percentage_3c = nan*zeros(soc_num_sim,pow_num_sim,deglifePartitions_num);
mean_resEnergyLoss_percentage_comsol = nan*zeros(soc_num_sim,pow_num_sim,deglifePartitions_num);
mean_conEnergyLoss_percentage_3c = nan*zeros(soc_num_sim,pow_num_sim,deglifePartitions_num);
mean_energyLoss_percentage_3c = nan*zeros(soc_num_sim,pow_num_sim,deglifePartitions_num);

ess_mean_energyLossInWh_map = nan*zeros(soc_num_sim,pow_num_sim,deglifePartitions_num);

z_kp1_idxs_prob_map = zeros(z_num,z_num,d_num,deglifePartitions_num);
mean_soc_kp1_3c = nan*zeros(z_num,d_num,deglifePartitions_num); %after quantization
relVar_soc_kp1_3c = nan*zeros(z_num,d_num,deglifePartitions_num); %after quantization


nrmse_energyLoss_percentage = nan*zeros(soc_num_sim,pow_num_sim,deglifePartitions_num);
nrmse_soc_kp1_percentage = nan*zeros(soc_num_sim,pow_num_sim,deglifePartitions_num); 

for partition_idx = 1:deglifePartitions_num
    cellSimData = cellSimData_all_partitions{partition_idx};
    relative_capacity_samples = cellSimData.relative_capacity_samples;
    cell_energy_applied_samples_comsol = cellSimData.energy_applied_samples;
    cell_energy_loss_samples_comsol = cellSimData.energy_loss_samples;
    cell_capacity_loss_factor_samples = cellSimData.capacity_loss_factor_samples;
    capacity_loss_factor_incl_calender_samples = cellSimData.capacity_loss_factor_incl_calender_samples;
    cell_internal_resistance_samples = cellSimData.internal_resistance_samples; % in Ohm
    cell_terminal_voltage_samples = cellSimData.terminal_voltage_samples; % in V
    simTimeRatio_samples = cellSimData.simTimeRatio_samples;
    soc_k_samples_3c = cellSimData.soc_k_samples_3c;
    soc_kp1_samples_3c = cellSimData.soc_kp1_samples_3c;
    soc_k_samples_coulombic = cellSimData.soc_k_samples_coulombic;
    soc_kp1_samples_coulombic = cellSimData.soc_kp1_samples_coulombic;
    
    simTimeThreshold = 0.5;
    
    for z_k_idx = 1:soc_num_sim
        for d_k_idx = 1:pow_num_sim
            simTimeRatioSamples_vec = reshape(simTimeRatio_samples(z_k_idx,d_k_idx,:),1,[]);
                        
            cell_capacity_loss_factor_samples_vec = reshape(cell_capacity_loss_factor_samples(z_k_idx,d_k_idx,:),1,[]);
            capacity_loss_factor_incl_calender_samples_vec = reshape(capacity_loss_factor_incl_calender_samples(z_k_idx,d_k_idx,:),1,[]);
            cell_resistance_samples_vec = reshape(cell_internal_resistance_samples(z_k_idx,d_k_idx,:),1,[]);
            cell_voltage_samples_vec = reshape(cell_terminal_voltage_samples(z_k_idx,d_k_idx,:),1,[]);
            ess_energyLoss_comsol_vec = reshape(cell_energy_loss_samples_comsol(z_k_idx,d_k_idx,:),1,[])*legsInParallel*cellsInSeries;
            relative_capacity_factor_samples_vec = reshape(relative_capacity_samples(z_k_idx,d_k_idx,:),1,[])/100;
            ess_bat_energy_applied_samples_comsol_vec = reshape(cell_energy_applied_samples_comsol(z_k_idx,d_k_idx,:),1,[])*legsInParallel*cellsInSeries;
            soc_k_samples_3c_vec = reshape(soc_k_samples_3c(z_k_idx,d_k_idx,:),1,[]);
            soc_kp1_samples_3c_vec = reshape(soc_kp1_samples_3c(z_k_idx,d_k_idx,:),1,[]);
            soc_k_samples_coulombic_vec = reshape(soc_k_samples_coulombic(z_k_idx,d_k_idx,:),1,[]);
            soc_kp1_samples_coulombic_vec = reshape(soc_kp1_samples_coulombic(z_k_idx,d_k_idx,:),1,[]);
            
            validSampleIdxs = (~isnan(cell_capacity_loss_factor_samples_vec)&~isnan(cell_resistance_samples_vec)&...
                ~isnan(cell_voltage_samples_vec)&~isnan(soc_k_samples_coulombic_vec)&simTimeRatioSamples_vec>=simTimeThreshold);
            
            cell_capacity_loss_factor_samples_vec = cell_capacity_loss_factor_samples_vec(validSampleIdxs);
            capacity_loss_factor_incl_calender_samples_vec = capacity_loss_factor_incl_calender_samples_vec(validSampleIdxs);
            cell_resistance_samples_vec = cell_resistance_samples_vec(validSampleIdxs);
            cell_voltage_samples_vec = cell_voltage_samples_vec(validSampleIdxs);
            ess_energyLoss_comsol_vec = abs(ess_energyLoss_comsol_vec(validSampleIdxs));
            relative_capacity_factor_samples_vec = relative_capacity_factor_samples_vec(validSampleIdxs);
            ess_bat_energy_applied_samples_comsol_vec = ess_bat_energy_applied_samples_comsol_vec(validSampleIdxs);
            soc_k_samples_3c_vec = soc_k_samples_3c_vec(validSampleIdxs);
            soc_kp1_samples_3c_vec = soc_kp1_samples_3c_vec(validSampleIdxs);
            soc_k_samples_coulombic_vec = soc_k_samples_coulombic_vec(validSampleIdxs);
            soc_kp1_samples_coulombic_vec = soc_kp1_samples_coulombic_vec(validSampleIdxs);
            simTimeRatioSamples_vec = simTimeRatioSamples_vec(validSampleIdxs);
            numValidSamples = length(cell_voltage_samples_vec);
                        
            if(numValidSamples>0)                
                mean_capacityLoss_factor_comsol(z_k_idx,d_k_idx,partition_idx) = abs(mean(cell_capacity_loss_factor_samples_vec.*simTimeRatioSamples_vec)); % actual data
                mean_capacityLoss_factor_incl_calender_ageing(z_k_idx,d_k_idx,partition_idx) = abs(mean(capacity_loss_factor_incl_calender_samples_vec.*simTimeRatioSamples_vec)); 
                cell_res_energyLossEstimate_3c = zeros(1,numValidSamples);
                cell_con_energyLossEstimate_3c = zeros(1,numValidSamples);
                z_kp1_estimate_3c = zeros(1,numValidSamples);
                z_kp1_comsol = zeros(1,numValidSamples);
                
                for validSampleIdx = 1:numValidSamples
                    cell_cur_en_cap = relative_capacity_factor_samples_vec(validSampleIdx)*cell_1C_capacityInAh*cell_nominalVoltage;
                    cell_cur_volt = cell_voltage_samples_vec(validSampleIdx);
                    cell_cur_pow = cell_cur_set_sim(d_k_idx)*cell_cur_volt;
                    cell_cur_res = cell_resistance_samples_vec(validSampleIdx);
                    if(cell_cur_pow>0)
                        cell_pow_in_w_conv = cell_cur_pow/converterEfficiency;
                    else
                        cell_pow_in_w_conv = cell_cur_pow*converterEfficiency;
                    end
                    cur_soc = soc_k_samples_3c_vec(validSampleIdx);
                    cur_soc_bin_idx = find(soc_grid_boundaries(2:end-1)>cur_soc,1);
                    if(isempty(cur_soc_bin_idx))
                        cur_soc_bin_idx =z_num;
                    end                    
                    if(cur_soc_bin_idx~=z_k_idx)
                        error('something is wrong!');
                    end
                    
                    cell_res_energyLossEstimate_3c(validSampleIdx) = gamma*cur_soc*cell_cur_en_cap +...
                        simTimeRatioSamples_vec(validSampleIdx)*(cell_cur_pow*slotIntervalInHours - ...
                        (gamma_tau*cell_cur_volt/2/cell_cur_res)*(sqrt(max((cell_cur_volt*cell_cur_volt)+4*cell_cur_res*cell_cur_pow,0))-cell_cur_volt));
                    
                    cell_con_energyLossEstimate_3c(validSampleIdx) = simTimeRatioSamples_vec(validSampleIdx)*((cell_pow_in_w_conv-cell_cur_pow)*slotIntervalInHours);
                                        
                    next_soc_bin_idx = find(soc_grid_boundaries(2:end-1)>soc_kp1_samples_3c_vec(validSampleIdx),1);
                    if(isempty(next_soc_bin_idx))
                        next_soc_bin_idx =z_num;
                    end
                    z_kp1_estimate_3c(validSampleIdx) = next_soc_bin_idx;           
                    
                    next_soc_bin_idx = find(soc_grid_boundaries(2:end-1)>soc_kp1_samples_coulombic_vec(validSampleIdx),1);
                    if(isempty(next_soc_bin_idx))
                        next_soc_bin_idx =z_num;
                    end
                    z_kp1_comsol(validSampleIdx) = next_soc_bin_idx;    
                end
                
                ess_res_energyLossEstimate_3c = max(cell_res_energyLossEstimate_3c*legsInParallel*cellsInSeries,0);
                ess_con_energyLossEstimate_3c = max(cell_con_energyLossEstimate_3c*legsInParallel*cellsInSeries,0);
                ess_energyLossEstimate_3c = ess_res_energyLossEstimate_3c + ess_con_energyLossEstimate_3c;
                ess_energyApplied = ess_bat_energy_applied_samples_comsol_vec +ess_con_energyLossEstimate_3c;
                
                ess_mean_energyLossInWh_map(z_k_idx,d_k_idx,partition_idx) = mean(abs(ess_energyLossEstimate_3c));
                                
                for idx = 1:z_num
                    z_kp1_idxs_prob_map(idx,z_k_idx,d_k_idx,partition_idx) = sum(z_kp1_estimate_3c==idx)/numValidSamples;
                end
                 
                nrmse_soc_kp1_percentage(z_k_idx,d_k_idx,partition_idx) = sqrt(mean(((soc_grid_bin_mean(z_kp1_comsol) - soc_grid_bin_mean(z_kp1_estimate_3c))).^2))/mean(soc_grid_bin_mean(z_kp1_comsol))*100; % including quantization   
%                 nrmse_soc_kp1_percentage(z_k_idx,d_k_idx,partition_idx) = sqrt(mean(((soc_kp1_samples_vec - soc_kp1_estimate_3c)).^2))/mean(soc_kp1_samples_vec)*100; % without quantization   
                                
                mean_soc_kp1_3c(z_k_idx,d_k_idx,partition_idx) = soc_grid_bin_mean'*z_kp1_idxs_prob_map(:,z_k_idx,d_k_idx,partition_idx);% after quantization
                soc_kp1_3c_var = (soc_grid_bin_mean.*soc_grid_bin_mean)'*z_kp1_idxs_prob_map(:,z_k_idx,d_k_idx,partition_idx) - ...
                    mean_soc_kp1_3c(z_k_idx,d_k_idx,partition_idx)^2;% after quantization
                relVar_soc_kp1_3c(z_k_idx,d_k_idx,partition_idx) = soc_kp1_3c_var/mean_soc_kp1_3c(z_k_idx,d_k_idx,partition_idx)*100;
                
                if(abs(cell_cur_set_sim(d_k_idx)) > 0)                    
                    mean_resEnergyLoss_percentage_3c(z_k_idx,d_k_idx,partition_idx) = mean(abs(ess_res_energyLossEstimate_3c./ess_energyApplied))*100;   
                    mean_conEnergyLoss_percentage_3c(z_k_idx,d_k_idx,partition_idx) = mean(abs(ess_con_energyLossEstimate_3c./ess_energyApplied))*100;
                    mean_energyLoss_percentage_3c(z_k_idx,d_k_idx,partition_idx) = mean(abs(ess_energyLossEstimate_3c./ess_energyApplied))*100;
                                                            
                    mean_resEnergyLoss_percentage_comsol(z_k_idx,d_k_idx,partition_idx) = mean(abs(ess_energyLoss_comsol_vec./ess_energyApplied))*100;
                    
                    nrmse_energyLoss_percentage(z_k_idx,d_k_idx,partition_idx) =...
                        sqrt(mean(((ess_res_energyLossEstimate_3c - ess_energyLoss_comsol_vec)).^2))/mean(ess_energyLoss_comsol_vec)*100;                    
                else
                    mean_resEnergyLoss_percentage_3c(z_k_idx,d_k_idx,partition_idx) = 0;
                    mean_conEnergyLoss_percentage_3c(z_k_idx,d_k_idx,partition_idx) = 0;
                    mean_energyLoss_percentage_3c(z_k_idx,d_k_idx,partition_idx) = 0;
                    
                    mean_resEnergyLoss_percentage_comsol(z_k_idx,d_k_idx,partition_idx) = 0;
                    
                    nrmse_energyLoss_percentage(z_k_idx,d_k_idx,partition_idx) = 0;
                end
            end
        end
    end
end

nrmse_energyLoss_percentage_vec = nrmse_energyLoss_percentage(:);
nrmse_energyLoss_percentage_vec = nrmse_energyLoss_percentage_vec(~isnan(nrmse_energyLoss_percentage_vec));
mean_nrmse_energyLoss_percentage = mean(nrmse_energyLoss_percentage_vec);

nrmse_soc_kp1_percentage_vec = nrmse_soc_kp1_percentage(:);
nrmse_soc_kp1_percentage_vec = nrmse_soc_kp1_percentage_vec(~isnan(nrmse_soc_kp1_percentage_vec));
mean_nrmse_soc_kp1_percentage = mean(nrmse_soc_kp1_percentage_vec);

%% Plot maps
savefiles = 0;
show_plots = 1;

plot_life_partitions_res_energy_loss_percentage_3c = 1;
if(show_plots&&plot_life_partitions_res_energy_loss_percentage_3c)
    filename = 'resEnLoss';
    partition_idxs = 1:deglifePartitions_num;
    positions = {[0.0924678097197633 0.246165556113454 0.20936667244291 0.615405872457975],...
        [0.373264911169039 0.238165556113454 0.20936667244291 0.623405872457975],...
        [0.654062012618313 0.238165556113454 0.209366672442911 0.627405872457975]};
    
    textpositions = {[0.232218011771345 1.10249990460461 0],...
        [0.181262597758605 1.10233123498622 0],...
        [0.24495686527453 1.09583826559003 0]...
        };
    
    y_axis_grid = cell_pow_set/(cell_1C_capacityInAh*cell_nominalVoltage);
    x_axis_grid = unique([0,soc_grid_boundaries,1]);
    x_axis_grid_mean = zeros(length(x_axis_grid)-1,1);
    for idx = 1:length(x_axis_grid)-1
        x_axis_grid_mean(idx) = 0.5*(x_axis_grid(idx) +x_axis_grid(idx+1));
    end
    
    [x,y] = meshgrid(x_axis_grid,y_axis_grid);
    
    plot_data = nan(length(y_axis_grid),length(x_axis_grid),deglifePartitions_num);
    for partition_idx = 1:deglifePartitions_num
        for soc_bin_idx = 1:z_num
            plot_data(:,soc_bin_idx+1,partition_idx) = mean_resEnergyLoss_percentage_3c(soc_bin_idx,:,partition_idx);
        end
    end
    plot_data(:,z_num+2,:) = plot_data(:,z_num+1,:);
    
    data_min = min(plot_data(:));
    data_max = max(plot_data(:));
    
    fontSize = 12;
    figure_xsize = 700;
    figure_ysize = 200;
    figure_position_offset = 50;
    figure_1 = figure('Color','w','Renderer','painters');
    set(gcf, 'Position',  [figure_position_offset,figure_position_offset,figure_position_offset+figure_xsize, figure_position_offset+figure_ysize]);
    
    for partition_idx = partition_idxs        
        axes1 = axes('Parent',figure_1,...
            'Position',positions{partition_idx});
        box(axes1,'on');
        hold(axes1,'on');
        h1 = surf(x,y,0*plot_data(:,:,partition_idx));
        xline(cell_SOC_low,'Alpha',1,'Color',[1 0 0],'LineStyle','--',...
            'LineWidth',1);
        xline(cell_SOC_high,'Alpha',1,'Color',[1 0 0],'LineStyle','--',...
            'LineWidth',1);
        ylabel({'Current (A)'},'Interpreter','latex','Units', 'normalized','position',[-0.113092878602192,0.525806928450061,0]);
        xlabel({'SOC'},'Interpreter','latex','Units','normalized');
        set(h1,'CData',plot_data(:,:,partition_idx),'FaceColor','interp','FaceAlpha',1)
        shading interp;
        axes1.YAxis.MinorTickValues = -1:0.5:1;
        xlim([0,1]);
        ylim([-1 ,1]);
        set(axes1,'FontSize',fontSize,'TickLabelInterpreter','latex',...
            'XMinorTick','on','YMinorTick','on','XTick',0:0.25:1,'YTick',[-1  0 1],...
            'YTickLabel',{'-1C','0','1C'});
        if(partition_idx == 1)
            text('Units','normalized','BackgroundColor',[1 1 1],...
                'EdgeColor',[0 0 0],...
                'Interpreter','latex',...
                'String',strcat('age $\in$[0,1/',num2str(deglifePartitions_num),']'),...
                'Position',textpositions{partition_idx},...
                'Color',[0 0 0],'FontSize',fontSize);
        elseif(partition_idx == deglifePartitions_num)
            text('Units','normalized','BackgroundColor',[1 1 1],...
                'EdgeColor',[0 0 0],...
                'Interpreter','latex',...
                'String',strcat('age $\in$[',num2str(deglifePartitions_num-1),'/',num2str(deglifePartitions_num),',1]'),...
                'Position',textpositions{partition_idx},...
                'Color',[0 0 0],'FontSize',fontSize);
        else
            text('Units','normalized','BackgroundColor',[1 1 1],...
                'EdgeColor',[0 0 0],...
                'Interpreter','latex',...
                'String',strcat('age $\in$[',num2str(partition_idx-1),'/',num2str(deglifePartitions_num),',',num2str(partition_idx),'/',num2str(deglifePartitions_num),']'),...
                'Position',textpositions{partition_idx},...
                'Color',[0 0 0],'FontSize',fontSize);
        end
        colormap parula;
        caxis([data_min data_max]);
    end
    labelTicks = round(linspace(data_min,data_max,5),2);
    h = colorbar('Position',[0.888 0.068 0.0266666666666669 0.892],...
        'TickLabelInterpreter','latex','FontSize',fontSize,'Ticks',labelTicks,'Limits',[0,round(data_max,2)]);
    caxis([data_min data_max]);
    ylabel(h, 'Rel. resistive en. loss estimate ($\%$)','Interpreter','latex','FontSize',fontSize);
    colormap parula;
    fig_pos = figure_1.PaperPosition;
    figure_1.PaperSize = [fig_pos(3) fig_pos(4)];
    if(savefiles)
        saveas(gcf,strcat(filename,'.fig'));
        saveas(gcf,strcat(filename,'.eps'),'epsc');
    end
end

plot_life_partitions_con_energy_loss_percentage_3c = 0;
if(plot_life_partitions_con_energy_loss_percentage_3c)
    filename = 'converterLoss';
    partition_idxs = 1:deglifePartitions_num;
    positions = {[0.0924678097197633 0.246165556113454 0.20936667244291 0.615405872457975],...
        [0.373264911169039 0.238165556113454 0.20936667244291 0.623405872457975],...
        [0.654062012618313 0.238165556113454 0.209366672442911 0.627405872457975]};
    
    textpositions = {[0.232218011771345 1.10249990460461 0],...
        [0.181262597758605 1.10233123498622 0],...
        [0.24495686527453 1.09583826559003 0]...
        };
    
    y_axis_grid = cell_pow_set/(cell_1C_capacityInAh*cell_nominalVoltage);
    x_axis_grid = unique([0,soc_grid_boundaries,1]);
    x_axis_grid_mean = zeros(length(x_axis_grid)-1,1);
    for idx = 1:length(x_axis_grid)-1
        x_axis_grid_mean(idx) = 0.5*(x_axis_grid(idx) +x_axis_grid(idx+1));
    end
    
    [x,y] = meshgrid(x_axis_grid,y_axis_grid);
    
    plot_data = nan(length(y_axis_grid),length(x_axis_grid),deglifePartitions_num);
    for partition_idx = 1:deglifePartitions_num
        for soc_bin_idx = 1:z_num
            plot_data(:,soc_bin_idx+1,partition_idx) = mean_conEnergyLoss_percentage_3c(soc_bin_idx,:,partition_idx);
        end
    end
    plot_data(:,z_num+2,:) = plot_data(:,z_num+1,:);
    
    data_min = min(plot_data(:));
    data_max = max(plot_data(:));
    
    fontSize = 12;
    figure_xsize = 700;
    figure_ysize = 200;
    figure_position_offset = 50;
    figure_1 = figure('Color','w','Renderer','painters');
    set(gcf, 'Position',  [figure_position_offset,figure_position_offset,figure_position_offset+figure_xsize, figure_position_offset+figure_ysize]);
    
    for partition_idx = partition_idxs
        axes1 = axes('Parent',figure_1,...
            'Position',positions{partition_idx});
        box(axes1,'on');
        hold(axes1,'on');
        h1 = surf(x,y,plot_data(:,:,partition_idx));
        xline(cell_SOC_low,'Alpha',1,'Color',[1 0 0],'LineStyle','--',...
            'LineWidth',1);
        xline(cell_SOC_high,'Alpha',1,'Color',[1 0 0],'LineStyle','--',...
            'LineWidth',1);
        ylabel({'Current (A)'},'Interpreter','latex','Units', 'normalized','position',[-0.113092878602192,0.525806928450061,0]);
        xlabel({'SOC'},'Interpreter','latex','Units','normalized');
        set(h1, 'FaceColor','interp','FaceAlpha',1)
        %         shading interp;
        axes1.YAxis.MinorTickValues = -1:0.5:1;
        xlim([0,1]);
        ylim([-1 ,1]);
        set(axes1,'FontSize',fontSize,'TickLabelInterpreter','latex',...
            'XMinorTick','on','YMinorTick','on','XTick',0:0.25:1,'YTick',[-1  0 1],...
            'YTickLabel',{'-1C','0','1C'});
        if(partition_idx == 1)
            text('Units','normalized','BackgroundColor',[1 1 1],...
                'EdgeColor',[0 0 0],...
                'Interpreter','latex',...
                'String',strcat('age $\in$[0,1/',num2str(deglifePartitions_num),']'),...
                'Position',textpositions{partition_idx},...
                'Color',[0 0 0],'FontSize',fontSize);
        elseif(partition_idx == deglifePartitions_num)
            text('Units','normalized','BackgroundColor',[1 1 1],...
                'EdgeColor',[0 0 0],...
                'Interpreter','latex',...
                'String',strcat('age $\in$[',num2str(deglifePartitions_num-1),'/',num2str(deglifePartitions_num),',1]'),...
                'Position',textpositions{partition_idx},...
                'Color',[0 0 0],'FontSize',fontSize);
        else
            text('Units','normalized','BackgroundColor',[1 1 1],...
                'EdgeColor',[0 0 0],...
                'Interpreter','latex',...
                'String',strcat('age $\in$[',num2str(partition_idx-1),'/',num2str(deglifePartitions_num),',',num2str(partition_idx),'/',num2str(deglifePartitions_num),']'),...
                'Position',textpositions{partition_idx},...
                'Color',[0 0 0],'FontSize',fontSize);
        end
        colormap jet;
        caxis([data_min data_max]);
    end
    labelTicks = round(linspace(data_min,data_max,5),2);
    h = colorbar('Position',[0.888 0.068 0.0266666666666669 0.892],...
        'TickLabelInterpreter','latex','FontSize',fontSize,'Ticks',labelTicks,'Limits',[0,round(data_max,2)]);
    caxis([data_min data_max]);
    ylabel(h, 'Rel. con. en. loss estimate ($\%$)','Interpreter','latex','FontSize',fontSize);
    colormap jet;
    fig_pos = figure_1.PaperPosition;
    figure_1.PaperSize = [fig_pos(3) fig_pos(4)];
    if(savefiles)
        saveas(gcf,strcat(filename,'.fig'));
        saveas(gcf,strcat(filename,'.eps'),'epsc');
    end
end

plot_life_partitions_energy_loss_percentage_3c = 0;
if(plot_life_partitions_energy_loss_percentage_3c)
    filename = 'totalEnergyLoss';
    partition_idxs = 1:deglifePartitions_num;
    positions = {[0.0924678097197633 0.246165556113454 0.20936667244291 0.615405872457975],...
        [0.373264911169039 0.238165556113454 0.20936667244291 0.623405872457975],...
        [0.654062012618313 0.238165556113454 0.209366672442911 0.627405872457975]};
    
    textpositions = {[0.232218011771345 1.10249990460461 0],...
        [0.181262597758605 1.10233123498622 0],...
        [0.24495686527453 1.09583826559003 0]...
        };
    
    y_axis_grid = cell_pow_set/(cell_1C_capacityInAh*cell_nominalVoltage);
    x_axis_grid = unique([0,soc_grid_boundaries,1]);
    x_axis_grid_mean = zeros(length(x_axis_grid)-1,1);
    for idx = 1:length(x_axis_grid)-1
        x_axis_grid_mean(idx) = 0.5*(x_axis_grid(idx) +x_axis_grid(idx+1));
    end
    
    [x,y] = meshgrid(x_axis_grid,y_axis_grid);
    
    plot_data = nan(length(y_axis_grid),length(x_axis_grid),deglifePartitions_num);
    for partition_idx = 1:deglifePartitions_num
        for soc_bin_idx = 1:z_num
            plot_data(:,soc_bin_idx+1,partition_idx) = mean_energyLoss_percentage_3c(soc_bin_idx,:,partition_idx);
        end
    end
    plot_data(:,z_num+2,:) = plot_data(:,z_num+1,:);
    
    data_min = min(plot_data(:));
    data_max = max(plot_data(:));
    
    fontSize = 12;
    figure_xsize = 700;
    figure_ysize = 200;
    figure_position_offset = 50;
    figure_1 = figure('Color','w','Renderer','painters');
    set(gcf, 'Position',  [figure_position_offset,figure_position_offset,figure_position_offset+figure_xsize, figure_position_offset+figure_ysize]);
    
    for partition_idx = partition_idxs
        axes1 = axes('Parent',figure_1,...
            'Position',positions{partition_idx});
        box(axes1,'on');
        hold(axes1,'on');
        h1 = surf(x,y,plot_data(:,:,partition_idx));
        xline(cell_SOC_low,'Alpha',1,'Color',[1 0 0],'LineStyle','--',...
            'LineWidth',1);
        xline(cell_SOC_high,'Alpha',1,'Color',[1 0 0],'LineStyle','--',...
            'LineWidth',1);
        ylabel({'Current (A)'},'Interpreter','latex','Units', 'normalized','position',[-0.113092878602192,0.525806928450061,0]);
        xlabel({'SOC'},'Interpreter','latex','Units','normalized');
        set(h1, 'FaceColor','interp','FaceAlpha',1)
        %         shading interp;
        axes1.YAxis.MinorTickValues = -1:0.5:1;
        xlim([0,1]);
        ylim([-1 ,1]);
        set(axes1,'FontSize',fontSize,'TickLabelInterpreter','latex',...
            'XMinorTick','on','YMinorTick','on','XTick',0:0.25:1,'YTick',[-1  0 1],...
            'YTickLabel',{'-1C','0','1C'});
        if(partition_idx == 1)
            text('Units','normalized','BackgroundColor',[1 1 1],...
                'EdgeColor',[0 0 0],...
                'Interpreter','latex',...
                'String',strcat('age $\in$[0,1/',num2str(deglifePartitions_num),']'),...
                'Position',textpositions{partition_idx},...
                'Color',[0 0 0],'FontSize',fontSize);
        elseif(partition_idx == deglifePartitions_num)
            text('Units','normalized','BackgroundColor',[1 1 1],...
                'EdgeColor',[0 0 0],...
                'Interpreter','latex',...
                'String',strcat('age $\in$[',num2str(deglifePartitions_num-1),'/',num2str(deglifePartitions_num),',1]'),...
                'Position',textpositions{partition_idx},...
                'Color',[0 0 0],'FontSize',fontSize);
        else
            text('Units','normalized','BackgroundColor',[1 1 1],...
                'EdgeColor',[0 0 0],...
                'Interpreter','latex',...
                'String',strcat('age $\in$[',num2str(partition_idx-1),'/',num2str(deglifePartitions_num),',',num2str(partition_idx),'/',num2str(deglifePartitions_num),']'),...
                'Position',textpositions{partition_idx},...
                'Color',[0 0 0],'FontSize',fontSize);
        end
        colormap jet;
        caxis([data_min data_max]);
    end
    labelTicks = round(linspace(data_min,data_max,5),2);
    h = colorbar('Position',[0.888 0.068 0.0266666666666669 0.892],...
        'TickLabelInterpreter','latex','FontSize',fontSize,'Ticks',labelTicks,'Limits',[0,round(data_max,2)]);
    caxis([data_min data_max]);
    ylabel(h, 'Rel. en. loss estimate ($\%$)','Interpreter','latex','FontSize',fontSize);
    colormap jet;
    fig_pos = figure_1.PaperPosition;
    figure_1.PaperSize = [fig_pos(3) fig_pos(4)];
    if(savefiles)
        saveas(gcf,strcat(filename,'.fig'));
        saveas(gcf,strcat(filename,'.eps'),'epsc');
    end
end

plot_life_partitions_res_energy_loss_percentage_comsol = 0;
if(plot_life_partitions_res_energy_loss_percentage_comsol)
    filename = 'energyDegradation';
    partition_idxs = 1:deglifePartitions_num;
    positions = {[0.0924678097197633 0.246165556113454 0.20936667244291 0.615405872457975],...
        [0.373264911169039 0.238165556113454 0.20936667244291 0.623405872457975],...
        [0.654062012618313 0.238165556113454 0.209366672442911 0.627405872457975]};
    
    textpositions = {[0.232218011771345 1.10249990460461 0],...
        [0.181262597758605 1.10233123498622 0],...
        [0.24495686527453 1.09583826559003 0]...
        };
    
    y_axis_grid = cell_pow_set/(cell_1C_capacityInAh*cell_nominalVoltage);
    x_axis_grid = unique([0,soc_grid_boundaries,1]);
    x_axis_grid_mean = zeros(length(x_axis_grid)-1,1);
    for idx = 1:length(x_axis_grid)-1
        x_axis_grid_mean(idx) = 0.5*(x_axis_grid(idx) +x_axis_grid(idx+1));
    end
    
    [x,y] = meshgrid(x_axis_grid,y_axis_grid);
    
    plot_data = nan(length(y_axis_grid),length(x_axis_grid),deglifePartitions_num);
    for partition_idx = 1:deglifePartitions_num
        for soc_bin_idx = 1:z_num
            plot_data(:,soc_bin_idx+1,partition_idx) = mean_resEnergyLoss_percentage_comsol(soc_bin_idx,:,partition_idx);
        end
    end
    plot_data(:,z_num+2,:) = mean_resEnergyLoss_percentage_comsol(z_num,:,:);
    
    data_min = min(mean_resEnergyLoss_percentage_comsol(:));
    data_max = max(mean_resEnergyLoss_percentage_comsol(:));
    
    fontSize = 12;
    figure_xsize = 700;
    figure_ysize = 200;
    figure_position_offset = 50;
    figure_1 = figure('Color','w','Renderer','painters');
    set(gcf, 'Position',  [figure_position_offset,figure_position_offset,figure_position_offset+figure_xsize, figure_position_offset+figure_ysize]);
    
    for partition_idx = partition_idxs
        axes1 = axes('Parent',figure_1,...
            'Position',positions{partition_idx});
        box(axes1,'on');
        hold(axes1,'on');
        h1 = surf(x,y,plot_data(:,:,partition_idx));
        xline(cell_SOC_low,'Alpha',1,'Color',[1 0 0],'LineStyle','--',...
            'LineWidth',1);
        xline(cell_SOC_high,'Alpha',1,'Color',[1 0 0],'LineStyle','--',...
            'LineWidth',1);
        ylabel({'Current (A)'},'Interpreter','latex','Units', 'normalized','position',[-0.113092878602192,0.525806928450061,0]);
        xlabel({'SOC'},'Interpreter','latex','Units','normalized');
        set(h1, 'FaceColor','interp','FaceAlpha',1)
        %         shading interp;
        axes1.YAxis.MinorTickValues = -1:0.5:1;
        xlim([0,1]);
        ylim([-1 ,1]);
        set(axes1,'FontSize',fontSize,'TickLabelInterpreter','latex',...
            'XMinorTick','on','YMinorTick','on','XTick',0:0.25:1,'YTick',[-1  0 1],...
            'YTickLabel',{'-1C','0','1C'});
        if(partition_idx == 1)
            text('Units','normalized','BackgroundColor',[1 1 1],...
                'EdgeColor',[0 0 0],...
                'Interpreter','latex',...
                'String',strcat('age $\in$[0,1/',num2str(deglifePartitions_num),']'),...
                'Position',textpositions{partition_idx},...
                'Color',[0 0 0],'FontSize',fontSize);
        elseif(partition_idx == deglifePartitions_num)
            text('Units','normalized','BackgroundColor',[1 1 1],...
                'EdgeColor',[0 0 0],...
                'Interpreter','latex',...
                'String',strcat('age $\in$[',num2str(deglifePartitions_num-1),'/',num2str(deglifePartitions_num),',1]'),...
                'Position',textpositions{partition_idx},...
                'Color',[0 0 0],'FontSize',fontSize);
        else
            text('Units','normalized','BackgroundColor',[1 1 1],...
                'EdgeColor',[0 0 0],...
                'Interpreter','latex',...
                'String',strcat('age $\in$[',num2str(partition_idx-1),'/',num2str(deglifePartitions_num),',',num2str(partition_idx),'/',num2str(deglifePartitions_num),']'),...
                'Position',textpositions{partition_idx},...
                'Color',[0 0 0],'FontSize',fontSize);
        end
        colormap jet;
        caxis([data_min data_max]);
    end
    labelTicks = round(linspace(data_min,data_max,5),2);
    h = colorbar('Position',[0.888 0.068 0.0266666666666669 0.892],...
        'TickLabelInterpreter','latex','FontSize',fontSize,'Ticks',labelTicks,'Limits',[0,round(data_max,2)]);
    caxis([data_min data_max]);
    ylabel(h, 'Rel. resistive en. loss ($\%$)','Interpreter','latex','FontSize',fontSize);
    colormap jet;
    fig_pos = figure_1.PaperPosition;
    figure_1.PaperSize = [fig_pos(3) fig_pos(4)];
    if(savefiles)
        saveas(gcf,strcat(filename,'.fig'));
        saveas(gcf,strcat(filename,'.eps'),'epsc');
    end
end

plot_life_partitions_capacity_loss_percentage = 1;
if(show_plots&&plot_life_partitions_capacity_loss_percentage)
    filename = 'capLoss';
    % plot figure--
    partition_idxs = 1:deglifePartitions_num;
    positions = {[0.0924678097197633 0.246165556113454 0.20936667244291 0.615405872457975],...
        [0.373264911169039 0.238165556113454 0.20936667244291 0.623405872457975],...
        [0.654062012618313 0.238165556113454 0.209366672442911 0.627405872457975]};
    
    textpositions = {[0.232218011771345 1.10249990460461 0],...
        [0.181262597758605 1.10233123498622 0],...
        [0.24495686527453 1.09583826559003 0]...
        };
        
    y_axis_grid = cell_pow_set/(cell_1C_capacityInAh*cell_nominalVoltage);
    x_axis_grid = unique([0,soc_grid_boundaries,1]);
    x_axis_grid_mean = zeros(length(x_axis_grid)-1,1);
    for idx = 1:length(x_axis_grid)-1
        x_axis_grid_mean(idx) = 0.5*(x_axis_grid(idx) +x_axis_grid(idx+1));
    end
    
    [x,y] = meshgrid(x_axis_grid,y_axis_grid);
    
    plot_data = nan(length(y_axis_grid),length(x_axis_grid),deglifePartitions_num);
    for partition_idx = 1:deglifePartitions_num
        for soc_bin_idx = 1:z_num
%             plot_data(:,soc_bin_idx+1,partition_idx) = mean_capacityLoss_factor_comsol(soc_bin_idx,:,partition_idx)*100*100;
            plot_data(:,soc_bin_idx+1,partition_idx) = mean_capacityLoss_factor_incl_calender_ageing(soc_bin_idx,:,partition_idx)*100*100;
        end
    end
%     plot_data(:,z_num+2,:) = mean_capacityLoss_factor_comsol(z_num,:,:)*100*100;
    plot_data(:,z_num+2,:) = mean_capacityLoss_factor_incl_calender_ageing(z_num,:,:)*100*100;
        
    data_min = min(plot_data(:));
    data_max = round(max(plot_data(:))*10^8)*10^-8;
    
    fontSize = 12;
    figure_xsize = 700;
    figure_ysize = 200;
    figure_position_offset = 50;
    figure_1 = figure('Color','w','Renderer','painters');
    set(gcf, 'Position',  [figure_position_offset,figure_position_offset,figure_position_offset+figure_xsize, figure_position_offset+figure_ysize]);
    
    for partition_idx = partition_idxs
        color_data_temp = cell(length(x_axis_grid)-1,1);
        for i = 1:length(x_axis_grid)-1
            color_data_temp{i} = [plot_data(:,i,partition_idx),plot_data(:,i,partition_idx)];
        end
        axes1 = axes('Parent',figure_1,...
            'Position',positions{partition_idx});
        box(axes1,'on');
        hold(axes1,'on');
        h1 = surf(x,y,0*plot_data(:,:,partition_idx));
        xline(cell_SOC_low,'Alpha',1,'Color',[1 0 0],'LineStyle','--',...
            'LineWidth',1);
        xline(cell_SOC_high,'Alpha',1,'Color',[1 0 0],'LineStyle','--',...
            'LineWidth',1);
        ylabel({'Current (A)'},'Interpreter','latex','Units', 'normalized','position',[-0.113092878602192,0.525806928450061,0]);
        xlabel({'SOC'},'Interpreter','latex','Units','normalized');
        set(h1,'CData',plot_data(:,:,partition_idx),'FaceColor','interp','FaceAlpha',1)
        shading interp;
        axes1.YAxis.MinorTickValues = -1:0.5:1;
        xlim([0,1]);
        ylim([-1 ,1]);
        set(axes1,'FontSize',fontSize,'TickLabelInterpreter','latex',...
            'XMinorTick','on','YMinorTick','on','XTick',0:0.25:1,'YTick',[-1  0 1],...
            'Colormap',...
            [0.2422 0.1504 0.6603;0.249283333333333 0.162769027777778 0.701203472222222;0.255735 0.177103055555556 0.738961944444444;0.261851111111111 0.1910475 0.776739166666667;0.267449444444444 0.205018333333333 0.813729444444444;0.272160416666667 0.220381944444444 0.847371527777778;0.2756675 0.237689166666667 0.87573;0.278331666666667 0.256401111111111 0.899388055555556;0.280133333333333 0.275595555555556 0.919547777777778;0.281129583333333 0.29485625 0.936625833333333;0.281240277777778 0.313958333333333 0.951533333333333;0.280573888888889 0.33292625 0.964617083333333;0.27847 0.351745555555556 0.975510555555556;0.274978472222222 0.370759027777778 0.984490972222222;0.269529444444444 0.390201111111111 0.990806111111111;0.261066666666667 0.410085416666667 0.994735416666667;0.248322222222222 0.430242222222222 0.998015555555556;0.229866527777778 0.450917916666667 0.998971527777778;0.208676666666666 0.472239166666667 0.993345;0.190192638888889 0.492728888888889 0.986256527777778;0.182322222222222 0.511811111111111 0.977522222222222;0.178614166666667 0.530272361111111 0.966969027777778;0.176563611111111 0.548336944444445 0.953331666666666;0.170852916666667 0.565923611111111 0.939202083333333;0.159056666666667 0.583216666666667 0.926856666666666;0.149183333333333 0.599913194444445 0.9145625;0.143464166666667 0.615898055555556 0.903748333333333;0.13510375 0.63177875 0.8953375;0.123628888888889 0.647465555555556 0.887368888888889;0.111881111111111 0.662524444444445 0.876855833333333;0.0986847222222221 0.676755555555556 0.863209722222222;0.0786555555555555 0.689891944444444 0.84639375;0.0486777777777777 0.701922222222222 0.827386666666667;0.0163675000000002 0.713000833333333 0.806969166666667;0.00277861111111108 0.723041388888889 0.785363888888889;0.0102388888888889 0.732311111111111 0.763129861111111;0.0416083333333331 0.740801111111111 0.740251111111111;0.0873354166666661 0.748534027777778 0.716729027777778;0.127959166666666 0.755790277777778 0.692834166666667;0.159025416666666 0.763085 0.668474166666667;0.181144444444444 0.770561111111111 0.64318888888889;0.199190555555555 0.778161805555555 0.616446944444446;0.219814999999999 0.7853025 0.587710000000001;0.247202638888887 0.791606666666666 0.55696402777778;0.283491111111109 0.796596666666666 0.524602222222224;0.326124999999997 0.799897222222222 0.490516666666669;0.369233611111108 0.802148611111111 0.454283333333336;0.413774444444441 0.803017777777778 0.416086250000003;0.462759999999995 0.8016 0.37836666666667;0.513110277777773 0.798542916666667 0.341453194444448;0.562234722222217 0.794126388888889 0.303961111111115;0.610257083333328 0.788480000000001 0.266809166666671;0.675033333333329 0.778621666666667 0.221158333333337;0.689467934782608 0.77620875 0.211704347826087;0.696378913043478 0.7749825 0.207469130434783;0.703221739130435 0.77375625 0.203336141304348;0.710035652173913 0.772519565217391 0.199262173913043;0.71680972826087 0.771243043478261 0.195307717391304;0.723566739130435 0.76996347826087 0.191390434782609;0.730284456521739 0.768683913043478 0.187551739130435;0.736956521739131 0.76738152173913 0.183758695652174;0.743570000000001 0.766050815217391 0.180027173913043;0.75012847826087 0.764717934782609 0.176426630434782;0.756632934782609 0.763385054347826 0.172961141304347;0.763069565217392 0.762052173913043 0.169699130434782;0.769464945652175 0.760719293478261 0.166575543478261;0.775824673913044 0.759386413043478 0.163772826086956;0.78215043478261 0.758053532608695 0.161182065217391;0.788372173913044 0.756720652173913 0.159007391304347;0.794556250000001 0.75538777173913 0.156988260869565;0.800687500000001 0.754054891304348 0.155708695652174;0.806780271739132 0.752734836956521 0.154595869565217;0.812768695652175 0.751449565217391 0.153935217391304;0.818691739130436 0.75017 0.153541793478261;0.824510108695654 0.748890434782608 0.153567065217391;0.830301086956523 0.747610869565217 0.153803804347826;0.836059130434784 0.746331304347826 0.154336956521739;0.841761630434784 0.745107282608695 0.155092282608696;0.847412173913045 0.743932826086956 0.156057826086957;0.853032336956523 0.742759891304347 0.15717375;0.85863043478261 0.741586956521739 0.1584;0.864147989130437 0.74045429347826 0.159827608695653;0.869633152173915 0.739337826086956 0.161342500000001;0.875037119565219 0.738261956521739 0.163141576086957;0.88040239130435 0.737205434782608 0.16510543478261;0.885699728260872 0.736182880434782 0.167408967391305;0.890964239130437 0.735183586956521 0.169869782608697;0.896135815217393 0.734277228260869 0.172748804347827;0.901244782608698 0.733402173913043 0.175815652173915;0.906211793478263 0.732598097826087 0.179308369565219;0.911126630434785 0.731841847826087 0.183022826086958;0.915981576086959 0.731145489130434 0.187006793478263;0.920753043478263 0.730532608695652 0.191219347826089;0.925444782608698 0.729999456521739 0.195644510869567;0.930036413043481 0.729533043478261 0.200236521739133;0.934560597826089 0.729117934782608 0.20493586956522;0.939020000000003 0.7288 0.209700000000003;0.943445163043481 0.728550380434782 0.214481413043481;0.947870326086959 0.728453369565217 0.219110217391307;0.952295489130438 0.728407771739131 0.223652880434785;0.956720652173916 0.728453260869565 0.227694565217394;0.961145815217394 0.728586467391305 0.231506793478263;0.965570978260873 0.728919239130435 0.23485336956522;0.969969728260873 0.72935847826087 0.237894293478263;0.974288260869569 0.729998260869566 0.240293478260871;0.978428260869568 0.730915760869566 0.24235543478261;0.982266956521742 0.73230195652174 0.243848260869566;0.985722173913046 0.733934673913045 0.244108478260869;0.988803913043481 0.735806086956523 0.243208478260869;0.991522500000002 0.737890108695654 0.241741576086955;0.993975000000002 0.740129347826089 0.239875543478259;0.995394565217392 0.74272614130435 0.237691684782607;0.996382608695653 0.745466086956524 0.235404347826085;0.996831793478261 0.748349728260872 0.233152934782606;0.997064456521739 0.751290543478264 0.230913695652172;0.997152989130435 0.754267391304351 0.228674456521737;0.997193695652174 0.75725934782609 0.226428913043476;0.997140380434783 0.760298315217395 0.22413635869565;0.997070217391304 0.763337282608699 0.221843804347824;0.9969675 0.766376250000003 0.219551249999997;0.996841304347826 0.769415217391308 0.217258695652171;0.996681358695652 0.772454184782612 0.214966141304345;0.996495217391304 0.775519347826091 0.212673586956519;0.996281956521739 0.778611630434786 0.210381032608693;0.996035652173913 0.781703913043482 0.208121521739128;0.995757065217391 0.784800271739134 0.205886141304345;0.995384347826086 0.787927934782613 0.203682282608693;0.994954402173912 0.7910735326087 0.201499782608693;0.99439304347826 0.794219130434787 0.199448695652171;0.993780054347825 0.797371358695656 0.19744260869565;0.993009239130434 0.800563043478265 0.195554891304345;0.992184076086955 0.803774673913048 0.193714293478258;0.99129108695652 0.807020217391309 0.191941521739128;0.990334076086955 0.810272445652179 0.190207445652171;0.989214456521738 0.813524673913048 0.188554673913041;0.988033695652172 0.816797282608701 0.186942663043476;0.986757391304346 0.82010173913044 0.185394347826085;0.985425978260867 0.823407282608701 0.183848206521737;0.984038913043476 0.826712826086962 0.18230163043478;0.982588478260867 0.830018369565223 0.180723369565215;0.981095652173911 0.833323913043484 0.179123913043476;0.979563641304345 0.836629456521745 0.177524456521736;0.978011304347823 0.83994119565218 0.175918804347823;0.976425543478258 0.843286358695658 0.174279728260867;0.974838260869562 0.846639130434788 0.172626956521736;0.973331521739128 0.849951630434788 0.170974184782606;0.971854782608693 0.853257173913049 0.169329456521736;0.970464021739128 0.85656271739131 0.167727717391301;0.969128695652172 0.85988521739131 0.166145217391301;0.967888913043476 0.863239456521745 0.164594565217388;0.966719565217389 0.866579347826093 0.163067391304345;0.965653260869563 0.869884891304354 0.161574565217388;0.964663913043476 0.873190434782615 0.160133043478258;0.963757554347824 0.876495978260876 0.158746847826084;0.962915760869564 0.879801521739137 0.157392934782606;0.962129076086955 0.883107065217398 0.156060054347823;0.961471956521738 0.886412608695659 0.154727173913041;0.960885489130434 0.88971815217392 0.153394293478258;0.960481195652173 0.89297815217392 0.152061413043475;0.960034782608695 0.897082608695658 0.150347826086954;0.959957336956521 0.898011956521745 0.149929619565215;0.959879891304347 0.898941304347831 0.149511413043476;0.959802445652173 0.899870652173918 0.149093206521737;0.959724999999999 0.900800000000005 0.148674999999998;0.959647554347826 0.901729347826092 0.148256793478258;0.959615760869565 0.902670108695658 0.147804347826084;0.959586503623188 0.903611503623194 0.147349999999997;0.959563043478261 0.904554347826093 0.146891304347823;0.959547554347826 0.905499184782615 0.14642663043478;0.959532065217391 0.906444021739136 0.145961956521736;0.959525815217391 0.907382699275368 0.145491123188403;0.959520652173913 0.908320652173919 0.145019565217388;0.959528804347826 0.909249728260875 0.144539130434779;0.959559782608696 0.910163586956528 0.144043478260866;0.959590760869566 0.91107744565218 0.143547826086953;0.959621739130435 0.911991304347832 0.14305217391304;0.959652717391305 0.912905163043484 0.142556521739127;0.959690217391305 0.913822282608702 0.142047826086953;0.959741847826087 0.91474646739131 0.141510869565214;0.95979347826087 0.915670652173919 0.140973913043475;0.959853260869566 0.916598913043485 0.140420652173909;0.959915217391305 0.917528260869571 0.139863043478257;0.959978532608696 0.918457155797108 0.139304528985503;0.960045652173914 0.919384782608702 0.138743478260866;0.960112771739131 0.920312409420296 0.138182427536228;0.960210326086957 0.921229891304354 0.137601086956518;0.960318750000001 0.922143750000006 0.137012499999996;0.960427173913044 0.923057608695659 0.136423913043474;0.960535597826088 0.923971467391311 0.135835326086952;0.960644021739131 0.924885326086963 0.135246739130431;0.960770108695653 0.925793297101456 0.13464637681159;0.960904347826088 0.926698550724644 0.134040579710141;0.961041847826088 0.927602717391311 0.133432608695648;0.961196739130436 0.928501086956528 0.132813043478256;0.961351630434784 0.929399456521746 0.132193478260865;0.961506521739132 0.930297826086963 0.131573913043474;0.961661413043479 0.931196195652181 0.130954347826082;0.96181793478261 0.932096195652181 0.130329891304343;0.961988315217393 0.933010054347833 0.129663858695647;0.962158695652175 0.933923913043485 0.128997826086952;0.962329076086958 0.934837771739137 0.128331793478256;0.96249945652174 0.935751630434789 0.12766576086956;0.962670561594204 0.936665126811601 0.126999003623183;0.962854710144929 0.937572101449282 0.126319202898546;0.963038858695654 0.938479076086963 0.125639402173908;0.963232065217393 0.939381521739137 0.124950543478256;0.963433423913045 0.940279891304355 0.12425353260869;0.963634782608697 0.941178260869572 0.123556521739125;0.963836141304349 0.94207663043479 0.12285951086956;0.964037500000002 0.942975000000007 0.122162499999995;0.964251902173915 0.943866847826094 0.121445923913038;0.964480797101451 0.944751449275369 0.120707608695646;0.964709692028987 0.945636050724645 0.119969293478255;0.964941847826089 0.946519021739137 0.119226086956516;0.965174184782611 0.94740190217392 0.118482608695646;0.965410869565219 0.948286956521746 0.117730434782603;0.965653532608698 0.949175000000007 0.11696630434782;0.965896195652176 0.950063043478268 0.116202173913037;0.966157336956524 0.950960326086964 0.115401086956515;0.966420652173915 0.951858695652181 0.114595652173906;0.966683967391307 0.952757065217399 0.113790217391298;0.966947282608698 0.953655434782616 0.112984782608689;0.967210597826089 0.954553804347834 0.11217934782608;0.967484057971017 0.955442028985515 0.111343478260863;0.967759420289857 0.956328351449283 0.110501902173906;0.96803586956522 0.957213586956529 0.109657065217384;0.968314673913046 0.958096467391312 0.108805163043471;0.968593478260872 0.958979347826095 0.107953260869558;0.968875000000002 0.959864945652182 0.107087771739123;0.969157246376814 0.96075126811595 0.106218659420282;0.969442663043481 0.961640760869573 0.105333695652166;0.969736956521742 0.96253913043479 0.104404347826079;0.970031250000003 0.963437500000008 0.103474999999992;0.970325543478263 0.964335869565225 0.102545652173905;0.970619836956524 0.965234239130443 0.101616304347818;0.970914130434785 0.966130434782617 0.1006804347826;0.971208423913046 0.967018478260878 0.0997201086956437;0.971502717391307 0.967906521739138 0.098759782608687;0.971797010869568 0.968791032608704 0.0977888586956434;0.972091304347829 0.969673913043486 0.0968130434782521;0.97238559782609 0.970557065217399 0.0958364130434694;0.972679891304351 0.971441666666675 0.0948554347825997;0.972974184782611 0.97232626811595 0.0938744565217299;0.973268478260872 0.9732195652174 0.0928673913043385;0.973562771739133 0.974117934782617 0.0918451086956428;0.973857065217394 0.975016304347834 0.0908228260869471;0.974151358695655 0.975914673913052 0.0898005434782514;0.974445652173916 0.976813043478269 0.0887782608695557;0.974744927536235 0.977706431159429 0.087750996376802;0.975047826086959 0.978596195652182 0.0867201086956425;0.975351086956525 0.979485597826095 0.0856888586956424;0.97566086956522 0.980368478260878 0.0846510869565119;0.975970652173916 0.98125135869566 0.0836133152173814;0.976280434782611 0.982134239130442 0.0825755434782525;0.976590217391306 0.983017119565222 0.0815377717391251;0.9769 0.9839 0.0805],...
    'YTickLabel',{'-1C','0','1C'});
        if(partition_idx == 1)
            text('Units','normalized','BackgroundColor',[1 1 1],...
                'EdgeColor',[0 0 0],...
                'Interpreter','latex',...
                'String',strcat('age $\in$[0,1/',num2str(deglifePartitions_num),']'),...
                'Position',textpositions{partition_idx},...
                'Color',[0 0 0],'FontSize',fontSize);
        elseif(partition_idx == deglifePartitions_num)
            text('Units','normalized','BackgroundColor',[1 1 1],...
                'EdgeColor',[0 0 0],...
                'Interpreter','latex',...
                'String',strcat('age $\in$[',num2str(deglifePartitions_num-1),'/',num2str(deglifePartitions_num),',1]'),...
                'Position',textpositions{partition_idx},...
                'Color',[0 0 0],'FontSize',fontSize);
        else
            text('Units','normalized','BackgroundColor',[1 1 1],...
                'EdgeColor',[0 0 0],...
                'Interpreter','latex',...
                'String',strcat('age $\in$[',num2str(partition_idx-1),'/',num2str(deglifePartitions_num),',',num2str(partition_idx),'/',num2str(deglifePartitions_num),']'),...
                'Position',textpositions{partition_idx},...
                'Color',[0 0 0],'FontSize',fontSize);
        end
%                 colormap jet;
        caxis([data_min data_max]);
    end
    precisionDigits_plot = 2;
    labelTicks = round(linspace(round(data_min,precisionDigits_plot),round(data_max,precisionDigits_plot),5),precisionDigits_plot); % 1 hour
    h = colorbar('peer',axes1,'Position',[0.888 0.064 0.0266666666666666 0.816],...
        'TickLabelInterpreter','latex',...
        'Colormap',...
        [0.2422 0.1504 0.6603;0.249283333333333 0.162769027777778 0.701203472222222;0.255735 0.177103055555556 0.738961944444444;0.261851111111111 0.1910475 0.776739166666667;0.267449444444444 0.205018333333333 0.813729444444444;0.272160416666667 0.220381944444444 0.847371527777778;0.2756675 0.237689166666667 0.87573;0.278331666666667 0.256401111111111 0.899388055555556;0.280133333333333 0.275595555555556 0.919547777777778;0.281129583333333 0.29485625 0.936625833333333;0.281240277777778 0.313958333333333 0.951533333333333;0.280573888888889 0.33292625 0.964617083333333;0.27847 0.351745555555556 0.975510555555556;0.274978472222222 0.370759027777778 0.984490972222222;0.269529444444444 0.390201111111111 0.990806111111111;0.261066666666667 0.410085416666667 0.994735416666667;0.248322222222222 0.430242222222222 0.998015555555556;0.229866527777778 0.450917916666667 0.998971527777778;0.208676666666666 0.472239166666667 0.993345;0.190192638888889 0.492728888888889 0.986256527777778;0.182322222222222 0.511811111111111 0.977522222222222;0.178614166666667 0.530272361111111 0.966969027777778;0.176563611111111 0.548336944444445 0.953331666666666;0.170852916666667 0.565923611111111 0.939202083333333;0.159056666666667 0.583216666666667 0.926856666666666;0.149183333333333 0.599913194444445 0.9145625;0.143464166666667 0.615898055555556 0.903748333333333;0.13510375 0.63177875 0.8953375;0.123628888888889 0.647465555555556 0.887368888888889;0.111881111111111 0.662524444444445 0.876855833333333;0.0986847222222221 0.676755555555556 0.863209722222222;0.0786555555555555 0.689891944444444 0.84639375;0.0486777777777777 0.701922222222222 0.827386666666667;0.0163675000000002 0.713000833333333 0.806969166666667;0.00277861111111108 0.723041388888889 0.785363888888889;0.0102388888888889 0.732311111111111 0.763129861111111;0.0416083333333331 0.740801111111111 0.740251111111111;0.0873354166666661 0.748534027777778 0.716729027777778;0.127959166666666 0.755790277777778 0.692834166666667;0.159025416666666 0.763085 0.668474166666667;0.181144444444444 0.770561111111111 0.64318888888889;0.199190555555555 0.778161805555555 0.616446944444446;0.219814999999999 0.7853025 0.587710000000001;0.247202638888887 0.791606666666666 0.55696402777778;0.283491111111109 0.796596666666666 0.524602222222224;0.326124999999997 0.799897222222222 0.490516666666669;0.369233611111108 0.802148611111111 0.454283333333336;0.413774444444441 0.803017777777778 0.416086250000003;0.462759999999995 0.8016 0.37836666666667;0.513110277777773 0.798542916666667 0.341453194444448;0.562234722222217 0.794126388888889 0.303961111111115;0.610257083333328 0.788480000000001 0.266809166666671;0.675033333333329 0.778621666666667 0.221158333333337;0.689467934782608 0.77620875 0.211704347826087;0.696378913043478 0.7749825 0.207469130434783;0.703221739130435 0.77375625 0.203336141304348;0.710035652173913 0.772519565217391 0.199262173913043;0.71680972826087 0.771243043478261 0.195307717391304;0.723566739130435 0.76996347826087 0.191390434782609;0.730284456521739 0.768683913043478 0.187551739130435;0.736956521739131 0.76738152173913 0.183758695652174;0.743570000000001 0.766050815217391 0.180027173913043;0.75012847826087 0.764717934782609 0.176426630434782;0.756632934782609 0.763385054347826 0.172961141304347;0.763069565217392 0.762052173913043 0.169699130434782;0.769464945652175 0.760719293478261 0.166575543478261;0.775824673913044 0.759386413043478 0.163772826086956;0.78215043478261 0.758053532608695 0.161182065217391;0.788372173913044 0.756720652173913 0.159007391304347;0.794556250000001 0.75538777173913 0.156988260869565;0.800687500000001 0.754054891304348 0.155708695652174;0.806780271739132 0.752734836956521 0.154595869565217;0.812768695652175 0.751449565217391 0.153935217391304;0.818691739130436 0.75017 0.153541793478261;0.824510108695654 0.748890434782608 0.153567065217391;0.830301086956523 0.747610869565217 0.153803804347826;0.836059130434784 0.746331304347826 0.154336956521739;0.841761630434784 0.745107282608695 0.155092282608696;0.847412173913045 0.743932826086956 0.156057826086957;0.853032336956523 0.742759891304347 0.15717375;0.85863043478261 0.741586956521739 0.1584;0.864147989130437 0.74045429347826 0.159827608695653;0.869633152173915 0.739337826086956 0.161342500000001;0.875037119565219 0.738261956521739 0.163141576086957;0.88040239130435 0.737205434782608 0.16510543478261;0.885699728260872 0.736182880434782 0.167408967391305;0.890964239130437 0.735183586956521 0.169869782608697;0.896135815217393 0.734277228260869 0.172748804347827;0.901244782608698 0.733402173913043 0.175815652173915;0.906211793478263 0.732598097826087 0.179308369565219;0.911126630434785 0.731841847826087 0.183022826086958;0.915981576086959 0.731145489130434 0.187006793478263;0.920753043478263 0.730532608695652 0.191219347826089;0.925444782608698 0.729999456521739 0.195644510869567;0.930036413043481 0.729533043478261 0.200236521739133;0.934560597826089 0.729117934782608 0.20493586956522;0.939020000000003 0.7288 0.209700000000003;0.943445163043481 0.728550380434782 0.214481413043481;0.947870326086959 0.728453369565217 0.219110217391307;0.952295489130438 0.728407771739131 0.223652880434785;0.956720652173916 0.728453260869565 0.227694565217394;0.961145815217394 0.728586467391305 0.231506793478263;0.965570978260873 0.728919239130435 0.23485336956522;0.969969728260873 0.72935847826087 0.237894293478263;0.974288260869569 0.729998260869566 0.240293478260871;0.978428260869568 0.730915760869566 0.24235543478261;0.982266956521742 0.73230195652174 0.243848260869566;0.985722173913046 0.733934673913045 0.244108478260869;0.988803913043481 0.735806086956523 0.243208478260869;0.991522500000002 0.737890108695654 0.241741576086955;0.993975000000002 0.740129347826089 0.239875543478259;0.995394565217392 0.74272614130435 0.237691684782607;0.996382608695653 0.745466086956524 0.235404347826085;0.996831793478261 0.748349728260872 0.233152934782606;0.997064456521739 0.751290543478264 0.230913695652172;0.997152989130435 0.754267391304351 0.228674456521737;0.997193695652174 0.75725934782609 0.226428913043476;0.997140380434783 0.760298315217395 0.22413635869565;0.997070217391304 0.763337282608699 0.221843804347824;0.9969675 0.766376250000003 0.219551249999997;0.996841304347826 0.769415217391308 0.217258695652171;0.996681358695652 0.772454184782612 0.214966141304345;0.996495217391304 0.775519347826091 0.212673586956519;0.996281956521739 0.778611630434786 0.210381032608693;0.996035652173913 0.781703913043482 0.208121521739128;0.995757065217391 0.784800271739134 0.205886141304345;0.995384347826086 0.787927934782613 0.203682282608693;0.994954402173912 0.7910735326087 0.201499782608693;0.99439304347826 0.794219130434787 0.199448695652171;0.993780054347825 0.797371358695656 0.19744260869565;0.993009239130434 0.800563043478265 0.195554891304345;0.992184076086955 0.803774673913048 0.193714293478258;0.99129108695652 0.807020217391309 0.191941521739128;0.990334076086955 0.810272445652179 0.190207445652171;0.989214456521738 0.813524673913048 0.188554673913041;0.988033695652172 0.816797282608701 0.186942663043476;0.986757391304346 0.82010173913044 0.185394347826085;0.985425978260867 0.823407282608701 0.183848206521737;0.984038913043476 0.826712826086962 0.18230163043478;0.982588478260867 0.830018369565223 0.180723369565215;0.981095652173911 0.833323913043484 0.179123913043476;0.979563641304345 0.836629456521745 0.177524456521736;0.978011304347823 0.83994119565218 0.175918804347823;0.976425543478258 0.843286358695658 0.174279728260867;0.974838260869562 0.846639130434788 0.172626956521736;0.973331521739128 0.849951630434788 0.170974184782606;0.971854782608693 0.853257173913049 0.169329456521736;0.970464021739128 0.85656271739131 0.167727717391301;0.969128695652172 0.85988521739131 0.166145217391301;0.967888913043476 0.863239456521745 0.164594565217388;0.966719565217389 0.866579347826093 0.163067391304345;0.965653260869563 0.869884891304354 0.161574565217388;0.964663913043476 0.873190434782615 0.160133043478258;0.963757554347824 0.876495978260876 0.158746847826084;0.962915760869564 0.879801521739137 0.157392934782606;0.962129076086955 0.883107065217398 0.156060054347823;0.961471956521738 0.886412608695659 0.154727173913041;0.960885489130434 0.88971815217392 0.153394293478258;0.960481195652173 0.89297815217392 0.152061413043475;0.960034782608695 0.897082608695658 0.150347826086954;0.959957336956521 0.898011956521745 0.149929619565215;0.959879891304347 0.898941304347831 0.149511413043476;0.959802445652173 0.899870652173918 0.149093206521737;0.959724999999999 0.900800000000005 0.148674999999998;0.959647554347826 0.901729347826092 0.148256793478258;0.959615760869565 0.902670108695658 0.147804347826084;0.959586503623188 0.903611503623194 0.147349999999997;0.959563043478261 0.904554347826093 0.146891304347823;0.959547554347826 0.905499184782615 0.14642663043478;0.959532065217391 0.906444021739136 0.145961956521736;0.959525815217391 0.907382699275368 0.145491123188403;0.959520652173913 0.908320652173919 0.145019565217388;0.959528804347826 0.909249728260875 0.144539130434779;0.959559782608696 0.910163586956528 0.144043478260866;0.959590760869566 0.91107744565218 0.143547826086953;0.959621739130435 0.911991304347832 0.14305217391304;0.959652717391305 0.912905163043484 0.142556521739127;0.959690217391305 0.913822282608702 0.142047826086953;0.959741847826087 0.91474646739131 0.141510869565214;0.95979347826087 0.915670652173919 0.140973913043475;0.959853260869566 0.916598913043485 0.140420652173909;0.959915217391305 0.917528260869571 0.139863043478257;0.959978532608696 0.918457155797108 0.139304528985503;0.960045652173914 0.919384782608702 0.138743478260866;0.960112771739131 0.920312409420296 0.138182427536228;0.960210326086957 0.921229891304354 0.137601086956518;0.960318750000001 0.922143750000006 0.137012499999996;0.960427173913044 0.923057608695659 0.136423913043474;0.960535597826088 0.923971467391311 0.135835326086952;0.960644021739131 0.924885326086963 0.135246739130431;0.960770108695653 0.925793297101456 0.13464637681159;0.960904347826088 0.926698550724644 0.134040579710141;0.961041847826088 0.927602717391311 0.133432608695648;0.961196739130436 0.928501086956528 0.132813043478256;0.961351630434784 0.929399456521746 0.132193478260865;0.961506521739132 0.930297826086963 0.131573913043474;0.961661413043479 0.931196195652181 0.130954347826082;0.96181793478261 0.932096195652181 0.130329891304343;0.961988315217393 0.933010054347833 0.129663858695647;0.962158695652175 0.933923913043485 0.128997826086952;0.962329076086958 0.934837771739137 0.128331793478256;0.96249945652174 0.935751630434789 0.12766576086956;0.962670561594204 0.936665126811601 0.126999003623183;0.962854710144929 0.937572101449282 0.126319202898546;0.963038858695654 0.938479076086963 0.125639402173908;0.963232065217393 0.939381521739137 0.124950543478256;0.963433423913045 0.940279891304355 0.12425353260869;0.963634782608697 0.941178260869572 0.123556521739125;0.963836141304349 0.94207663043479 0.12285951086956;0.964037500000002 0.942975000000007 0.122162499999995;0.964251902173915 0.943866847826094 0.121445923913038;0.964480797101451 0.944751449275369 0.120707608695646;0.964709692028987 0.945636050724645 0.119969293478255;0.964941847826089 0.946519021739137 0.119226086956516;0.965174184782611 0.94740190217392 0.118482608695646;0.965410869565219 0.948286956521746 0.117730434782603;0.965653532608698 0.949175000000007 0.11696630434782;0.965896195652176 0.950063043478268 0.116202173913037;0.966157336956524 0.950960326086964 0.115401086956515;0.966420652173915 0.951858695652181 0.114595652173906;0.966683967391307 0.952757065217399 0.113790217391298;0.966947282608698 0.953655434782616 0.112984782608689;0.967210597826089 0.954553804347834 0.11217934782608;0.967484057971017 0.955442028985515 0.111343478260863;0.967759420289857 0.956328351449283 0.110501902173906;0.96803586956522 0.957213586956529 0.109657065217384;0.968314673913046 0.958096467391312 0.108805163043471;0.968593478260872 0.958979347826095 0.107953260869558;0.968875000000002 0.959864945652182 0.107087771739123;0.969157246376814 0.96075126811595 0.106218659420282;0.969442663043481 0.961640760869573 0.105333695652166;0.969736956521742 0.96253913043479 0.104404347826079;0.970031250000003 0.963437500000008 0.103474999999992;0.970325543478263 0.964335869565225 0.102545652173905;0.970619836956524 0.965234239130443 0.101616304347818;0.970914130434785 0.966130434782617 0.1006804347826;0.971208423913046 0.967018478260878 0.0997201086956437;0.971502717391307 0.967906521739138 0.098759782608687;0.971797010869568 0.968791032608704 0.0977888586956434;0.972091304347829 0.969673913043486 0.0968130434782521;0.97238559782609 0.970557065217399 0.0958364130434694;0.972679891304351 0.971441666666675 0.0948554347825997;0.972974184782611 0.97232626811595 0.0938744565217299;0.973268478260872 0.9732195652174 0.0928673913043385;0.973562771739133 0.974117934782617 0.0918451086956428;0.973857065217394 0.975016304347834 0.0908228260869471;0.974151358695655 0.975914673913052 0.0898005434782514;0.974445652173916 0.976813043478269 0.0887782608695557;0.974744927536235 0.977706431159429 0.087750996376802;0.975047826086959 0.978596195652182 0.0867201086956425;0.975351086956525 0.979485597826095 0.0856888586956424;0.97566086956522 0.980368478260878 0.0846510869565119;0.975970652173916 0.98125135869566 0.0836133152173814;0.976280434782611 0.982134239130442 0.0825755434782525;0.976590217391306 0.983017119565222 0.0815377717391251;0.9769 0.9839 0.0805],...
    'FontSize',fontSize,'Ticks',labelTicks,'Limits',[0,round(data_max,precisionDigits_plot)]); % 1 hour
    caxis([data_min data_max]);
%         colormap jet;
    ylabel(h, 'Rel. capacity loss ($\times 10^{-2}\,\%$)','Interpreter','latex','FontSize',fontSize);
    fig_pos = figure_1.PaperPosition;
    figure_1.PaperSize = [fig_pos(3) fig_pos(4)];
    if(savefiles)
        saveas(gcf,strcat(filename,'.fig'));
        saveas(gcf,strcat(filename,'.eps'),'epsc');
    end
end

plot_life_partitions_soc_kp1_mean = 1;
if(show_plots&&plot_life_partitions_soc_kp1_mean)
    filename = 'socMean';
    partition_idxs = 1:deglifePartitions_num;
    positions = {[0.0924678097197633 0.246165556113454 0.20936667244291 0.615405872457975],...
        [0.373264911169039 0.238165556113454 0.20936667244291 0.623405872457975],...
        [0.654062012618313 0.238165556113454 0.209366672442911 0.627405872457975]};
    
    textpositions = {[0.232218011771345 1.10249990460461 0],...
        [0.181262597758605 1.10233123498622 0],...
        [0.24495686527453 1.09583826559003 0]...
        };
    
    data_min = min(mean_soc_kp1_3c(:));
    data_max = round(max(mean_soc_kp1_3c(:))*10^6)*10^-6;
    
    y_axis_grid = cell_pow_set/(cell_1C_capacityInAh*cell_nominalVoltage);
    x_axis_grid = unique([0,soc_grid_boundaries,1]);
    x_axis_grid_mean = zeros(length(x_axis_grid)-1,1);
    for idx = 1:length(x_axis_grid)-1
        x_axis_grid_mean(idx) = 0.5*(x_axis_grid(idx) +x_axis_grid(idx+1));
    end
    
    [x,y] = meshgrid(x_axis_grid,y_axis_grid);
    
    plot_data = nan(length(y_axis_grid),length(x_axis_grid),deglifePartitions_num);
    for partition_idx = 1:deglifePartitions_num
        for soc_bin_idx = 1:z_num
            plot_data(:,soc_bin_idx+1,partition_idx) = mean_soc_kp1_3c(soc_bin_idx,:,partition_idx);
        end
    end
    plot_data(:,z_num+2,:) = mean_soc_kp1_3c(z_num,:,:);
    
    fontSize = 12;
    figure_xsize = 700;
    figure_ysize = 200;
    figure_position_offset = 50;
    figure_1 = figure('Color','w','Renderer','painters');
    set(gcf, 'Position',  [figure_position_offset,figure_position_offset,figure_position_offset+figure_xsize, figure_position_offset+figure_ysize]);
    
    for partition_idx = partition_idxs
        color_data_temp = cell(length(x_axis_grid)-1,1);
        for i = 1:length(x_axis_grid)-1
            color_data_temp{i} = [plot_data(:,i,partition_idx),plot_data(:,i,partition_idx)];
        end
        axes1 = axes('Parent',figure_1,...
            'Position',positions{partition_idx});
        box(axes1,'on');
        hold(axes1,'on');
        h1 = surf(x,y,0*plot_data(:,:,partition_idx));
        xline(cell_SOC_low,'Alpha',1,'Color',[1 0 0],'LineStyle','--',...
            'LineWidth',1);
        xline(cell_SOC_high,'Alpha',1,'Color',[1 0 0],'LineStyle','--',...
            'LineWidth',1);
        ylabel({'Current (A)'},'Interpreter','latex','Units', 'normalized','position',[-0.113092878602192,0.525806928450061,0]);
        xlabel({'SOC'},'Interpreter','latex','Units','normalized');
        set(h1,'CData',plot_data(:,:,partition_idx),'FaceColor','interp','FaceAlpha',1)
        shading interp;
        axes1.YAxis.MinorTickValues = -1:0.5:1;
        xlim([0,1]);
        ylim([-1 ,1]);
        set(axes1,'FontSize',fontSize,'TickLabelInterpreter','latex',...
            'XMinorTick','on','YMinorTick','on','XTick',0:0.25:1,'YTick',[-1  0 1],...
            'YTickLabel',{'-1C','0','1C'});
        if(partition_idx == 1)
            text('Units','normalized','BackgroundColor',[1 1 1],...
                'EdgeColor',[0 0 0],...
                'Interpreter','latex',...
                'String',strcat('age $\in$[0,1/',num2str(deglifePartitions_num),']'),...
                'Position',textpositions{partition_idx},...
                'Color',[0 0 0],'FontSize',fontSize);
        elseif(partition_idx == deglifePartitions_num)
            text('Units','normalized','BackgroundColor',[1 1 1],...
                'EdgeColor',[0 0 0],...
                'Interpreter','latex',...
                'String',strcat('age $\in$[',num2str(deglifePartitions_num-1),'/',num2str(deglifePartitions_num),',1]'),...
                'Position',textpositions{partition_idx},...
                'Color',[0 0 0],'FontSize',fontSize);
        else
            text('Units','normalized','BackgroundColor',[1 1 1],...
                'EdgeColor',[0 0 0],...
                'Interpreter','latex',...
                'String',strcat('age $\in$[',num2str(partition_idx-1),'/',num2str(deglifePartitions_num),',',num2str(partition_idx),'/',num2str(deglifePartitions_num),']'),...
                'Position',textpositions{partition_idx},...
                'Color',[0 0 0],'FontSize',fontSize);
        end
        colormap parula;
        caxis([data_min data_max]);
    end
    labelTicks = round(linspace(round(data_min,1),round(data_max,1),5),1); % 1 hour
    h = colorbar('peer',axes1,'Position',[0.888 0.064 0.0266666666666666 0.816],...
        'TickLabelInterpreter','latex',...
        'FontSize',fontSize,'Ticks',labelTicks,'Limits',[round(data_min,1),round(data_max,1)]); % 1 hour
    caxis([round(data_min,1) round(data_max,1)]);
    colormap parula;
    ylabel(h, 'Mean value of SOC estimate','Interpreter','latex','FontSize',fontSize,'Position',[2.91555871963501 0.500000381469727 0]);
    fig_pos = figure_1.PaperPosition;
    figure_1.PaperSize = [fig_pos(3) fig_pos(4)];
    if(savefiles)
        saveas(gcf,strcat(filename,'.fig'));
        saveas(gcf,strcat(filename,'.eps'),'epsc');
    end
end

plot_life_partitions_soc_kp1_relVar = 1;
if(show_plots&&plot_life_partitions_soc_kp1_relVar)
    filename = 'socRelVar';
    partition_idxs = 1:deglifePartitions_num;
    positions = {[0.0924678097197633 0.246165556113454 0.20936667244291 0.615405872457975],...
        [0.373264911169039 0.238165556113454 0.20936667244291 0.623405872457975],...
        [0.654062012618313 0.238165556113454 0.209366672442911 0.627405872457975]};
    
    textpositions = {[0.232218011771345 1.10249990460461 0],...
        [0.181262597758605 1.10233123498622 0],...
        [0.24495686527453 1.09583826559003 0]...
        };
    
    data_min = min(relVar_soc_kp1_3c(:));
    data_max = round(max(relVar_soc_kp1_3c(:))*10^6)*10^-6+0.2*10^-6;
    
    y_axis_grid = cell_pow_set/(cell_1C_capacityInAh*cell_nominalVoltage);
    x_axis_grid = unique([0,soc_grid_boundaries,1]);
    x_axis_grid_mean = zeros(length(x_axis_grid)-1,1);
    for idx = 1:length(x_axis_grid)-1
        x_axis_grid_mean(idx) = 0.5*(x_axis_grid(idx) +x_axis_grid(idx+1));
    end
    
    [x,y] = meshgrid(x_axis_grid,y_axis_grid);
    
    plot_data = nan(length(y_axis_grid),length(x_axis_grid),deglifePartitions_num);
    for partition_idx = 1:deglifePartitions_num
        for soc_bin_idx = 1:z_num
            plot_data(:,soc_bin_idx+1,partition_idx) = relVar_soc_kp1_3c(soc_bin_idx,:,partition_idx);
        end
    end
    plot_data(:,z_num+2,:) = relVar_soc_kp1_3c(z_num,:,:);
    
    fontSize = 12;
    figure_xsize = 700;
    figure_ysize = 200;
    figure_position_offset = 50;
    figure_1 = figure('Color','w','Renderer','painters');
    set(gcf, 'Position',  [figure_position_offset,figure_position_offset,figure_position_offset+figure_xsize, figure_position_offset+figure_ysize]);
    
    for partition_idx = partition_idxs
        color_data_temp = cell(length(x_axis_grid)-1,1);
        for i = 1:length(x_axis_grid)-1
            color_data_temp{i} = [plot_data(:,i,partition_idx),plot_data(:,i,partition_idx)];
        end
        axes1 = axes('Parent',figure_1,...
            'Position',positions{partition_idx});
        box(axes1,'on');
        hold(axes1,'on');
        h1 = surf(x,y,0*plot_data(:,:,partition_idx));
        xline(cell_SOC_low,'Alpha',1,'Color',[1 0 0],'LineStyle','--',...
            'LineWidth',1);
        xline(cell_SOC_high,'Alpha',1,'Color',[1 0 0],'LineStyle','--',...
            'LineWidth',1);
        ylabel({'Current (A)'},'Interpreter','latex','Units', 'normalized','position',[-0.113092878602192,0.525806928450061,0]);
        xlabel({'SOC'},'Interpreter','latex','Units','normalized');
        set(h1,'CData',plot_data(:,:,partition_idx),'FaceColor','interp','FaceAlpha',1)
        shading interp;
        axes1.YAxis.MinorTickValues = -1:0.5:1;
        xlim([0,1]);
        ylim([-1 ,1]);
        set(axes1,'FontSize',fontSize,'TickLabelInterpreter','latex',...
            'XMinorTick','on','YMinorTick','on','XTick',0:0.25:1,'YTick',[-1  0 1],'Colormap',...
    [0.2422 0.1504 0.6603;0.249022413793103 0.162043103448276 0.699703448275862;0.255179310344828 0.175886206896552 0.735610344827586;0.261070689655172 0.18923275862069 0.771712068965517;0.266565517241379 0.202372413793103 0.807558620689655;0.271229310344828 0.216734482758621 0.84053275862069;0.274868965517241 0.232827586206897 0.868924137931034;0.277605172413793 0.25048275862069 0.892705172413793;0.279648275862069 0.268713793103448 0.912875862068966;0.280884482758621 0.287210344827586 0.930241379310345;0.2814 0.30551724137931 0.945258620689655;0.280979310344828 0.323718965517241 0.95853275862069;0.279834482758621 0.341727586206897 0.969979310344828;0.277341379310345 0.359722413793103 0.979772413793103;0.273220689655172 0.37808275862069 0.987313793103448;0.267181034482759 0.396808620689655 0.992296551724138;0.258255172413793 0.4159 0.995741379310345;0.244512068965517 0.435256896551724 0.998755172413793;0.225631034482759 0.455258620689655 0.998603448275862;0.205255172413793 0.475644827586207 0.992151724137931;0.188424137931034 0.494989655172414 0.985431034482759;0.181956896551724 0.513148275862069 0.976925862068965;0.178355172413793 0.530696551724138 0.966955172413793;0.176586206896552 0.547979310344828 0.953648275862069;0.171931034482759 0.564810344827586 0.939837931034483;0.160298275862069 0.581429310344828 0.928236206896552;0.150213793103448 0.597431034482759 0.916389655172414;0.14483275862069 0.612786206896552 0.905465517241379;0.1378 0.627975862068966 0.897120689655172;0.12705 0.64305 0.88985;0.115768965517241 0.657658620689655 0.88091724137931;0.104194827586207 0.671620689655173 0.868848275862069;0.0885517241379306 0.684586206896552 0.85383448275862;0.0644120689655166 0.696555172413793 0.836524137931034;0.0315034482758614 0.707541379310345 0.817486206896551;0.00700517241379276 0.71765 0.797413793103448;0.00161379310344835 0.726941379310345 0.776503448275862;0.0172982758620697 0.735489655172414 0.755003448275862;0.0564068965517253 0.743331034482759 0.732896551724137;0.0995206896551735 0.750529310344828 0.710324137931034;0.136372413793104 0.757417241379311 0.687386206896551;0.163874137931035 0.764405172413793 0.663996551724138;0.183996551724138 0.771568965517242 0.639710344827586;0.2008 0.778881034482759 0.613901724137931;0.220531034482759 0.785727586206897 0.586341379310344;0.246772413793104 0.791760344827586 0.556927586206896;0.281727586206897 0.796475862068965 0.525937931034483;0.322556896551724 0.79971724137931 0.493574137931034;0.363706896551724 0.801962068965517 0.459224137931035;0.405693103448276 0.803094827586207 0.422725862068966;0.4518 0.802172413793103 0.386341379310345;0.500149999999999 0.799563793103448 0.351106896551725;0.547479310344827 0.795651724137931 0.315524137931035;0.593596551724137 0.790684482758621 0.279425862068966;0.639079310344827 0.784520689655173 0.245341379310345;0.683524137931033 0.777260344827586 0.215362068965518;0.726144827586206 0.769472413793104 0.189917241379311;0.7945 0.7554 0.157;0.798159090909091 0.754604545454545 0.156236363636364;0.801818181818182 0.753809090909091 0.155472727272727;0.805477272727273 0.753013636363636 0.154709090909091;0.809054545454545 0.752245454545455 0.1543;0.812618181818182 0.751481818181818 0.15395;0.816181818181818 0.750718181818182 0.1536;0.819677272727272 0.749954545454546 0.153522727272727;0.823145454545454 0.749190909090909 0.153554545454545;0.826613636363636 0.748427272727273 0.153586363636364;0.830063636363636 0.747663636363637 0.153781818181818;0.833499999999999 0.7469 0.1541;0.836936363636363 0.746136363636364 0.154418181818182;0.840345454545454 0.7454 0.154845454545454;0.843718181818181 0.7447 0.155418181818182;0.847090909090908 0.744 0.155990909090909;0.850454545454545 0.7433 0.156609090909091;0.853795454545454 0.7426 0.157340909090909;0.857136363636363 0.7419 0.158072727272727;0.860468181818181 0.741204545454546 0.158827272727272;0.863745454545453 0.740536363636364 0.159718181818182;0.867022727272726 0.739868181818182 0.160609090909091;0.870299999999999 0.7392 0.1615;0.873513636363635 0.738563636363637 0.162613636363636;0.876727272727271 0.737927272727273 0.163727272727272;0.879940909090908 0.737290909090909 0.164840909090909;0.883099999999999 0.736681818181818 0.166227272727272;0.886249999999999 0.736077272727273 0.16765909090909;0.889399999999999 0.735472727272728 0.169090909090908;0.892504545454544 0.734913636363637 0.170727272727272;0.895590909090908 0.734372727272728 0.172445454545454;0.898677272727271 0.733831818181818 0.174163636363635;0.901690909090908 0.733327272727273 0.176099999999999;0.904649999999998 0.73285 0.178199999999999;0.907609090909089 0.732372727272728 0.180299999999999;0.910540909090907 0.731922727272727 0.182522727272726;0.913436363636362 0.731509090909091 0.18490909090909;0.916331818181816 0.731095454545455 0.187295454545453;0.919199999999998 0.730709090909091 0.189754545454544;0.921999999999998 0.730390909090909 0.192395454545453;0.924799999999998 0.730072727272727 0.195036363636362;0.927586363636362 0.729763636363637 0.197699999999998;0.930290909090907 0.729509090909091 0.200499999999998;0.932995454545453 0.729254545454546 0.203299999999998;0.935699999999998 0.729 0.206099999999998;0.938340909090907 0.728840909090909 0.208963636363634;0.940981818181816 0.728681818181818 0.211827272727271;0.943622727272725 0.728522727272727 0.214690909090907;0.946263636363634 0.728472727272727 0.217445454545452;0.948904545454543 0.728440909090909 0.220181818181816;0.951545454545452 0.728409090909091 0.22291818181818;0.954186363636361 0.728422727272727 0.225404545454543;0.956827272727271 0.728454545454545 0.227790909090907;0.95946818181818 0.728486363636364 0.230177272727271;0.962109090909089 0.728627272727273 0.232309090909089;0.964749999999998 0.72885 0.234249999999998;0.967390909090907 0.729072727272727 0.236190909090907;0.970004545454543 0.729363636363636 0.237913636363635;0.972581818181816 0.729745454545454 0.239345454545453;0.975159090909088 0.730127272727272 0.240777272727271;0.977654545454543 0.730636363636363 0.242054545454545;0.979945454545452 0.731463636363636 0.242945454545454;0.982236363636361 0.732290909090908 0.243836363636363;0.984463636363634 0.73315909090909 0.244522727272728;0.986309090909089 0.734272727272726 0.243981818181819;0.988154545454544 0.735386363636362 0.24344090909091;0.989999999999998 0.736499999999999 0.242900000000001;0.991463636363635 0.737836363636362 0.241786363636365;0.992927272727271 0.739172727272726 0.240672727272729;0.994390909090907 0.740509090909089 0.239559090909092;0.995145454545454 0.742090909090907 0.238227272727274;0.995781818181817 0.743713636363635 0.236859090909092;0.996418181818181 0.745336363636362 0.235490909090911;0.996713636363636 0.747049999999998 0.234145454545456;0.996872727272727 0.748799999999998 0.232809090909093;0.997031818181818 0.750549999999998 0.231472727272729;0.997118181818182 0.75231818181818 0.230136363636365;0.99715 0.754099999999998 0.228800000000002;0.997181818181818 0.755881818181816 0.227463636363638;0.997186363636364 0.75767727272727 0.226113636363638;0.997154545454546 0.759490909090907 0.224745454545456;0.997122727272727 0.761304545454543 0.223377272727275;0.997081818181818 0.763118181818179 0.222009090909093;0.997018181818182 0.764931818181816 0.220640909090911;0.996954545454546 0.766745454545452 0.219272727272729;0.996886363636364 0.768559090909088 0.217904545454547;0.996790909090909 0.770372727272725 0.216536363636366;0.996695454545455 0.772186363636361 0.215168181818184;0.9966 0.773999999999997 0.213800000000002;0.996472727272728 0.775845454545452 0.21243181818182;0.996345454545455 0.777690909090906 0.211063636363638;0.996218181818182 0.779536363636361 0.209695454545457;0.996063636363637 0.781381818181815 0.208354545454547;0.995904545454546 0.78322727272727 0.207018181818184;0.995745454545455 0.785072727272724 0.20568181818182;0.995518181818182 0.786940909090906 0.204368181818184;0.995263636363637 0.788818181818179 0.203063636363638;0.995009090909091 0.790695454545452 0.201759090909093;0.994700000000001 0.792572727272724 0.200509090909093;0.994350000000001 0.794449999999997 0.199300000000002;0.994000000000001 0.79632727272727 0.198090909090911;0.993595454545455 0.798218181818179 0.196922727272729;0.993118181818183 0.80012727272727 0.195809090909093;0.99264090909091 0.80203636363636 0.194695454545456;0.992145454545455 0.803954545454542 0.193600000000002;0.991604545454546 0.805895454545451 0.192550000000002;0.991063636363637 0.80783636363636 0.191500000000002;0.990504545454547 0.809777272727269 0.190459090909093;0.989836363636365 0.811718181818178 0.189472727272729;0.989168181818183 0.813659090909088 0.188486363636365;0.988500000000001 0.815599999999997 0.187500000000002;0.987736363636365 0.817572727272724 0.186577272727274;0.986972727272729 0.819545454545451 0.185654545454547;0.986209090909092 0.821518181818178 0.18473181818182;0.985390909090911 0.823490909090905 0.183809090909093;0.984563636363638 0.825463636363633 0.182886363636365;0.983736363636365 0.82743636363636 0.181963636363638;0.982863636363638 0.829409090909087 0.181018181818184;0.981972727272729 0.831381818181814 0.180063636363638;0.98108181818182 0.833354545454542 0.179109090909093;0.980172727272729 0.835327272727269 0.178154545454547;0.979250000000002 0.837299999999996 0.177200000000002;0.978327272727275 0.839272727272723 0.176245454545456;0.977390909090911 0.841259090909087 0.175277272727275;0.976436363636366 0.843263636363632 0.174290909090911;0.97548181818182 0.845268181818178 0.173304545454547;0.974545454545456 0.847263636363632 0.172318181818184;0.973654545454547 0.84923636363636 0.17133181818182;0.972763636363638 0.851209090909087 0.170345454545457;0.97188181818182 0.853181818181814 0.169363636363638;0.971054545454547 0.855154545454541 0.168409090909093;0.970227272727274 0.857127272727269 0.167454545454547;0.969400000000002 0.859099999999996 0.166500000000002;0.968668181818183 0.861104545454541 0.165577272727275;0.967936363636365 0.863109090909086 0.164654545454547;0.967204545454547 0.865113636363632 0.16373181818182;0.966554545454547 0.867090909090905 0.162836363636366;0.965918181818183 0.869063636363632 0.161945454545457;0.96528181818182 0.871036363636359 0.161054545454547;0.964713636363638 0.873009090909087 0.160209090909093;0.964172727272729 0.874981818181814 0.15938181818182;0.963631818181819 0.876954545454541 0.158554545454547;0.963127272727274 0.878927272727268 0.157745454545456;0.962650000000001 0.880899999999995 0.156950000000002;0.962172727272728 0.882872727272723 0.156154545454547;0.961750000000001 0.88484545454545 0.155359090909093;0.961400000000001 0.886818181818177 0.154563636363638;0.961050000000001 0.888790909090904 0.153768181818184;0.960736363636364 0.890754545454541 0.152972727272729;0.960513636363637 0.89269545454545 0.152177272727275;0.96029090909091 0.894636363636359 0.15138181818182;0.960077272727273 0.896572727272723 0.150577272727275;0.959918181818182 0.898481818181813 0.149718181818184;0.959759090909091 0.900390909090904 0.148859090909093;0.9596 0.902299999999995 0.148000000000002;0.959568181818182 0.904240909090904 0.147045454545457;0.959536363636364 0.906181818181813 0.146090909090911;0.959504545454546 0.908122727272722 0.145136363636366;0.959554545454545 0.910009090909086 0.144127272727275;0.959618181818182 0.911886363636359 0.143109090909094;0.959681818181818 0.913763636363631 0.142090909090912;0.959790909090909 0.915663636363631 0.140981818181821;0.959918181818181 0.917572727272722 0.139836363636367;0.960045454545454 0.919481818181813 0.138690909090912;0.960227272727272 0.921372727272722 0.137509090909094;0.960449999999999 0.923249999999995 0.136300000000003;0.960672727272727 0.925127272727268 0.135090909090912;0.960936363636363 0.926990909090904 0.133854545454549;0.961254545454545 0.928836363636359 0.132581818181822;0.961572727272726 0.930681818181813 0.131309090909094;0.961899999999999 0.932536363636359 0.130009090909095;0.962249999999999 0.934413636363631 0.128640909090913;0.962599999999999 0.936290909090904 0.127272727272731;0.96295909090909 0.938163636363631 0.125895454545459;0.963372727272726 0.940009090909086 0.12446363636364;0.963786363636362 0.94185454545454 0.123031818181822;0.964199999999999 0.943699999999995 0.121600000000004;0.964677272727271 0.945513636363631 0.120072727272732;0.965154545454544 0.947327272727267 0.118545454545459;0.965631818181817 0.949140909090904 0.117018181818186;0.966163636363635 0.950981818181813 0.115381818181823;0.966704545454544 0.952827272727267 0.113727272727278;0.967245454545453 0.954672727272722 0.112072727272732;0.967809090909089 0.956495454545449 0.110350000000005;0.968381818181816 0.958309090909085 0.108600000000005;0.968954545454544 0.960122727272722 0.106850000000005;0.969545454545453 0.96195454545454 0.105009090909097;0.970149999999998 0.963799999999994 0.103100000000006;0.970754545454544 0.965645454545449 0.101190909090915;0.971359090909089 0.967477272727267 0.0992409090909152;0.971963636363634 0.969290909090904 0.0972363636363698;0.97256818181818 0.97110454545454 0.0952318181818244;0.973172727272725 0.972927272727267 0.0932000000000065;0.973777272727271 0.974772727272721 0.0911000000000066;0.974381818181816 0.976618181818176 0.0890000000000066;0.974990909090907 0.978459090909085 0.0868954545454613;0.975627272727271 0.980272727272722 0.0847636363636431;0.976263636363634 0.982086363636358 0.082631818181825;0.976899999999998 0.983899999999994 0.0805000000000069],...
    'YTickLabel',{'-1C','0','1C'});
        if(partition_idx == 1)
            text('Units','normalized','BackgroundColor',[1 1 1],...
                'EdgeColor',[0 0 0],...
                'Interpreter','latex',...
                'String',strcat('age $\in$[0,1/',num2str(deglifePartitions_num),']'),...
                'Position',textpositions{partition_idx},...
                'Color',[0 0 0],'FontSize',fontSize);
        elseif(partition_idx == deglifePartitions_num)
            text('Units','normalized','BackgroundColor',[1 1 1],...
                'EdgeColor',[0 0 0],...
                'Interpreter','latex',...
                'String',strcat('age $\in$[',num2str(deglifePartitions_num-1),'/',num2str(deglifePartitions_num),',1]'),...
                'Position',textpositions{partition_idx},...
                'Color',[0 0 0],'FontSize',fontSize);
        else
            text('Units','normalized','BackgroundColor',[1 1 1],...
                'EdgeColor',[0 0 0],...
                'Interpreter','latex',...
                'String',strcat('age $\in$[',num2str(partition_idx-1),'/',num2str(deglifePartitions_num),',',num2str(partition_idx),'/',num2str(deglifePartitions_num),']'),...
                'Position',textpositions{partition_idx},...
                'Color',[0 0 0],'FontSize',fontSize);
        end
%         colormap parula;
        caxis([data_min data_max]);
    end
    precisionDigits_plot = 2;
    labelTicks = round(linspace(round(data_min,precisionDigits_plot),round(data_max,precisionDigits_plot),5),precisionDigits_plot); % 1 hour
    h = colorbar('peer',axes1,'Position',[0.888 0.064 0.0266666666666666 0.816],...
        'TickLabelInterpreter','latex','Colormap',...
    [0.2422 0.1504 0.6603;0.249022413793103 0.162043103448276 0.699703448275862;0.255179310344828 0.175886206896552 0.735610344827586;0.261070689655172 0.18923275862069 0.771712068965517;0.266565517241379 0.202372413793103 0.807558620689655;0.271229310344828 0.216734482758621 0.84053275862069;0.274868965517241 0.232827586206897 0.868924137931034;0.277605172413793 0.25048275862069 0.892705172413793;0.279648275862069 0.268713793103448 0.912875862068966;0.280884482758621 0.287210344827586 0.930241379310345;0.2814 0.30551724137931 0.945258620689655;0.280979310344828 0.323718965517241 0.95853275862069;0.279834482758621 0.341727586206897 0.969979310344828;0.277341379310345 0.359722413793103 0.979772413793103;0.273220689655172 0.37808275862069 0.987313793103448;0.267181034482759 0.396808620689655 0.992296551724138;0.258255172413793 0.4159 0.995741379310345;0.244512068965517 0.435256896551724 0.998755172413793;0.225631034482759 0.455258620689655 0.998603448275862;0.205255172413793 0.475644827586207 0.992151724137931;0.188424137931034 0.494989655172414 0.985431034482759;0.181956896551724 0.513148275862069 0.976925862068965;0.178355172413793 0.530696551724138 0.966955172413793;0.176586206896552 0.547979310344828 0.953648275862069;0.171931034482759 0.564810344827586 0.939837931034483;0.160298275862069 0.581429310344828 0.928236206896552;0.150213793103448 0.597431034482759 0.916389655172414;0.14483275862069 0.612786206896552 0.905465517241379;0.1378 0.627975862068966 0.897120689655172;0.12705 0.64305 0.88985;0.115768965517241 0.657658620689655 0.88091724137931;0.104194827586207 0.671620689655173 0.868848275862069;0.0885517241379306 0.684586206896552 0.85383448275862;0.0644120689655166 0.696555172413793 0.836524137931034;0.0315034482758614 0.707541379310345 0.817486206896551;0.00700517241379276 0.71765 0.797413793103448;0.00161379310344835 0.726941379310345 0.776503448275862;0.0172982758620697 0.735489655172414 0.755003448275862;0.0564068965517253 0.743331034482759 0.732896551724137;0.0995206896551735 0.750529310344828 0.710324137931034;0.136372413793104 0.757417241379311 0.687386206896551;0.163874137931035 0.764405172413793 0.663996551724138;0.183996551724138 0.771568965517242 0.639710344827586;0.2008 0.778881034482759 0.613901724137931;0.220531034482759 0.785727586206897 0.586341379310344;0.246772413793104 0.791760344827586 0.556927586206896;0.281727586206897 0.796475862068965 0.525937931034483;0.322556896551724 0.79971724137931 0.493574137931034;0.363706896551724 0.801962068965517 0.459224137931035;0.405693103448276 0.803094827586207 0.422725862068966;0.4518 0.802172413793103 0.386341379310345;0.500149999999999 0.799563793103448 0.351106896551725;0.547479310344827 0.795651724137931 0.315524137931035;0.593596551724137 0.790684482758621 0.279425862068966;0.639079310344827 0.784520689655173 0.245341379310345;0.683524137931033 0.777260344827586 0.215362068965518;0.726144827586206 0.769472413793104 0.189917241379311;0.7945 0.7554 0.157;0.798159090909091 0.754604545454545 0.156236363636364;0.801818181818182 0.753809090909091 0.155472727272727;0.805477272727273 0.753013636363636 0.154709090909091;0.809054545454545 0.752245454545455 0.1543;0.812618181818182 0.751481818181818 0.15395;0.816181818181818 0.750718181818182 0.1536;0.819677272727272 0.749954545454546 0.153522727272727;0.823145454545454 0.749190909090909 0.153554545454545;0.826613636363636 0.748427272727273 0.153586363636364;0.830063636363636 0.747663636363637 0.153781818181818;0.833499999999999 0.7469 0.1541;0.836936363636363 0.746136363636364 0.154418181818182;0.840345454545454 0.7454 0.154845454545454;0.843718181818181 0.7447 0.155418181818182;0.847090909090908 0.744 0.155990909090909;0.850454545454545 0.7433 0.156609090909091;0.853795454545454 0.7426 0.157340909090909;0.857136363636363 0.7419 0.158072727272727;0.860468181818181 0.741204545454546 0.158827272727272;0.863745454545453 0.740536363636364 0.159718181818182;0.867022727272726 0.739868181818182 0.160609090909091;0.870299999999999 0.7392 0.1615;0.873513636363635 0.738563636363637 0.162613636363636;0.876727272727271 0.737927272727273 0.163727272727272;0.879940909090908 0.737290909090909 0.164840909090909;0.883099999999999 0.736681818181818 0.166227272727272;0.886249999999999 0.736077272727273 0.16765909090909;0.889399999999999 0.735472727272728 0.169090909090908;0.892504545454544 0.734913636363637 0.170727272727272;0.895590909090908 0.734372727272728 0.172445454545454;0.898677272727271 0.733831818181818 0.174163636363635;0.901690909090908 0.733327272727273 0.176099999999999;0.904649999999998 0.73285 0.178199999999999;0.907609090909089 0.732372727272728 0.180299999999999;0.910540909090907 0.731922727272727 0.182522727272726;0.913436363636362 0.731509090909091 0.18490909090909;0.916331818181816 0.731095454545455 0.187295454545453;0.919199999999998 0.730709090909091 0.189754545454544;0.921999999999998 0.730390909090909 0.192395454545453;0.924799999999998 0.730072727272727 0.195036363636362;0.927586363636362 0.729763636363637 0.197699999999998;0.930290909090907 0.729509090909091 0.200499999999998;0.932995454545453 0.729254545454546 0.203299999999998;0.935699999999998 0.729 0.206099999999998;0.938340909090907 0.728840909090909 0.208963636363634;0.940981818181816 0.728681818181818 0.211827272727271;0.943622727272725 0.728522727272727 0.214690909090907;0.946263636363634 0.728472727272727 0.217445454545452;0.948904545454543 0.728440909090909 0.220181818181816;0.951545454545452 0.728409090909091 0.22291818181818;0.954186363636361 0.728422727272727 0.225404545454543;0.956827272727271 0.728454545454545 0.227790909090907;0.95946818181818 0.728486363636364 0.230177272727271;0.962109090909089 0.728627272727273 0.232309090909089;0.964749999999998 0.72885 0.234249999999998;0.967390909090907 0.729072727272727 0.236190909090907;0.970004545454543 0.729363636363636 0.237913636363635;0.972581818181816 0.729745454545454 0.239345454545453;0.975159090909088 0.730127272727272 0.240777272727271;0.977654545454543 0.730636363636363 0.242054545454545;0.979945454545452 0.731463636363636 0.242945454545454;0.982236363636361 0.732290909090908 0.243836363636363;0.984463636363634 0.73315909090909 0.244522727272728;0.986309090909089 0.734272727272726 0.243981818181819;0.988154545454544 0.735386363636362 0.24344090909091;0.989999999999998 0.736499999999999 0.242900000000001;0.991463636363635 0.737836363636362 0.241786363636365;0.992927272727271 0.739172727272726 0.240672727272729;0.994390909090907 0.740509090909089 0.239559090909092;0.995145454545454 0.742090909090907 0.238227272727274;0.995781818181817 0.743713636363635 0.236859090909092;0.996418181818181 0.745336363636362 0.235490909090911;0.996713636363636 0.747049999999998 0.234145454545456;0.996872727272727 0.748799999999998 0.232809090909093;0.997031818181818 0.750549999999998 0.231472727272729;0.997118181818182 0.75231818181818 0.230136363636365;0.99715 0.754099999999998 0.228800000000002;0.997181818181818 0.755881818181816 0.227463636363638;0.997186363636364 0.75767727272727 0.226113636363638;0.997154545454546 0.759490909090907 0.224745454545456;0.997122727272727 0.761304545454543 0.223377272727275;0.997081818181818 0.763118181818179 0.222009090909093;0.997018181818182 0.764931818181816 0.220640909090911;0.996954545454546 0.766745454545452 0.219272727272729;0.996886363636364 0.768559090909088 0.217904545454547;0.996790909090909 0.770372727272725 0.216536363636366;0.996695454545455 0.772186363636361 0.215168181818184;0.9966 0.773999999999997 0.213800000000002;0.996472727272728 0.775845454545452 0.21243181818182;0.996345454545455 0.777690909090906 0.211063636363638;0.996218181818182 0.779536363636361 0.209695454545457;0.996063636363637 0.781381818181815 0.208354545454547;0.995904545454546 0.78322727272727 0.207018181818184;0.995745454545455 0.785072727272724 0.20568181818182;0.995518181818182 0.786940909090906 0.204368181818184;0.995263636363637 0.788818181818179 0.203063636363638;0.995009090909091 0.790695454545452 0.201759090909093;0.994700000000001 0.792572727272724 0.200509090909093;0.994350000000001 0.794449999999997 0.199300000000002;0.994000000000001 0.79632727272727 0.198090909090911;0.993595454545455 0.798218181818179 0.196922727272729;0.993118181818183 0.80012727272727 0.195809090909093;0.99264090909091 0.80203636363636 0.194695454545456;0.992145454545455 0.803954545454542 0.193600000000002;0.991604545454546 0.805895454545451 0.192550000000002;0.991063636363637 0.80783636363636 0.191500000000002;0.990504545454547 0.809777272727269 0.190459090909093;0.989836363636365 0.811718181818178 0.189472727272729;0.989168181818183 0.813659090909088 0.188486363636365;0.988500000000001 0.815599999999997 0.187500000000002;0.987736363636365 0.817572727272724 0.186577272727274;0.986972727272729 0.819545454545451 0.185654545454547;0.986209090909092 0.821518181818178 0.18473181818182;0.985390909090911 0.823490909090905 0.183809090909093;0.984563636363638 0.825463636363633 0.182886363636365;0.983736363636365 0.82743636363636 0.181963636363638;0.982863636363638 0.829409090909087 0.181018181818184;0.981972727272729 0.831381818181814 0.180063636363638;0.98108181818182 0.833354545454542 0.179109090909093;0.980172727272729 0.835327272727269 0.178154545454547;0.979250000000002 0.837299999999996 0.177200000000002;0.978327272727275 0.839272727272723 0.176245454545456;0.977390909090911 0.841259090909087 0.175277272727275;0.976436363636366 0.843263636363632 0.174290909090911;0.97548181818182 0.845268181818178 0.173304545454547;0.974545454545456 0.847263636363632 0.172318181818184;0.973654545454547 0.84923636363636 0.17133181818182;0.972763636363638 0.851209090909087 0.170345454545457;0.97188181818182 0.853181818181814 0.169363636363638;0.971054545454547 0.855154545454541 0.168409090909093;0.970227272727274 0.857127272727269 0.167454545454547;0.969400000000002 0.859099999999996 0.166500000000002;0.968668181818183 0.861104545454541 0.165577272727275;0.967936363636365 0.863109090909086 0.164654545454547;0.967204545454547 0.865113636363632 0.16373181818182;0.966554545454547 0.867090909090905 0.162836363636366;0.965918181818183 0.869063636363632 0.161945454545457;0.96528181818182 0.871036363636359 0.161054545454547;0.964713636363638 0.873009090909087 0.160209090909093;0.964172727272729 0.874981818181814 0.15938181818182;0.963631818181819 0.876954545454541 0.158554545454547;0.963127272727274 0.878927272727268 0.157745454545456;0.962650000000001 0.880899999999995 0.156950000000002;0.962172727272728 0.882872727272723 0.156154545454547;0.961750000000001 0.88484545454545 0.155359090909093;0.961400000000001 0.886818181818177 0.154563636363638;0.961050000000001 0.888790909090904 0.153768181818184;0.960736363636364 0.890754545454541 0.152972727272729;0.960513636363637 0.89269545454545 0.152177272727275;0.96029090909091 0.894636363636359 0.15138181818182;0.960077272727273 0.896572727272723 0.150577272727275;0.959918181818182 0.898481818181813 0.149718181818184;0.959759090909091 0.900390909090904 0.148859090909093;0.9596 0.902299999999995 0.148000000000002;0.959568181818182 0.904240909090904 0.147045454545457;0.959536363636364 0.906181818181813 0.146090909090911;0.959504545454546 0.908122727272722 0.145136363636366;0.959554545454545 0.910009090909086 0.144127272727275;0.959618181818182 0.911886363636359 0.143109090909094;0.959681818181818 0.913763636363631 0.142090909090912;0.959790909090909 0.915663636363631 0.140981818181821;0.959918181818181 0.917572727272722 0.139836363636367;0.960045454545454 0.919481818181813 0.138690909090912;0.960227272727272 0.921372727272722 0.137509090909094;0.960449999999999 0.923249999999995 0.136300000000003;0.960672727272727 0.925127272727268 0.135090909090912;0.960936363636363 0.926990909090904 0.133854545454549;0.961254545454545 0.928836363636359 0.132581818181822;0.961572727272726 0.930681818181813 0.131309090909094;0.961899999999999 0.932536363636359 0.130009090909095;0.962249999999999 0.934413636363631 0.128640909090913;0.962599999999999 0.936290909090904 0.127272727272731;0.96295909090909 0.938163636363631 0.125895454545459;0.963372727272726 0.940009090909086 0.12446363636364;0.963786363636362 0.94185454545454 0.123031818181822;0.964199999999999 0.943699999999995 0.121600000000004;0.964677272727271 0.945513636363631 0.120072727272732;0.965154545454544 0.947327272727267 0.118545454545459;0.965631818181817 0.949140909090904 0.117018181818186;0.966163636363635 0.950981818181813 0.115381818181823;0.966704545454544 0.952827272727267 0.113727272727278;0.967245454545453 0.954672727272722 0.112072727272732;0.967809090909089 0.956495454545449 0.110350000000005;0.968381818181816 0.958309090909085 0.108600000000005;0.968954545454544 0.960122727272722 0.106850000000005;0.969545454545453 0.96195454545454 0.105009090909097;0.970149999999998 0.963799999999994 0.103100000000006;0.970754545454544 0.965645454545449 0.101190909090915;0.971359090909089 0.967477272727267 0.0992409090909152;0.971963636363634 0.969290909090904 0.0972363636363698;0.97256818181818 0.97110454545454 0.0952318181818244;0.973172727272725 0.972927272727267 0.0932000000000065;0.973777272727271 0.974772727272721 0.0911000000000066;0.974381818181816 0.976618181818176 0.0890000000000066;0.974990909090907 0.978459090909085 0.0868954545454613;0.975627272727271 0.980272727272722 0.0847636363636431;0.976263636363634 0.982086363636358 0.082631818181825;0.976899999999998 0.983899999999994 0.0805000000000069],...
   'FontSize',fontSize,'Ticks',labelTicks,'Limits',[round(data_min,precisionDigits_plot),round(data_max,precisionDigits_plot)]); % 1 hour
    caxis([round(data_min,precisionDigits_plot) round(data_max,precisionDigits_plot)]);
%     colormap parula;
    ylabel(h, 'Rel. var. of SOC estimate ($\%$)','Interpreter','latex','FontSize',fontSize);
    fig_pos = figure_1.PaperPosition;
    figure_1.PaperSize = [fig_pos(3) fig_pos(4)];
    if(savefiles)
        saveas(gcf,strcat(filename,'.fig'));
        saveas(gcf,strcat(filename,'.eps'),'epsc');
    end
end

plot_life_partitions_energy_loss_percentage_NRMSE = 0;
if(plot_life_partitions_energy_loss_percentage_NRMSE)
    filename = 'enLossNRMSE';
    partition_idxs = 1:deglifePartitions_num;
    positions = {[0.0924678097197633 0.246165556113454 0.20936667244291 0.615405872457975],...
        [0.373264911169039 0.238165556113454 0.20936667244291 0.623405872457975],...
        [0.654062012618313 0.238165556113454 0.209366672442911 0.627405872457975]};
    
    textpositions = {[0.232218011771345 1.10249990460461 0],...
        [0.181262597758605 1.10233123498622 0],...
        [0.24495686527453 1.09583826559003 0]...
        };
    
    y_axis_grid = cell_pow_set/(cell_1C_capacityInAh*cell_nominalVoltage);
    x_axis_grid = unique([0,soc_grid_boundaries,1]);
    x_axis_grid_mean = zeros(length(x_axis_grid)-1,1);
    for idx = 1:length(x_axis_grid)-1
        x_axis_grid_mean(idx) = 0.5*(x_axis_grid(idx) +x_axis_grid(idx+1));
    end
    
    [x,y] = meshgrid(x_axis_grid,y_axis_grid);
    
    plot_data = nan(length(y_axis_grid),length(x_axis_grid),deglifePartitions_num);
    for partition_idx = 1:deglifePartitions_num
        for soc_bin_idx = 1:z_num
            plot_data(:,soc_bin_idx+1,partition_idx) = nrmse_energyLoss_percentage(soc_bin_idx,:,partition_idx);
        end
    end
    plot_data(:,z_num+2,:) = plot_data(:,z_num+1,:);
    
    data_min = min(plot_data(:));
    data_max = max(plot_data(:));
    
    fontSize = 12;
    figure_xsize = 700;
    figure_ysize = 200;
    figure_position_offset = 50;
    figure_1 = figure('Color','w','Renderer','painters');
    set(gcf, 'Position',  [figure_position_offset,figure_position_offset,figure_position_offset+figure_xsize, figure_position_offset+figure_ysize]);
    
    for partition_idx = partition_idxs
        axes1 = axes('Parent',figure_1,...
            'Position',positions{partition_idx});
        box(axes1,'on');
        hold(axes1,'on');
        h1 = surf(x,y,plot_data(:,:,partition_idx));
        xline(cell_SOC_low,'Alpha',1,'Color',[1 0 0],'LineStyle','--',...
            'LineWidth',1);
        xline(cell_SOC_high,'Alpha',1,'Color',[1 0 0],'LineStyle','--',...
            'LineWidth',1);
        ylabel({'Current (A)'},'Interpreter','latex','Units', 'normalized','position',[-0.113092878602192,0.525806928450061,0]);
        xlabel({'SOC'},'Interpreter','latex','Units','normalized');
        set(h1, 'FaceColor','interp','FaceAlpha',1)
        %         shading interp;
        axes1.YAxis.MinorTickValues = -1:0.5:1;
        xlim([0,1]);
        ylim([-1 ,1]);
        set(axes1,'FontSize',fontSize,'TickLabelInterpreter','latex',...
            'XMinorTick','on','YMinorTick','on','XTick',0:0.25:1,'YTick',[-1  0 1],...
            'YTickLabel',{'-1C','0','1C'});
        if(partition_idx == 1)
            text('Units','normalized','BackgroundColor',[1 1 1],...
                'EdgeColor',[0 0 0],...
                'Interpreter','latex',...
                'String',strcat('age $\in$[0,1/',num2str(deglifePartitions_num),']'),...
                'Position',textpositions{partition_idx},...
                'Color',[0 0 0],'FontSize',fontSize);
        elseif(partition_idx == deglifePartitions_num)
            text('Units','normalized','BackgroundColor',[1 1 1],...
                'EdgeColor',[0 0 0],...
                'Interpreter','latex',...
                'String',strcat('age $\in$[',num2str(deglifePartitions_num-1),'/',num2str(deglifePartitions_num),',1]'),...
                'Position',textpositions{partition_idx},...
                'Color',[0 0 0],'FontSize',fontSize);
        else
            text('Units','normalized','BackgroundColor',[1 1 1],...
                'EdgeColor',[0 0 0],...
                'Interpreter','latex',...
                'String',strcat('age $\in$[',num2str(partition_idx-1),'/',num2str(deglifePartitions_num),',',num2str(partition_idx),'/',num2str(deglifePartitions_num),']'),...
                'Position',textpositions{partition_idx},...
                'Color',[0 0 0],'FontSize',fontSize);
        end
        colormap jet;
        caxis([data_min data_max]);
    end
    labelTicks = round(linspace(data_min,data_max,5),2);
    h = colorbar('Position',[0.888 0.068 0.0266666666666669 0.892],...
        'TickLabelInterpreter','latex','FontSize',fontSize,'Ticks',labelTicks,'Limits',[0,round(data_max,2)]);
    caxis([data_min data_max]);
    ylabel(h, 'NRMSE of en. loss estimate ($\%$)','Interpreter','latex','FontSize',fontSize);
    colormap jet;
    fig_pos = figure_1.PaperPosition;
    figure_1.PaperSize = [fig_pos(3) fig_pos(4)];
    if(savefiles)
        saveas(gcf,strcat(filename,'.fig'));
        saveas(gcf,strcat(filename,'.eps'),'epsc');
    end
end

plot_life_partitions_soc_kp1_NRMSE = 0;
if(plot_life_partitions_soc_kp1_NRMSE)
    filename = 'socNRMSE';
    partition_idxs = 1:deglifePartitions_num;
    positions = {[0.0924678097197633 0.246165556113454 0.20936667244291 0.615405872457975],...
        [0.373264911169039 0.238165556113454 0.20936667244291 0.623405872457975],...
        [0.654062012618313 0.238165556113454 0.209366672442911 0.627405872457975]};
    
    textpositions = {[0.232218011771345 1.10249990460461 0],...
        [0.181262597758605 1.10233123498622 0],...
        [0.24495686527453 1.09583826559003 0]...
        };
    
    data_min = min(nrmse_soc_kp1_percentage(:));
    data_max = round(max(nrmse_soc_kp1_percentage(:))*10^6)*10^-6;
    
    y_axis_grid = cell_pow_set/(cell_1C_capacityInAh*cell_nominalVoltage);
    x_axis_grid = unique([0,soc_grid_boundaries,1]);
    x_axis_grid_mean = zeros(length(x_axis_grid)-1,1);
    for idx = 1:length(x_axis_grid)-1
        x_axis_grid_mean(idx) = 0.5*(x_axis_grid(idx) +x_axis_grid(idx+1));
    end
    
    [x,y] = meshgrid(x_axis_grid,y_axis_grid);
    
    plot_data = nan(length(y_axis_grid),length(x_axis_grid),deglifePartitions_num);
    for partition_idx = 1:deglifePartitions_num
        for soc_bin_idx = 1:z_num
            plot_data(:,soc_bin_idx+1,partition_idx) = nrmse_soc_kp1_percentage(soc_bin_idx,:,partition_idx);
        end
    end
    plot_data(:,z_num+2,:) = nrmse_soc_kp1_percentage(z_num,:,:);
    
    fontSize = 12;
    figure_xsize = 700;
    figure_ysize = 200;
    figure_position_offset = 50;
    figure_1 = figure('Color','w','Renderer','painters');
    set(gcf, 'Position',  [figure_position_offset,figure_position_offset,figure_position_offset+figure_xsize, figure_position_offset+figure_ysize]);
    
    for partition_idx = partition_idxs
        color_data_temp = cell(length(x_axis_grid)-1,1);
        for i = 1:length(x_axis_grid)-1
            color_data_temp{i} = [plot_data(:,i,partition_idx),plot_data(:,i,partition_idx)];
        end
        axes1 = axes('Parent',figure_1,...
            'Position',positions{partition_idx});
        box(axes1,'on');
        hold(axes1,'on');
        h1 = surf(x,y,plot_data(:,:,partition_idx));
        xline(cell_SOC_low,'Alpha',1,'Color',[1 0 0],'LineStyle','--',...
            'LineWidth',1);
        xline(cell_SOC_high,'Alpha',1,'Color',[1 0 0],'LineStyle','--',...
            'LineWidth',1);
        ylabel({'Current (A)'},'Interpreter','latex','Units', 'normalized','position',[-0.113092878602192,0.525806928450061,0]);
        xlabel({'SOC'},'Interpreter','latex','Units','normalized');
        set(h1, 'FaceColor','interp','FaceAlpha',1)
        %         shading interp;
        axes1.YAxis.MinorTickValues = -1:0.5:1;
        xlim([0,1]);
        ylim([-1 ,1]);
        set(axes1,'FontSize',fontSize,'TickLabelInterpreter','latex',...
            'XMinorTick','on','YMinorTick','on','XTick',0:0.25:1,'YTick',[-1  0 1],...
            'Colormap',...
            [0 0 0.515625;0 0 0.552364864864865;0 0 0.58910472972973;0 0 0.625844594594595;0 0 0.662584459459459;0 0 0.699324324324324;0 0 0.736064189189189;0 0 0.772804054054054;0 0 0.809543918918919;0 0 0.846283783783784;0 0 0.883023648648649;0 0 0.919763513513513;0 0 0.956503378378378;0 0 0.993243243243243;0 0.0299831081081081 1;0 0.066722972972973 1;0 0.103462837837838 1;0 0.140202702702703 1;0 0.176942567567568 1;0 0.213682432432433 1;0 0.250422297297298 1;0 0.287162162162163 1;0 0.323902027027027 1;0 0.360641891891892 1;0 0.397381756756757 1;0 0.434121621621622 1;0 0.470861486486487 1;0 0.507601351351352 1;0 0.544341216216217 1;0 0.581081081081082 1;0 0.617820945945947 1;0 0.654560810810812 1;0 0.691300675675677 1;0 0.728040540540541 1;0 0.764780405405406 1;0 0.801520270270271 1;0 0.859375 1;0 0.871432648401826 1;0 0.883490296803653 1;0 0.895547945205479 1;0 0.907605593607306 1;0 0.919663242009132 1;0 0.931720890410959 1;0 0.943778538812785 1;0 0.955836187214611 1;0 0.967893835616438 1;0 0.979951484018264 1;0 0.992009132420091 1;0.00406678082191725 1 0.995933219178083;0.0161244292237437 1 0.983875570776256;0.0281820776255701 1 0.97181792237443;0.0402397260273966 1 0.959760273972603;0.052297374429223 1 0.947702625570777;0.0643550228310494 1 0.935644977168951;0.0764126712328759 1 0.923587328767124;0.0884703196347023 1 0.911529680365298;0.100527968036529 1 0.899472031963471;0.112585616438355 1 0.887414383561645;0.124643264840182 1 0.875356735159818;0.136700913242008 1 0.863299086757992;0.148758561643834 1 0.851241438356166;0.160816210045661 1 0.839183789954339;0.172873858447487 1 0.827126141552513;0.184931506849314 1 0.815068493150686;0.19698915525114 1 0.80301084474886;0.209046803652967 1 0.790953196347033;0.221104452054793 1 0.778895547945207;0.23316210045662 1 0.76683789954338;0.245219748858446 1 0.754780251141554;0.257277397260272 1 0.742722602739728;0.269335045662099 1 0.730664954337901;0.281392694063925 1 0.718607305936075;0.293450342465752 1 0.706549657534248;0.305507990867578 1 0.694492009132422;0.317565639269405 1 0.682434360730595;0.329623287671231 1 0.670376712328769;0.341680936073057 1 0.658319063926943;0.353738584474884 1 0.646261415525116;0.36579623287671 1 0.63420376712329;0.377853881278537 1 0.622146118721463;0.389911529680363 1 0.610088470319637;0.40196917808219 1 0.59803082191781;0.414026826484016 1 0.585973173515984;0.426084474885843 1 0.573915525114157;0.438142123287669 1 0.561857876712331;0.450199771689495 1 0.549800228310505;0.462257420091322 1 0.537742579908678;0.474315068493148 1 0.525684931506852;0.486372716894975 1 0.513627283105025;0.498430365296801 1 0.501569634703199;0.510488013698628 1 0.489511986301372;0.522545662100454 1 0.477454337899546;0.53460331050228 1 0.46539668949772;0.546660958904107 1 0.453339041095893;0.558718607305933 1 0.441281392694067;0.57077625570776 1 0.42922374429224;0.582833904109586 1 0.417166095890414;0.594891552511413 1 0.405108447488587;0.606949200913239 1 0.393050799086761;0.619006849315066 1 0.380993150684934;0.631064497716892 1 0.368935502283108;0.643122146118718 1 0.356877853881282;0.655179794520545 1 0.344820205479455;0.667237442922371 1 0.332762557077629;0.679295091324198 1 0.320704908675802;0.691352739726024 1 0.308647260273976;0.703410388127851 1 0.296589611872149;0.715468036529677 1 0.284531963470323;0.727525684931503 1 0.272474315068497;0.73958333333333 1 0.26041666666667;0.751640981735156 1 0.248359018264844;0.763698630136983 1 0.236301369863017;0.775756278538809 1 0.224243721461191;0.787813926940636 1 0.212186073059364;0.799871575342462 1 0.200128424657538;0.811929223744289 1 0.188070776255711;0.823986872146115 1 0.176013127853885;0.836044520547941 1 0.163955479452059;0.848102168949768 1 0.151897831050232;0.860159817351594 1 0.139840182648406;0.872217465753421 1 0.127782534246579;0.884275114155247 1 0.115724885844753;0.896332762557074 1 0.103667237442926;0.9083904109589 1 0.0916095890410999;0.920448059360726 1 0.0795519406392735;0.932505707762553 1 0.0674942922374471;0.944563356164379 1 0.0554366438356206;0.956621004566206 1 0.0433789954337942;0.968678652968032 1 0.0313213470319678;0.980736301369859 1 0.0192636986301413;0.992793949771685 1 0.00720605022831489;1 0.995148401826488 0;1 0.983090753424662 0;1 0.971033105022836 0;1 0.958975456621009 0;1 0.946917808219183 0;1 0.934860159817356 0;1 0.92280251141553 0;1 0.910744863013703 0;1 0.898687214611877 0;1 0.886629566210051 0;1 0.874571917808224 0;1 0.862514269406398 0;1 0.850456621004571 0;1 0.838398972602745 0;1 0.826341324200918 0;1 0.814283675799092 0;1 0.802226027397265 0;1 0.790168378995439 0;1 0.778110730593613 0;1 0.766053082191786 0;1 0.75399543378996 0;1 0.741937785388133 0;1 0.729880136986307 0;1 0.71782248858448 0;1 0.705764840182654 0;1 0.693707191780828 0;1 0.681649543379001 0;1 0.669591894977175 0;1 0.657534246575348 0;1 0.645476598173522 0;1 0.633418949771695 0;1 0.621361301369869 0;1 0.609303652968042 0;1 0.597246004566216 0;1 0.58518835616439 0;1 0.573130707762563 0;1 0.561073059360737 0;1 0.54901541095891 0;1 0.536957762557084 0;1 0.524900114155257 0;1 0.512842465753431 0;1 0.500784817351605 0;1 0.488727168949778 0;1 0.476669520547952 0;1 0.464611872146125 0;1 0.452554223744299 0;1 0.440496575342472 0;1 0.428438926940646 0;1 0.416381278538819 0;1 0.404323630136993 0;1 0.392265981735167 0;1 0.38020833333334 0;1 0.368150684931514 0;1 0.356093036529687 0;1 0.344035388127861 0;1 0.331977739726034 0;1 0.319920091324208 0;1 0.307862442922382 0;1 0.295804794520555 0;1 0.283747146118729 0;1 0.271689497716902 0;1 0.259631849315076 0;1 0.247574200913249 0;1 0.235516552511423 0;1 0.223458904109596 0;1 0.21140125570777 0;1 0.199343607305944 0;1 0.187285958904117 0;1 0.175228310502291 0;1 0.163170662100464 0;1 0.151113013698638 0;1 0.139055365296811 0;1 0.126997716894985 0;1 0.114940068493159 0;1 0.102882420091332 0;1 0.0908247716895056 0;1 0.0787671232876792 0;1 0.0667094748858528 0;1 0.0546518264840263 0;1 0.0425941780821999 0;1 0.0305365296803735 0;1 0.018478881278547 0;1 0.00642123287672058 0;0.994363584474894 0 0;0.982305936073068 0 0;0.970248287671241 0 0;0.958190639269415 0 0;0.946132990867588 0 0;0.934075342465762 0 0;0.922017694063936 0 0;0.909960045662109 0 0;0.897902397260283 0 0;0.885844748858456 0 0;0.87378710045663 0 0;0.861729452054803 0 0;0.849671803652977 0 0;0.83761415525115 0 0;0.825556506849324 0 0;0.813498858447498 0 0;0.801441210045671 0 0;0.789383561643845 0 0;0.777325913242018 0 0;0.765268264840192 0 0;0.753210616438365 0 0;0.741152968036539 0 0;0.729095319634713 0 0;0.717037671232886 0 0;0.70498002283106 0 0;0.692922374429233 0 0;0.680864726027407 0 0;0.66880707762558 0 0;0.656749429223754 0 0;0.644691780821927 0 0;0.632634132420101 0 0;0.620576484018275 0 0;0.608518835616448 0 0;0.596461187214622 0 0;0.584403538812795 0 0;0.572345890410969 0 0;0.560288242009142 0 0;0.548230593607316 0 0;0.53617294520549 0 0;0.524115296803663 0 0;0.512057648401837 0 0;0.50000000000001 0 0],...
            'YTickLabel',{'-1C','0','1C'});
        if(partition_idx == 1)
            text('Units','normalized','BackgroundColor',[1 1 1],...
                'EdgeColor',[0 0 0],...
                'Interpreter','latex',...
                'String',strcat('age $\in$[0,1/',num2str(deglifePartitions_num),']'),...
                'Position',textpositions{partition_idx},...
                'Color',[0 0 0],'FontSize',fontSize);
        elseif(partition_idx == deglifePartitions_num)
            text('Units','normalized','BackgroundColor',[1 1 1],...
                'EdgeColor',[0 0 0],...
                'Interpreter','latex',...
                'String',strcat('age $\in$[',num2str(deglifePartitions_num-1),'/',num2str(deglifePartitions_num),',1]'),...
                'Position',textpositions{partition_idx},...
                'Color',[0 0 0],'FontSize',fontSize);
        else
            text('Units','normalized','BackgroundColor',[1 1 1],...
                'EdgeColor',[0 0 0],...
                'Interpreter','latex',...
                'String',strcat('age $\in$[',num2str(partition_idx-1),'/',num2str(deglifePartitions_num),',',num2str(partition_idx),'/',num2str(deglifePartitions_num),']'),...
                'Position',textpositions{partition_idx},...
                'Color',[0 0 0],'FontSize',fontSize);
        end
        %         colormap jet;
        caxis([round(data_min,1) round(data_max,1)]);
    end
    labelTicks = round(linspace(round(data_min,1),round(data_max,1),5),1); % 1 hour
    h = colorbar('peer',axes1,'Position',[0.888 0.064 0.0266666666666666 0.816],...
        'TickLabelInterpreter','latex','Colormap',...
        [0 0 0.515625;0 0 0.552364864864865;0 0 0.58910472972973;0 0 0.625844594594595;0 0 0.662584459459459;0 0 0.699324324324324;0 0 0.736064189189189;0 0 0.772804054054054;0 0 0.809543918918919;0 0 0.846283783783784;0 0 0.883023648648649;0 0 0.919763513513513;0 0 0.956503378378378;0 0 0.993243243243243;0 0.0299831081081081 1;0 0.066722972972973 1;0 0.103462837837838 1;0 0.140202702702703 1;0 0.176942567567568 1;0 0.213682432432433 1;0 0.250422297297298 1;0 0.287162162162163 1;0 0.323902027027027 1;0 0.360641891891892 1;0 0.397381756756757 1;0 0.434121621621622 1;0 0.470861486486487 1;0 0.507601351351352 1;0 0.544341216216217 1;0 0.581081081081082 1;0 0.617820945945947 1;0 0.654560810810812 1;0 0.691300675675677 1;0 0.728040540540541 1;0 0.764780405405406 1;0 0.801520270270271 1;0 0.859375 1;0 0.871432648401826 1;0 0.883490296803653 1;0 0.895547945205479 1;0 0.907605593607306 1;0 0.919663242009132 1;0 0.931720890410959 1;0 0.943778538812785 1;0 0.955836187214611 1;0 0.967893835616438 1;0 0.979951484018264 1;0 0.992009132420091 1;0.00406678082191725 1 0.995933219178083;0.0161244292237437 1 0.983875570776256;0.0281820776255701 1 0.97181792237443;0.0402397260273966 1 0.959760273972603;0.052297374429223 1 0.947702625570777;0.0643550228310494 1 0.935644977168951;0.0764126712328759 1 0.923587328767124;0.0884703196347023 1 0.911529680365298;0.100527968036529 1 0.899472031963471;0.112585616438355 1 0.887414383561645;0.124643264840182 1 0.875356735159818;0.136700913242008 1 0.863299086757992;0.148758561643834 1 0.851241438356166;0.160816210045661 1 0.839183789954339;0.172873858447487 1 0.827126141552513;0.184931506849314 1 0.815068493150686;0.19698915525114 1 0.80301084474886;0.209046803652967 1 0.790953196347033;0.221104452054793 1 0.778895547945207;0.23316210045662 1 0.76683789954338;0.245219748858446 1 0.754780251141554;0.257277397260272 1 0.742722602739728;0.269335045662099 1 0.730664954337901;0.281392694063925 1 0.718607305936075;0.293450342465752 1 0.706549657534248;0.305507990867578 1 0.694492009132422;0.317565639269405 1 0.682434360730595;0.329623287671231 1 0.670376712328769;0.341680936073057 1 0.658319063926943;0.353738584474884 1 0.646261415525116;0.36579623287671 1 0.63420376712329;0.377853881278537 1 0.622146118721463;0.389911529680363 1 0.610088470319637;0.40196917808219 1 0.59803082191781;0.414026826484016 1 0.585973173515984;0.426084474885843 1 0.573915525114157;0.438142123287669 1 0.561857876712331;0.450199771689495 1 0.549800228310505;0.462257420091322 1 0.537742579908678;0.474315068493148 1 0.525684931506852;0.486372716894975 1 0.513627283105025;0.498430365296801 1 0.501569634703199;0.510488013698628 1 0.489511986301372;0.522545662100454 1 0.477454337899546;0.53460331050228 1 0.46539668949772;0.546660958904107 1 0.453339041095893;0.558718607305933 1 0.441281392694067;0.57077625570776 1 0.42922374429224;0.582833904109586 1 0.417166095890414;0.594891552511413 1 0.405108447488587;0.606949200913239 1 0.393050799086761;0.619006849315066 1 0.380993150684934;0.631064497716892 1 0.368935502283108;0.643122146118718 1 0.356877853881282;0.655179794520545 1 0.344820205479455;0.667237442922371 1 0.332762557077629;0.679295091324198 1 0.320704908675802;0.691352739726024 1 0.308647260273976;0.703410388127851 1 0.296589611872149;0.715468036529677 1 0.284531963470323;0.727525684931503 1 0.272474315068497;0.73958333333333 1 0.26041666666667;0.751640981735156 1 0.248359018264844;0.763698630136983 1 0.236301369863017;0.775756278538809 1 0.224243721461191;0.787813926940636 1 0.212186073059364;0.799871575342462 1 0.200128424657538;0.811929223744289 1 0.188070776255711;0.823986872146115 1 0.176013127853885;0.836044520547941 1 0.163955479452059;0.848102168949768 1 0.151897831050232;0.860159817351594 1 0.139840182648406;0.872217465753421 1 0.127782534246579;0.884275114155247 1 0.115724885844753;0.896332762557074 1 0.103667237442926;0.9083904109589 1 0.0916095890410999;0.920448059360726 1 0.0795519406392735;0.932505707762553 1 0.0674942922374471;0.944563356164379 1 0.0554366438356206;0.956621004566206 1 0.0433789954337942;0.968678652968032 1 0.0313213470319678;0.980736301369859 1 0.0192636986301413;0.992793949771685 1 0.00720605022831489;1 0.995148401826488 0;1 0.983090753424662 0;1 0.971033105022836 0;1 0.958975456621009 0;1 0.946917808219183 0;1 0.934860159817356 0;1 0.92280251141553 0;1 0.910744863013703 0;1 0.898687214611877 0;1 0.886629566210051 0;1 0.874571917808224 0;1 0.862514269406398 0;1 0.850456621004571 0;1 0.838398972602745 0;1 0.826341324200918 0;1 0.814283675799092 0;1 0.802226027397265 0;1 0.790168378995439 0;1 0.778110730593613 0;1 0.766053082191786 0;1 0.75399543378996 0;1 0.741937785388133 0;1 0.729880136986307 0;1 0.71782248858448 0;1 0.705764840182654 0;1 0.693707191780828 0;1 0.681649543379001 0;1 0.669591894977175 0;1 0.657534246575348 0;1 0.645476598173522 0;1 0.633418949771695 0;1 0.621361301369869 0;1 0.609303652968042 0;1 0.597246004566216 0;1 0.58518835616439 0;1 0.573130707762563 0;1 0.561073059360737 0;1 0.54901541095891 0;1 0.536957762557084 0;1 0.524900114155257 0;1 0.512842465753431 0;1 0.500784817351605 0;1 0.488727168949778 0;1 0.476669520547952 0;1 0.464611872146125 0;1 0.452554223744299 0;1 0.440496575342472 0;1 0.428438926940646 0;1 0.416381278538819 0;1 0.404323630136993 0;1 0.392265981735167 0;1 0.38020833333334 0;1 0.368150684931514 0;1 0.356093036529687 0;1 0.344035388127861 0;1 0.331977739726034 0;1 0.319920091324208 0;1 0.307862442922382 0;1 0.295804794520555 0;1 0.283747146118729 0;1 0.271689497716902 0;1 0.259631849315076 0;1 0.247574200913249 0;1 0.235516552511423 0;1 0.223458904109596 0;1 0.21140125570777 0;1 0.199343607305944 0;1 0.187285958904117 0;1 0.175228310502291 0;1 0.163170662100464 0;1 0.151113013698638 0;1 0.139055365296811 0;1 0.126997716894985 0;1 0.114940068493159 0;1 0.102882420091332 0;1 0.0908247716895056 0;1 0.0787671232876792 0;1 0.0667094748858528 0;1 0.0546518264840263 0;1 0.0425941780821999 0;1 0.0305365296803735 0;1 0.018478881278547 0;1 0.00642123287672058 0;0.994363584474894 0 0;0.982305936073068 0 0;0.970248287671241 0 0;0.958190639269415 0 0;0.946132990867588 0 0;0.934075342465762 0 0;0.922017694063936 0 0;0.909960045662109 0 0;0.897902397260283 0 0;0.885844748858456 0 0;0.87378710045663 0 0;0.861729452054803 0 0;0.849671803652977 0 0;0.83761415525115 0 0;0.825556506849324 0 0;0.813498858447498 0 0;0.801441210045671 0 0;0.789383561643845 0 0;0.777325913242018 0 0;0.765268264840192 0 0;0.753210616438365 0 0;0.741152968036539 0 0;0.729095319634713 0 0;0.717037671232886 0 0;0.70498002283106 0 0;0.692922374429233 0 0;0.680864726027407 0 0;0.66880707762558 0 0;0.656749429223754 0 0;0.644691780821927 0 0;0.632634132420101 0 0;0.620576484018275 0 0;0.608518835616448 0 0;0.596461187214622 0 0;0.584403538812795 0 0;0.572345890410969 0 0;0.560288242009142 0 0;0.548230593607316 0 0;0.53617294520549 0 0;0.524115296803663 0 0;0.512057648401837 0 0;0.50000000000001 0 0],...
        'FontSize',fontSize,'Ticks',labelTicks,'Limits',[0,round(data_max,1)]); % 1 hour
    caxis([round(data_min,1) round(data_max,1)]);
    %     colormap jet;
    ylabel(h, 'NRMSE of SOC estimate ($\%$)','Interpreter','latex','FontSize',fontSize);
    fig_pos = figure_1.PaperPosition;
    figure_1.PaperSize = [fig_pos(3) fig_pos(4)];
    if(savefiles)
        saveas(gcf,strcat(filename,'.fig'));
        saveas(gcf,strcat(filename,'.eps'),'epsc');
    end
end

end