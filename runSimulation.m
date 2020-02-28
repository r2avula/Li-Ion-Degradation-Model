clear;
simStartup();
comsol_model_filename = [pwd filesep 'util' filesep 'li_degradation_core.mph'];

config_filename = 'degradation_20percent_1hour.yaml';

config = ReadYaml(config_filename);
config.comsol_model_filename = comsol_model_filename;

batteryRatedCapacityInAh = 50;
config.batteryRatedCapacityInAh = batteryRatedCapacityInAh;
cellSimData_out = getDegradationData(config);
[mean_nrmse_energyLoss_percentage,mean_nrmse_soc_kp1_percentage] = plotMeanDegradationMaps_1hour(config,cellSimData_out);

