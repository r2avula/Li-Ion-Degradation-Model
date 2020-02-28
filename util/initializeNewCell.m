function [cellSimParams] = initializeNewCell(cellSimParams)
model = cellSimParams.model;
init_soc = cellSimParams.SOC_init;

model.sol('sol1').clearSolutionData;
model.sol('sol2').clearSolutionData;

model.param.set('period',cellSimParams.slotIntervalInSeconds);
model.param.set('dt',cellSimParams.slotIntervalInSeconds/100);
model.param.set('SOC_init',init_soc);
model.param.set('P_ref',0);
model.param.set('SOC_high',cellSimParams.SOC_high);
model.param.set('SOC_low',cellSimParams.SOC_low);
model.param.set('des_SOC',0.9);
model.param.set('E_max',cellSimParams.cell_voltage_high);
model.param.set('E_min',cellSimParams.cell_voltage_low);
model.sol('sol2').feature('t1').feature('st1').set('storestopcondsol', 'stepbefore');
model.sol('sol2').feature('t1').feature('st1').setIndex('stopcondActive', true, 0); % SOC out of range
model.sol('sol2').feature('t1').feature('st1').setIndex('stopcondActive', true, 1); % Ecell out of range
model.sol('sol2').feature('t1').feature('st1').setIndex('stopcondActive', false, 2); % charging
model.sol('sol2').feature('t1').feature('st1').setIndex('stopcondActive', false, 3); % discharging

model.study('std1').run;% cell initialization with init_soc

model.study('std2').feature('time1').set('initstudy', 'std1'); % for COMSOL v5.4

model.study('std2').run;% keep idle for one time slot

model.study('std2').feature('time1').set('initstudy', 'std2'); % for COMSOL v5.4

cellSimParams.cur_voltage = mphglobal(model,'Ecell','dataset','dset2','solnum','end');
cellSimParams.model = model;
cellSimParams.cur_soc = mphglobal(model,'SOCcell_load','dataset','dset2','solnum','end');
cellSimParams.cellsIterated = cellSimParams.cellsIterated + 1;
cellSimParams.currentCellRelativeCapacity = 100;


desCellRelativeCapacity = cellSimParams.initialRelCap;
[cellSimParams] = driveToRelCap(cellSimParams,desCellRelativeCapacity);
end