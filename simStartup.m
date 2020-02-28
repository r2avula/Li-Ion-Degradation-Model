function simStartup()
pathCell = regexp(path, pathsep, 'split');
test_dir = [pwd filesep 'config'];
onPath = any(strcmpi(test_dir, pathCell));

if (~onPath)        
    path(pathdef);
    addpath(genpath('util'));
    addpath(genpath('config'));    
        
    if(ispc)
        addpath(genpath('C:\Program Files\COMSOL\COMSOL54\Multiphysics\mli\'));
    else     
        addpath('/opt/comsol/5.4/mli/');
    end   
end
rng(1,'twister');
end
