clear 
close all
addpath(genpath('.'));
datasets ={ 'Dino 319' 
 'Dino 4983' 
 'Corridor' 
 'House' 
 'Gustav Vasa' 
 'Folke Filbyter' 
 'Park Gate' 
 'Nijo' 
'Drinking Fountain' 
 'Golden Statue' 
'Jonas Ahls' 
 'De Guerre' 
 'Dome' 
 'Alcatraz Courtyard' 
'Alcatraz Water Tower' 
'Cherub'
'Pumpkin'
'Sphinx'
 'Toronto University'  
 'Sri Thendayuthapani' 
'Porta san Donato'
 'Buddah Tooth' 
'Tsar Nikolai I'
'Smolny Cathedral'
'Skansen Kronan'};


result=struct;
for i=1:length(datasets)
    [ repError,elapsedTime,sizes ] = runProjective( datasets{i},1 );
    result(i).name= datasets{i};
    result(i).repError= repError;
    result(i).elapsedTime= elapsedTime;
    
    result(i).sizes= sizes;
end


%%
t=table({result.name}',[result.repError]',[result.elapsedTime]', reshape([result.sizes],2,[])', 'VariableNames',{'Dataset', 'Reprojection_Error' ,'Runtime', 'cameras_points'});
writetable(t,'table.csv');