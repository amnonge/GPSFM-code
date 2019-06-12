clear
addpath(genpath('.'));
% mexopencv.make

%Solved Data Sets
% dataset = 'Dino 319' 
% dataset = 'Dino 4983' 
% dataset = 'Corridor' 
% dataset = 'House' 
% dataset = 'Gustav Vasa' 
% dataset = 'Folke Filbyter' 
% dataset = 'Park Gate' 
% dataset = 'Nijo' 
% dataset = 'Drinking Fountain' 
dataset = 'Golden Statue' 
% dataset = 'Jonas Ahls' 
% dataset = 'De Guerre' 
% dataset = 'Dome' 
% dataset = 'Alcatraz Courtyard' 
% dataset = 'Alcatraz Water Tower' 
% dataset = 'Cherub';
% dataset = 'Pumpkin'
% dataset = 'Sphinx'; 
% dataset = 'Toronto University'  
% dataset = 'Sri Thendayuthapani' 
% dataset = 'Porta san Donato'
% dataset = 'Buddah Tooth' 
% dataset = 'Tsar Nikolai I'
% dataset = 'Smolny Cathedral'
% dataset = 'Skansen Kronan'



runProjective(dataset,1);
