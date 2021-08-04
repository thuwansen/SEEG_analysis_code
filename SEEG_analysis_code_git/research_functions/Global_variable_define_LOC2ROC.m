% Global variables define for the current study
% By Sen Wan.

subject_info_path = 'D:\SEEG_preoject\BIDS_data\bids\metadata.csv';
dataset_path = 'D:\SEEG_preoject\BIDS_data\bids\';
function_path = 'D:\SEEG_preoject\BIDS_data\code\research_functions';

subs_name = {
     'sub-0031';
     'sub-0032';
     'sub-0033';
     'sub-0034';
     'sub-0035';
     'sub-0036';
     'sub-0037';
     'sub-0038';
     'sub-0039';};

data_MO_name_ana = {
     [dataset_path,'sub-0031\ieeg\SEEG.mat'];
     [dataset_path,'sub-0032\ieeg\SEEG.mat'];
     [dataset_path,'sub-0033\ieeg\SEEG.mat'];
     [dataset_path,'sub-0034\ieeg\SEEG.mat'];
     [dataset_path,'sub-0035\ieeg\SEEG.mat'];
     [dataset_path,'sub-0036\ieeg\SEEG.mat'];
     [dataset_path,'sub-0037\ieeg\SEEG.mat'];
     [dataset_path,'sub-0038\ieeg\SEEG.mat'];
     [dataset_path,'sub-0039\ieeg\SEEG.mat'];};

SEEG_coordinate_path = {
     [dataset_path,'sub-0031\ieeg\SEEG_COORDINATE.mat'];
     [dataset_path,'sub-0032\ieeg\SEEG_COORDINATE.mat'];
     [dataset_path,'sub-0033\ieeg\SEEG_COORDINATE.mat'];
     [dataset_path,'sub-0034\ieeg\SEEG_COORDINATE.mat'];
     [dataset_path,'sub-0035\ieeg\SEEG_COORDINATE.mat'];
     [dataset_path,'sub-0036\ieeg\SEEG_COORDINATE.mat'];
     [dataset_path,'sub-0037\ieeg\SEEG_COORDINATE.mat'];
     [dataset_path,'sub-0038\ieeg\SEEG_COORDINATE.mat'];
     [dataset_path,'sub-0039\ieeg\SEEG_COORDINATE.mat'];};

cla_selected = {
    'FI''5';
    'CI1';
    'CI2';
    'FI4';
    'PI3';
    'FI''4';
    'TI2';
    'FI''1';
    'FI''4';};

ob_selected = {
    'OB''8';
    'OB8';
    '';
    'OB4';  
    'OB7';
    'OB''12';
    'OB11';
    '';
    'OB''7';};

fi_selected = {
    'FI''13';
    'FI10';
    'FI11';
    'FI6';
    'FI11';
    'FI''14';
    'FI10';
    'FI''11'; 
    'FI''11';};

ana_time = {
    300,400;
    200,300;
    560,660;
    1000,1100;
    350,400;
    300,400;
    200,300;
    100,200;
    360,460;};

con_time = {
    770,870;
    686,786;
    1535,1635;
    3000,3100;
    1757,1807;
    800,900;
    1200,1300;
    900,1000;
    1730,1830;};





