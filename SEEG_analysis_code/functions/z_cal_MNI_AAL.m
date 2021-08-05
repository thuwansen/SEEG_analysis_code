function [MNI_coor,aal_idx] = z_cal_MNI_AAL(savepath,lead_recon_path,data_MO_label,MNI_atlas_nii,raal_atlas_nii,CLA_name)
% This function calculated the coordinates from the results of Lead-DBS
%INPUTS
% - savepath          : savepath of the reuslts SEEG_COORDINATE.mat
% - lead_recon_path   : path of the lead reconstructed in Lead-DBS
% - data_MO_label     : name for each of the electrode
% - MNI_atlas_nii     : the path of MNI atlas
% - raal_atlas_nii    : the path of the AAL atlas
% - CLA_name          : name of the electrode contact in claustrum
%OUTPUTS
% - MNI_coor          : coordinate for each of the electrode contact in MNI space 
% - spindle_de        : the start time point,end time point and time window of the detected spindles
% - aal_idx           : the AAL brain area for each of the electrode contact


leads = dir([lead_recon_path,'\ea_reconstruction_*']);
mni_all_lead = [];
native_all_lead=[];
name_all_lead=[];
for j = 1:length(leads)
    disp(j);
    tmp = load([lead_recon_path,'\',leads(j).name]);
    mni_coor = tmp.reco.mni.coords_mm{1};
    mni_all_lead = [mni_all_lead;mni_coor];
    native_coor = tmp.reco.native.coords_mm{1};
    native_all_lead = [native_all_lead;native_coor];  
    tmp_name = strsplit(leads(j).name,{'_','.'});
    for k=1:length(mni_coor)
        temp_lead_name = cellstr([tmp_name{4},num2str(k)]);
        name_all_lead = [name_all_lead;temp_lead_name];
    end

end
lead_name = data_MO_label;
mni_coor_final = zeros(length(lead_name),3);
native_coor_final = zeros(length(lead_name),3);
for i=1:length(lead_name)
    idx = find(ismember(name_all_lead,lead_name{i}));
    if ~isempty(idx)
        mni_coor_final(i,:) = mni_all_lead(idx,:);
    end
end

% calculate the AAL brain area index for each of the electrode contact
aal_hdr = spm_vol(raal_atlas_nii);
aal_img = spm_read_vols(aal_hdr);
hdr = spm_vol(MNI_atlas_nii);
img = spm_read_vols(hdr);
trans_M = hdr.mat;

[lead_num,~]=size(mni_coor_final);
aal_idx = zeros(lead_num,1);
for lead_idx=1:lead_num
    lead_point = [mni_coor_final(lead_idx,:),1]';
    coor = round(trans_M\lead_point);
    try
        aal_idx(lead_idx,1)=aal_img(coor(1),coor(2),coor(3));
    catch
    end
end

for i=1:length(CLA_name)
    cla_channel =CLA_name{i};
    idx_cla = find(ismember(data_MO_label,cla_channel));
    if ismember('''',cla_channel)
        aal_idx(idx_cla,1)=117;
    else
        aal_idx(idx_cla,1)=118;
    end
end
lead_coor_info = {mni_coor_final,aal_idx,native_coor_final,lead_name}; 
MNI_coordinate = mni_coor_final;
AAL_index = aal_idx;
SEEG_leadlabel = lead_name;
SEEG_COORDINATE = struct('MNI_coordinate',{MNI_coordinate},...
                          'AAL_index',{AAL_index},...
                          'SEEG_leadlabel',{SEEG_leadlabel});
save([savepath,'/SEEG_COORDINATE.mat'],'SEEG_COORDINATE');
MNI_coor = mni_coor_final;



