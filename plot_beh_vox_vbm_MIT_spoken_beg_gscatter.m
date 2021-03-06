%% Grey Matter Volume in Baseline Chronic Stroke with Aphasia: data visualization for MIT of phrases in spoken format 
% (C) Noelia Martinez-Molina, MIT License
clear all
% Specify directories and patients
data_path='G:\Aphasia_project\VBM_v4\data_v4';
covbeh_path='G:\Aphasia_project\VBM_v4\\analysis\covariates\behaviour\MIT\';
peak_path='G:\Aphasia_project\VBM_v4\analysis\L2\MIT\GMV_WPM_MIT_Spoken_Total_Dur_Beg_TIV_Age_TP1\';
output_path='G:\Aphasia_project\VBM_v4\quality_checks\MIT\';
group_path='G:\Aphasia_project\Behaviour\baseline\';
names= dir(data_path);
names(ismember({names.name},{'.','..'}))=[];
ses='ses-001';
anat='anat';
prep_folder='spm_us_cfm_maskingoutlesion_cleanup_NEW_TPM_med_reg_old_NORM'; 
% Create variable with smwrp1 images
n=1;
for sub=1:size(names,1)
    % Exclude patients with no lesions: sub-24(ID143); sub-31 (ID154); sub-32(ID155); sub-33(ID157); sub-35(ID159)
    if ~strcmp(names(sub).name, 'sub-24')  && ~strcmp(names(sub).name, 'sub-31') && ~strcmp(names(sub).name, 'sub-32')  && ~strcmp(names(sub).name, 'sub-33') && ~strcmp(names(sub).name, 'sub-35')
%         display(sub)
        sub_path=fullfile(data_path, names(sub).name, ses, anat, prep_folder);
        smwrc1{n,1}=spm_select('List', fullfile(sub_path), '^smwc1msk_sub.*\.nii$'); %grey matter
        smwrc1_names{n,1}=fullfile(sub_path,smwrc1{n,1});
        n=n+1;
    end
end
%% Get peak voxel intensity for all patients
%%
% Save peak Voxel in mm and loop for the clusters from different L2 analysis
        cd(peak_path)
        peaks_xyz=xlsread('WPM_MIT_Spoken_Total_Dur_Beg_TIV_Age_p001unc_k50', 2, 'I2:K6'); %Whole-brain L2 results
        % Loop for patients
        for c=1:size(peaks_xyz,1)
            for i=1:size(smwrc1,1)
%% Get image information
%%
                V=spm_vol(smwrc1_names{i,1});
                [Y, XYZmm]=spm_read_vols(V);
%% Get the volume index
%%
                idx=find(XYZmm(1,:)==peaks_xyz(c,1) & XYZmm(2,:)==peaks_xyz(c,2) & XYZmm(3,:)==peaks_xyz(c,3));
                XYZmm(:,idx);
%% Get the peak coordinates in Voxels
%%
                [Xvox, Yvox, Zvox]=ind2sub(size(Y(:,:,:,1)),idx);
%% Get the intensity value for the peak
%%
                int=Y(Xvox,Yvox,Zvox,:);
%% Save intensity value for the peak from all patients
%%
                pa_int{i,1}=names(i,1).name;
                pa_int{i,c+1}=int;
            end
        end
%% Get behavioural value for all patients
%%
WPM_MIT_Spoken_Total_Dur_Beg_TP1_N45= load(fullfile(covbeh_path,'Beh_MIT_WPM_Spoken_Total_Dur_Beg_TP1_N45.mat'));
Group=xlsread(fullfile(group_path,'LASA_aphasia_group.xlsx'));
% Remove patients with NaN values in behaviour
if sum(isnan(WPM_MIT_Spoken_Total_Dur_Beg_TP1_N45.beh_cov_MIT(:,2)))~=0
    NAN_Index=find(isnan(WPM_MIT_Spoken_Total_Dur_Beg_TP1_N45.beh_cov_MIT(:,2)));
    WPM_MIT_Spoken_Total_Dur_Beg_TP1_N45.beh_cov_MIT(NAN_Index,:)=[];
    Group(NAN_Index,:)=[];
    for j=1:size(NAN_Index,1)
        pa_int(NAN_Index(j,1),:)=[];
    end
    pa_int_fixed=pa_int;
end
%% Plot Intensity against behaviour for the peak coordinate
%%
for p=1:size(peaks_xyz,1)
    if exist('pa_int_fixed')
        x=WPM_MIT_Spoken_Total_Dur_Beg_TP1_N45.beh_cov_MIT(:,2);
        y=cell2mat(pa_int_fixed(:,p+1));
        g=Group(:,5);
        labels=pa_int_fixed(:,1);
        figure
        gscatter(x,y,g)
        title ('VBM: GMV vs MIT SPOKEN'), ylabel(['GM Intensity Peak Voxel Cluster ' num2str(p)]), xlabel ('WPM MIT SPOKEN')
        cd(output_path)
        saveas(gcf,['gscatter_VBM_WPM_MIT_Spoken_Total_Dur_Beg_Cluster_' num2str(p) '.jpg'])
        % label data points
        figure
        gscatter(x,y,g)
        title ('VBM: GMV vs MIT SPOKEN WITH PATIENTS ID'), ylabel(['GM Intensity Peak Voxel Cluster ' num2str(p)]), xlabel ('WPM MIT SPOKEN')
        for k=1:length(labels)
            text (x(k,1), y(k,1), labels {k,1}(5:6),'fontSize',8,'VerticalAlignment','top', 'HorizontalAlignment','left');
        end
        saveas(gcf,['gscatter_VBM_WPM_MIT_Spoken_Total_Dur_Beg_Cluster_' num2str(p) '_labels.jpg'])
    else
%% Plot Intensity against behaviour for the peak coordinate
%%
        x=WPM_MIT_Spoken_Total_Dur_Beg_TP1_N45.beh_cov_MIT(:,2);
        y=cell2mat(pa_int(:,p+1));
        g=Group(:,5);
        labels=pa_int(:,1);
        figure
        gscatter(x,y,g)
        title ('VBM: GMV vs MIT SPOKEN'), ylabel(['GM Intensity Peak Voxel Cluster ' num2str(p)]), xlabel ('WPM MIT Spoken')
        cd(output_path)
        saveas(gcf,['gscatter_VBM_WPM_MIT_Spoken_Total_Dur_Beg_Cluster_' num2str(p) '.jpg'])
        % label data points
        figure
        gscatter(x,y,g)
        title ('VBM: GMV vs MIT SPOKEN WITH PATIENTS ID'), ylabel(['GM Intensity Peak Voxel Cluster ' num2str(p)]), xlabel ('WPM MIT Spoken')
        for k=1:length(labels)
            text (x(k,1), y(k,1), labels {k,1}(5:6),'fontSize',8,'VerticalAlignment','top', 'HorizontalAlignment','left');
        end
        saveas(gcf,['gscatter_VBM_WPM_MIT_Spoken_Total_Dur_Beg_Cluster_' num2str(p) '_labels.jpg'])
    end
end