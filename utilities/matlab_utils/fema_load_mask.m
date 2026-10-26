function mask = fema_load_mask(dirname_imaging, datatype)
%FEMA_LOAD_MASK Load analysis mask the same way FEMA_process_data does (voxel).
%
% Voxel: vol_mask_sub from vol_info.mat in dirname_imaging (FEMA_process_data.m).

datatype = lower(strtrim(char(datatype)));

switch datatype
    case {'voxel', 'voxelwise'}
        fname_volinfo = fullfile(dirname_imaging, 'vol_info.mat');
        if ~isfile(fname_volinfo)
            fname_volinfo = fullfile(dirname_imaging, 'volinfo.mat');
        end
        if ~isfile(fname_volinfo)
            error('fema_load_mask: vol_info.mat not found under %s', dirname_imaging);
        end
        tmp = load(fname_volinfo, 'vol_mask_sub');
        mask = tmp.vol_mask_sub;
    case {'vertex', 'vertexwise'}
        error(['fema_load_mask: vertex mask is not stored in vol_info.mat. ', ...
               'Re-run with mask saved in legacy results or extend this helper.']);
    otherwise
        error('fema_load_mask: unsupported datatype ''%s''.', datatype);
end
end
