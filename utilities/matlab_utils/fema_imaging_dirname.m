function dirname_imaging = fema_imaging_dirname(datatype, fstem_imaging, abcd_sync_path)
%FEMA_IMAGING_DIRNAME Imaging data directory passed to FEMA_wrapper (same as s6).
%
% modality / fstem_imaging examples: rni, b0, rnt_gm_sm1000, area_sm1000, ...

if nargin < 3 || isempty(abcd_sync_path)
    abcd_sync_path = fema_default_abcd_sync_path();
end

datatype = lower(strtrim(char(datatype)));
fstem_imaging = char(fstem_imaging);

switch datatype
    case {'voxel', 'voxelwise'}
        switch lower(fstem_imaging)
            case 'b0'
                dirname_imaging = fullfile(abcd_sync_path, 'concatenated/imaging/voxelwise/smri');
            otherwise
                dirname_imaging = fullfile(abcd_sync_path, 'concatenated/imaging/voxelwise/rsi');
        end
    case {'vertex', 'vertexwise'}
        switch lower(fstem_imaging)
            case {'area_sm1000', 'thk_sm1000', 'sulc_sm1000', 't1_gwc_sm1000'}
                dirname_imaging = fullfile(abcd_sync_path, 'concatenated/imaging/vertexwise/smri');
            otherwise
                dirname_imaging = fullfile(abcd_sync_path, 'concatenated/imaging/vertexwise/rsi');
        end
    otherwise
        error('fema_imaging_dirname: unsupported datatype ''%s''.', datatype);
end
end
