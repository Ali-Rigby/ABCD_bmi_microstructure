function abcd_sync_path = fema_default_abcd_sync_path()
%FEMA_DEFAULT_ABCD_SYNC_PATH ABCD 6.0 concat root (same as s6_run_FEMA.m).
abcd_sync_path = getenv('BMI_IMAGING_ABCD_SYNC');
if isempty(abcd_sync_path)
    abcd_sync_path = fullfile('/space/ceph/1/ABCD/users/parekh/2025-11-19_ABCD_concat_uncompressed/6_0');
end
end
