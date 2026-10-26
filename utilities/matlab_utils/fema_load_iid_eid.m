function [iid, eid, src_file] = fema_load_iid_eid(dirname_out, datatype, fstem_imaging)
%FEMA_LOAD_IID_EID Load analytic visit IDs from FEMA 4.x / legacy outputs.
%
% FEMA 4.x saves iid/eid in *_designMatrix.mat (not *_estimates.mat).

src_file = '';
iid = [];
eid = [];

% 1) Design-matrix sidecar (FEMA 4.x default)
try
    fpath = fema_design_matrix_path(dirname_out, datatype, fstem_imaging);
    [iid, eid] = load_iid_eid_from_mat(fpath);
    src_file = fpath;
    return;
catch
    % fall through
end

% 2) Legacy combined estimates file (pre-split save)
try
    fpath = fema_results_path(dirname_out, datatype, fstem_imaging);
    est_names = {whos('-file', fpath).name};
    if all(ismember({'iid', 'eid'}, est_names))
        S = load(fpath, 'iid', 'eid');
        iid = S.iid(:);
        eid = S.eid(:);
        src_file = fpath;
        return;
    end
catch
    % fall through
end

error('fema_load_iid_eid: could not find iid/eid under %s (%s / %s).', ...
    dirname_out, datatype, fstem_imaging);
end

function [iid, eid] = load_iid_eid_from_mat(fpath)
vars = {whos('-file', fpath).name};

if all(ismember({'iid', 'eid'}, vars))
    S = load(fpath, 'iid', 'eid');
    iid = S.iid(:);
    eid = S.eid(:);
    return;
end

if ismember('designMatrix', vars)
    S = load(fpath, 'designMatrix');
    dm = S.designMatrix;
    if istable(dm)
        vn = dm.Properties.VariableNames;
        if all(ismember({'participant_id', 'session_id'}, vn))
            iid = dm.participant_id;
            eid = dm.session_id;
            return;
        end
        if all(ismember({'iid', 'eid'}, vn))
            iid = dm.iid;
            eid = dm.eid;
            return;
        end
    end
end

error('fema_load_iid_eid:%s has no iid/eid or designMatrix IDs.', fpath);
end
