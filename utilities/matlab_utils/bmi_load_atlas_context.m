function ctx = bmi_load_atlas_context(dirname_atlas, atlas_type)
%BMI_LOAD_ATLAS_CONTEXT Load atlas labels and probability volumes for ROI masks.

fname_atlas = fullfile(dirname_atlas, [atlas_type '.mat']);
if ~isfile(fname_atlas)
    error('Atlas not found: %s', fname_atlas);
end

ctx = struct('atlas_type', atlas_type, 'fname_atlas', fname_atlas);

switch atlas_type
    case 'thalamus_abcd3_cor10'
        atlas = load(fname_atlas, 'thalamus_lut', 'vol_labels_reg');
        ctx.roinames = atlas.thalamus_lut.LabelName;
        ctx.volmat_prob = atlas.vol_labels_reg;
    case 'pauli_subcortnuclei_abcd3_cor10'
        atlas = load(fname_atlas, 'pauli_lut', 'vol_labels_reg');
        ctx.roinames = atlas.pauli_lut.Name;
        ctx.volmat_prob = atlas.vol_labels_reg;
    otherwise
        atlas = load(fname_atlas, 'roinames', 'volmat_prob');
        ctx.roinames = atlas.roinames;
        ctx.volmat_prob = atlas.volmat_prob;
end

end
