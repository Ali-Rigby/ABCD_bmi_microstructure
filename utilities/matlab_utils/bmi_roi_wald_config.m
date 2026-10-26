function cfg = bmi_roi_wald_config()
%BMI_ROI_WALD_CONFIG ROI definitions for mean-Wald ranking (s9b).
%
% Groups:
%   bilateral  - aseg left+right combined (same pairs as s9 bilateral output)
%   fiber      - forceps + combined CST (L_CST | R_CST)
%   pauli      - Pauli subcortical nuclei (non-lateralized)

cfg.dirname_atlas = '/space/ceph/1/ABCD/users/parekh/2025-11-19_ABCD_concat_uncompressed/6_0/concatenated/imaging/voxelwise/roi';
cfg.prob_thresh = 0.7;

% Bilateral aseg structures (base name without Left-/Right- prefix)
cfg.bilateral_aseg = {
    'Thalamus-Proper'
    'Caudate'
    'Putamen'
    'Pallidum'
    'Accumbens-area'
    'Amygdala'
    'VentralDC'
    'Cerebral-White-Matter'
    'Cerebral-Cortex'
    'Hippocampus'
};

% Fiber tracts (fiber.mat) — keep in sync with s9 roi_mod fiber list
cfg.fiber_single = {'Fmaj', 'Fmin', 'CC'};
cfg.fiber_bilateral = {
    'Fx', 'CgC', 'CgH', 'CST', 'ATR', 'Unc', 'ILF', 'IFO', ...
    'SLF', 'tSLF', 'pSLF', 'SCS', 'fSCS', 'pSCS', 'SIFC', 'IFSFC'};

% Pauli nuclei (pauli_subcortnuclei_abcd3_cor10.mat)
cfg.pauli = {
    'Putamen'
    'Caudate'
    'Nucleus Accumbens'
    'Extended Amygdala'
    'Globus Pallidus external'
    'Globus Pallidus internal'
    'Substantia Nigra pars compacta'
    'Red Nucleus'
    'Substantia Nigra pars reticulata'
    'Parabrachial Pigmented Nucleus'
    'Ventral Tegmental Area'
    'Ventral Pallidum'
    'Habenular Nuclei'
    'Hypothalamus'
    'Mammillary Nucleus'
    'Subthalamic Nucleus'
};

end
