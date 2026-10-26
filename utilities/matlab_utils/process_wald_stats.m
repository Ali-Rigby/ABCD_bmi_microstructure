function vols = process_wald_stats(vols, dirname_results, modality, xname, designmat_label, pwaldlimits, pwaldlimits_Bon, pwaldlimits_FDR, showPval, cmap, bg_vol, CSFprob, interpMethod, atlasVersion, varargin)
%PROCESS_WALD_STATS Append Wald -log10(p) maps to a showVol volume cell array.
%
% Optional name-value:
%   'Corrections' - cell of correction types to include (default: all):
%       'uncorr', 'bonferroni', 'fdr'
%
% Example (FDR only):
%   vols = process_wald_stats(vols, ..., 'Corrections', {'fdr'});

    p = inputParser;
    addParameter(p, 'Corrections', {'uncorr', 'bonferroni', 'fdr'});
    parse(p, varargin{:});
    corrections = lower(string(p.Results.Corrections));

    splines_file = sprintf('%s/%s_Wald_splines_CI_output_%s.mat', dirname_results, modality, xname);
    if ~isfile(splines_file)
        warning('Splines output file not found for %s', xname);
        return;
    end
    S = load(splines_file);

    if any(corrections == "uncorr")
        stat_labels_uncorrected = {'vol_logp_l', 'vol_logp_s', 'vol_logp_ls'};
        for s = 1:numel(stat_labels_uncorrected)
            fld = stat_labels_uncorrected{s};
            if isfield(S, fld)
                statname = sprintf('%s %s, %s', xname, strrep(fld, 'vol_', ''), designmat_label);
                vols = convert_and_append_vols(vols, S.(fld), modality, statname, {xname}, 1, pwaldlimits, showPval, cmap, bg_vol, CSFprob, interpMethod, atlasVersion);
            end
        end
    end

    if any(corrections == "bonferroni")
        stat_labels_bonferroni = {'vol_logp_l_Bon', 'vol_logp_s_Bon', 'vol_logp_ls_Bon'};
        for s = 1:numel(stat_labels_bonferroni)
            fld = stat_labels_bonferroni{s};
            if isfield(S, fld)
                statname = sprintf('%s %s, %s', xname, strrep(fld, 'vol_', ''), designmat_label);
                vols = convert_and_append_vols(vols, S.(fld), modality, statname, {xname}, 1, pwaldlimits_Bon, showPval, cmap, bg_vol, CSFprob, interpMethod, atlasVersion);
            end
        end
    end

    if any(corrections == "fdr")
        stat_labels_fdr = {'vol_logp_l_FDR', 'vol_logp_s_FDR', 'vol_logp_ls_FDR'};
        for s = 1:numel(stat_labels_fdr)
            fld = stat_labels_fdr{s};
            if isfield(S, fld)
                statname = sprintf('%s %s, %s', xname, strrep(fld, 'vol_', ''), designmat_label);
                vols = convert_and_append_vols(vols, S.(fld), modality, statname, {xname}, 1, pwaldlimits_FDR, showPval, cmap, bg_vol, CSFprob, interpMethod, atlasVersion);
            end
        end
    end
end
