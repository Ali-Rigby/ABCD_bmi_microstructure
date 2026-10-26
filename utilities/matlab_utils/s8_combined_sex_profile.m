function tf = s8_combined_sex_profile(analysis_profile)
%S8_COMBINED_SEX_PROFILE True for profiles with designmat1_* (not 1f/1m).
tf = ~isempty(analysis_profile) && strncmp(analysis_profile, 've11_nove11_', 12);
end
