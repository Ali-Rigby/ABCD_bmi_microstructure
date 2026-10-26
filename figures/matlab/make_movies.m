%% Make AVI/MP4 movies from showVol derivative screenshot PNGs.
% Unified BMI imaging pipeline (bmip / perc_median).
%
% Prerequisite: run showvol_screenshots.m (after showvol_results_for_screenshots.m)
% so PNG frames exist under:
%   <results>/screenshots_no_colorbar/<designmat>/<modality>/<orientation>/<coords>/
%   (or screenshots_with_colorbar/, depending on hideColorbar)
%
% GUI-friendly: set knobs below (no .csh required).
% MP4 uses VideoWriter profile 'MPEG-4' when that profile exists. Otherwise the
% AVI is converted with ffmpeg (odd frame sizes are padded to even dimensions).

%% User knobs
pheno_profile = 'bmip';                 % 'bmip' | 'perc_median'
results_run   = 'results_2026_06_29';  % '' = env / auto date
desmatnum     = '1f';                  % '1f' | '1m' (matches designmat tag in filenames)
summary_modality = 'rni_derivatives';  % screenshot modality folder name
orientation   = 'axial';               % 'coronal' | 'sagittal' | 'axial'
coord_str     = 'X069_Y100_Z147';      % '' = use first available coord folder
frame_rate    = 10;
frame_dup     = 3;                     % each frame written this many times
video_formats = {'avi', 'mp4'};        % any of: 'avi', 'mp4'
% Optional hard override (leave '' to auto-build from knobs above):
dirname_png   = '';

%% Paths (this file lives in figures/matlab/)
this_dir = fileparts(mfilename('fullpath'));
repo_root = fileparts(fileparts(this_dir));
run(fullfile(repo_root, 'scripts', 'matlab_add_paths.m'));

if ~isempty(pheno_profile)
    setenv('BMI_PHENO_VARIABLE', pheno_profile);
    setenv('BMI_IMAGING_ANALYSIS_DIR', pheno_profile);
end
if ~isempty(results_run)
    setenv('BMI_IMAGING_RESULTS_RUN', results_run);
end

paths = get_project_paths();
xname = paths.pheno_xname;   % 'bmi' or 'perc_median'
fprintf('make_movies: pheno=%s | xname=%s | results=%s\n', ...
    paths.pheno_variable, xname, paths.results_dir);

%% Resolve PNG directory
if isempty(dirname_png)
    shot_root = fullfile(paths.results_dir, 'screenshots_with_colorbar');
    if ~exist(shot_root, 'dir')
        error('Screenshot root not found: %s\nRun showvol_screenshots.m first.', shot_root);
    end

    designmat_dirs = dir(fullfile(shot_root, sprintf('designmat%s_*', desmatnum)));
    designmat_dirs = designmat_dirs([designmat_dirs.isdir]);
    if isempty(designmat_dirs)
        error('No designmat%s_* folders under %s', desmatnum, shot_root);
    end
    designmat_name = designmat_dirs(1).name;
    if numel(designmat_dirs) > 1
        warning('Multiple designmat%s_* folders; using %s', desmatnum, designmat_name);
    end

    mod_dir = fullfile(shot_root, designmat_name, summary_modality, orientation);
    if ~exist(mod_dir, 'dir')
        error('Missing orientation folder: %s', mod_dir);
    end

    if isempty(coord_str)
        coord_dirs = dir(fullfile(mod_dir, 'X*'));
        coord_dirs = coord_dirs([coord_dirs.isdir]);
        if isempty(coord_dirs)
            error('No coord folders under %s', mod_dir);
        end
        coord_str = coord_dirs(1).name;
        fprintf('Using first coord folder: %s\n', coord_str);
    end

    dirname_png = fullfile(mod_dir, coord_str);
end

if ~exist(dirname_png, 'dir')
    error('PNG directory not found: %s', dirname_png);
end
fprintf('PNG dir: %s\n', dirname_png);

%% Collect and order spline-derivative frames
pattern = sprintf('*_spline_deriv_%s_%s_*_*.png', desmatnum, xname);
allSpline = dir(fullfile(dirname_png, pattern));
if isempty(allSpline)
    error('No PNGs matching %s in %s', pattern, dirname_png);
end

getPheno = @(nm) parsePhenoFromFilename(nm, xname);
pheno_vals = arrayfun(@(f) getPheno(f.name), allSpline);

validMask = ~isnan(pheno_vals);
allSpline = allSpline(validMask);
pheno_vals = pheno_vals(validMask);
if isempty(allSpline)
    error('No PNGs with parseable %s values in %s', xname, dirname_png);
end

[pheno_vals, ord] = sort(pheno_vals);
allSpline = allSpline(ord);

% Drop duplicate values (keep first after sort)
[pheno_vals, ia] = unique(pheno_vals, 'stable');
allSpline = allSpline(ia);

fprintf('Frames: %d | %s range: [%g, %g]\n', ...
    numel(pheno_vals), xname, pheno_vals(1), pheno_vals(end));
disp('Order going into video:');
disp(pheno_vals(:).');

%% Build labeled frames once, then write each requested format
dirname_video = fullfile(dirname_png, 'videos');
if ~exist(dirname_video, 'dir')
    mkdir(dirname_video);
end

out_basename = sprintf('%s_%s_%s_%s_%s', ...
    xname, orientation, desmatnum, summary_modality, coord_str);

frames = cell(numel(allSpline), 1);
for k = 1:numel(allSpline)
    Ik = imread(fullfile(allSpline(k).folder, allSpline(k).name));
    img_size = size(Ik);
    text_pos = [max(1, img_size(2) - 220), 10];
    label_text = formatPhenoOverlay(pheno_vals(k), xname);
    Ik = insertText(Ik, text_pos, label_text, ...
        'FontSize', 24, 'BoxColor', 'black', 'BoxOpacity', 0, 'TextColor', 'white');
    if size(Ik, 3) == 1
        Ik = repmat(Ik, [1 1 3]);
    end
    frames{k} = Ik;
end

video_formats = unique(lower(strtrim(video_formats)), 'stable');
fname_avi = fullfile(dirname_video, [out_basename '.avi']);
fname_mp4 = fullfile(dirname_video, [out_basename '.mp4']);
want_avi = any(strcmp(video_formats, 'avi'));
want_mp4 = any(strcmp(video_formats, 'mp4'));

if want_avi || want_mp4
    writeFrameVideo(fname_avi, frames, frame_rate, frame_dup);
end

if want_mp4
    wrote_mp4 = false;
    vw = [];
    try
        vw = VideoWriter(fname_mp4, 'MPEG-4');
        vw.FrameRate = frame_rate;
        if isprop(vw, 'Quality')
            vw.Quality = 95;
        end
        open(vw);
        for k = 1:numel(frames)
            for r = 1:frame_dup
                writeVideo(vw, frames{k});
            end
        end
        close(vw);
        wrote_mp4 = true;
        fprintf('Wrote video: %s\n', fname_mp4);
    catch ME
        if ~isempty(vw) && isvalid(vw)
            close(vw);
        end
        warning('MATLAB MPEG-4 writer unavailable (%s). Converting the AVI with ffmpeg.', ME.message);
    end
    if ~wrote_mp4
        convertAviToMp4(fname_avi, fname_mp4);
    end
end

%% Helpers
function v = parsePhenoFromFilename(filename, xname)
%PARSEPHENOFROMFILENAME Extract phenotype value from screenshot fstem.
% Supports:
%   ..._bmi_22_...
%   ..._perc_median_10_...
%   ..._perc_median_m28p2_...   (-28.2)
    token = regexp(filename, ['_' regexptranslate('escape', xname) '_(m?\d+(?:p\d+)?)_'], ...
        'tokens', 'once');
    if isempty(token)
        v = NaN;
        return;
    end
    s = token{1};
    sign = 1;
    if startsWith(s, 'm')
        sign = -1;
        s = s(2:end);
    end
    s = strrep(s, 'p', '.');
    v = sign * str2double(s);
end

function writeFrameVideo(fname_video, frames, frame_rate, frame_dup)
    vw = VideoWriter(fname_video);  % Motion JPEG AVI
    vw.FrameRate = frame_rate;
    if isprop(vw, 'Quality')
        vw.Quality = 95;
    end
    open(vw);
    for k = 1:numel(frames)
        for r = 1:frame_dup
            writeVideo(vw, frames{k});
        end
    end
    close(vw);
    fprintf('Wrote video: %s\n', fname_video);
end

function convertAviToMp4(fname_avi, fname_mp4)
    ffmpeg = findFfmpeg();
    if isempty(ffmpeg)
        warning(['ffmpeg was not found, so the MP4 was not written. ', ...
            'Set BMI_FFMPEG to the ffmpeg executable and rerun.']);
        return;
    end
    cmd = sprintf(['"%s" -y -i "%s" -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2" ', ...
        '-c:v libx264 -pix_fmt yuv420p "%s"'], ffmpeg, fname_avi, fname_mp4);
    status = system(cmd);
    if status ~= 0 || ~isfile(fname_mp4)
        warning('ffmpeg failed to write %s', fname_mp4);
        return;
    end
    fprintf('Wrote video: %s\n', fname_mp4);
end

function ffmpeg = findFfmpeg()
    ffmpeg = strtrim(getenv('BMI_FFMPEG'));
    if ~isempty(ffmpeg) && isfile(ffmpeg)
        return;
    end
    candidates = { ...
        '/home/arigby/programs/anaconda3/envs/sci/bin/ffmpeg', ...
        '/home/arigby/programs/anaconda3/envs/r/bin/ffmpeg'};
    for i = 1:numel(candidates)
        if isfile(candidates{i})
            ffmpeg = candidates{i};
            return;
        end
    end
    [status, out] = system('command -v ffmpeg');
    if status == 0
        ffmpeg = strtrim(out);
        return;
    end
    ffmpeg = '';
end

function label_text = formatPhenoOverlay(v, xname)
    switch xname
        case 'bmi'
            label_text = sprintf('BMI %g%%', v);
        case 'perc_median'
            label_text = sprintf('%% from median\n%g', v);
        otherwise
            label_text = sprintf('%s %g', xname, v);
    end
end
