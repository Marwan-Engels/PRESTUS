function combined_paths = combine_sequential_niftis(run_params_list, run_affixes, options)
% COMBINE_SEQUENTIAL_NIFTIS  Voxelwise-sum NIfTI outputs across a sequential chain
%
% After a sequential chain of PRESTUS runs completes, this sums each data type's
% NIfTI across all runs in the chain, grouped by hemisphere, and writes the result
% to <dir_output>/nii/sequential_nii/. Replaces the old workflow of running
% C_CtrlTUS_combine_subject_sims_v5.m by hand with fslmaths -add against a
% hardcoded affix list: here the runs and their affixes come directly from the
% chain that was just simulated (first_config + options.sequential_configs, in
% the ascending order sequential_pipeline already dispatches them in), so no
% affix list needs to be maintained separately.
%
% Grouping: hemisphere is not a field on the parameters struct — it only exists
% as an 'L'/'R' token embedded in each run's io.output_affix (e.g.
% '_L_pgACC_1_pos1_F45_I600_r2mm_CtrlTUS_posthoc'). Runs are grouped by the first
% matching token in options.sequential_combine_hemisphere_tokens found in their
% output_affix; runs matching no token fall into a single 'all' group. If no run
% in the whole chain matches any token, the entire chain is treated as one
% ungrouped chain (so this also works for sequential chains that don't follow the
% CtrlTUS L/R convention at all).
%
% Use as:
%   combined_paths = combine_sequential_niftis(run_params_list, run_affixes)
%   combined_paths = combine_sequential_niftis(run_params_list, run_affixes, options)
%
% Input:
%   run_params_list - (1xN) cell array of PRESTUS parameters structs, one per
%                      run in the chain, in ascending order (as accumulated by
%                      sequential_pipeline.m into options.sequential_report_runs)
%   run_affixes     - (1xN) cell array of structs with field .output_affix, one
%                      per run, same order as run_params_list (as accumulated by
%                      sequential_pipeline.m into options.sequential_run_affixes)
%   options         - (optional) struct with fields:
%     .sequential_combine_datatypes        - cell array of data type name strings.
%                                             Default: {'intensity','MI','pressure',
%                                             'heating_end','CEM43_end','CEM43_iso_end'}
%     .sequential_combine_hemisphere_tokens - cell array of tokens to look for as
%                                             '_<token>_' in output_affix.
%                                             Default: {'L','R'}
%
% Output:
%   combined_paths - struct array with fields .hemisphere, .data_type, .space,
%                     .path — one entry per combined NIfTI actually written
%
% See also: GENERATE_SEQUENTIAL_REPORT, SEQUENTIAL_PIPELINE, NIFTI_ACOUSTIC, NIFTI_THERMAL

arguments
    run_params_list (1,:) cell
    run_affixes     (1,:) cell
    options         (1,1) struct = struct()
end

combined_paths = struct('hemisphere', {}, 'data_type', {}, 'space', {}, 'path', {});

n_runs = numel(run_params_list);
if n_runs == 0
    warn('combine_sequential_niftis:empty', 'run_params_list is empty — nothing to combine.');
    return
end
if numel(run_affixes) ~= n_runs
    warn('combine_sequential_niftis:mismatch', ...
        'run_affixes must have the same number of entries as run_params_list.');
    return
end

if isfield(options, 'sequential_combine_datatypes') && ~isempty(options.sequential_combine_datatypes)
    data_types = options.sequential_combine_datatypes;
else
    data_types = {'intensity', 'MI', 'pressure', 'heating_end', 'CEM43_end', 'CEM43_iso_end'};
end

if isfield(options, 'sequential_combine_hemisphere_tokens') && ~isempty(options.sequential_combine_hemisphere_tokens)
    hemi_tokens = options.sequential_combine_hemisphere_tokens;
else
    hemi_tokens = {'L', 'R'};
end

base_p     = run_params_list{1};
subject_id = base_p.subject_id;
medium     = base_p.simulation.medium;

% ---- group runs by hemisphere token found in their output_affix ----
run_hemi = cell(1, n_runs);
for ri = 1:n_runs
    affix = run_affixes{ri}.output_affix;
    run_hemi{ri} = '';
    for ti = 1:numel(hemi_tokens)
        tok = hemi_tokens{ti};
        if contains(affix, ['_' tok '_'])
            run_hemi{ri} = tok;
            break
        end
    end
end

if all(cellfun(@isempty, run_hemi))
    % No run matched any hemisphere token: treat the whole chain as one group.
    run_hemi(:) = {'all'};
end
[run_hemi{cellfun(@isempty, run_hemi)}] = deal('all');

hemi_groups = unique(run_hemi, 'stable');

% ---- output directory ----
if isfield(base_p.io, 'dir_output') && ~isempty(base_p.io.dir_output)
    out_dir = fullfile(base_p.io.dir_output, 'nii', 'sequential_nii');
else
    warn('combine_sequential_niftis:noOutputDir', ...
        'base run has no io.dir_output set — cannot determine where to write combined NIfTIs.');
    return
end
if ~isfolder(out_dir); mkdir(out_dir); end

spaces = {'T1w', 'MNI'};

for hi = 1:numel(hemi_groups)
    hemi = hemi_groups{hi};
    group_idx = find(strcmp(run_hemi, hemi));
    if numel(group_idx) < 2
        continue % nothing to sum for a lone run
    end

    common_suffix = derive_common_suffix(run_affixes{group_idx(1)}.output_affix);

    for si = 1:numel(spaces)
        space = spaces{si};
        for di = 1:numel(data_types)
            dtype = data_types{di};

            vols_found = 0;
            vol_sum    = [];
            hdr_ref    = [];
            for gi = 1:numel(group_idx)
                ri = group_idx(gi);
                p  = run_params_list{ri};
                affix = run_affixes{ri}.output_affix;

                if strcmp(space, 'T1w')
                    dir_nii = p.io.dir_nii_T1w;
                else
                    dir_nii = p.io.dir_nii_MNI;
                end
                nii_path = fullfile(dir_nii, ...
                    sprintf('sub-%03d_%s_%s%s_%s.nii.gz', subject_id, medium, space, affix, dtype));

                if ~isfile(nii_path)
                    continue
                end

                try
                    vol = double(niftiread(nii_path));
                catch ME
                    warn('combine_sequential_niftis:read', ...
                        'Could not read %s: %s', nii_path, ME.message);
                    continue
                end

                if isempty(vol_sum)
                    vol_sum = vol;
                    hdr_ref = niftiinfo(nii_path);
                else
                    vol_sum = vol_sum + vol;
                end
                vols_found = vols_found + 1;
            end

            if vols_found < 2
                continue % need at least 2 runs' worth of data to call it "combined"
            end

            out_name = sprintf('sub-%03d_%s_%s_sequential_%s_combined%s_%s.nii.gz', ...
                subject_id, medium, space, hemi, common_suffix, dtype);
            out_path = fullfile(out_dir, out_name);

            try
                hdr_ref.Filename = out_path;
                hdr_ref.Datatype = 'single';
                niftiwrite(single(vol_sum), strrep(out_path, '.nii.gz', '.nii'), hdr_ref, 'Compressed', true);
                fprintf('Combined sequential NIfTI saved to: %s\n', out_path);

                combined_paths(end+1) = struct( ... %#ok<AGROW>
                    'hemisphere', hemi, 'data_type', dtype, 'space', space, 'path', out_path);
            catch ME
                warn('combine_sequential_niftis:write', ...
                    'Could not save combined NIfTI %s: %s', out_path, ME.message);
            end
        end
    end
end

end

function suffix = derive_common_suffix(affix)
% Everything from the first '_pos<N>_' onward is identical across all runs in a
% chain (only the target/position token varies per run). Falls back to the raw
% affix with a warning if that convention isn't present.
    tok = regexp(affix, '_pos\d+_(.+)$', 'tokens', 'once');
    if ~isempty(tok)
        suffix = ['_' tok{1}];
    else
        warn('combine_sequential_niftis:noPosToken', ...
            ['output_affix "%s" has no "_pos<N>_" token — falling back to the raw affix ' ...
             'for the combined filename. Collision-safety across parallel chains is not ' ...
             'guaranteed in this fallback case.'], affix);
        suffix = [affix '_combined'];
    end
end
