function combined_paths = combine_sequential_niftis(run_params_list, run_affixes, options)
% COMBINE_SEQUENTIAL_NIFTIS  Voxelwise-sum per-run NIfTIs across a sequential chain
%
% After a chain of sequential simulations completes, this adds together the
% NIfTI output of every run into a single "combined" map per data type and per
% coordinate space, for one subject.  It is the PRESTUS-native equivalent of an
% `fslmaths <run1> -add <run2> -add … <combined>` sweep: each run is identified
% by its own io.output_affix, and the runs that share a data type / space are
% summed voxelwise.
%
% Unlike generate_sequential_report's voxelwise-MAX intensity map (worst-case
% single exposure), this produces the voxelwise SUM (cumulative deposition
% across the whole sequence).
%
% Combined files are written to a 'combined' subfolder of the base run's output
% directory (never into nii/), so they are not re-detected as per-run outputs
% and survive options.sequential_cleanup_intermediate:
%     <dir_output>/combined/sub-XXX_<medium>_<space>_seqcombined_<unit>.nii.gz
%
% Per-run filenames are located through the explicit affix registry
% (run_affixes) — the actual output affixes submitted for runs 1..N — exactly
% like generate_sequential_report.  This guarantees every submitted config is
% combined by its own affix, rather than re-inferring the name from each
% parameter struct.  When the registry is omitted, the affix falls back to each
% run's io.output_affix.
%
% Use as:
%   combined_paths = combine_sequential_niftis(run_params_list)
%   combined_paths = combine_sequential_niftis(run_params_list, run_affixes)
%   combined_paths = combine_sequential_niftis(run_params_list, run_affixes, options)
%
% Input:
%   run_params_list - (1xN) cell array of PRESTUS parameters structs, one per
%                     run (base run first, then each sequential follow-up).
%   run_affixes     - (optional) (1xN) cell array of structs, each with field
%                     .output_affix, giving the affix submitted for that run.
%                     When its length matches run_params_list it is the source
%                     of truth for filenames; otherwise per-run io.output_affix
%                     is used.
%   options         - (optional) struct.  Recognised fields:
%                     .sequential_combine_units - cellstr of data types to
%                        combine.  Default: physically-additive maps only
%                        (see DEFAULT_UNITS below).
%
% Output:
%   combined_paths  - cellstr of the combined NIfTI files that were written.
%
% See also: GENERATE_SEQUENTIAL_REPORT, NIFTI_ACOUSTIC, NIFTI_THERMAL,
%           SEQUENTIAL_PIPELINE

arguments
    run_params_list (1,:) cell
    run_affixes     (1,:) cell   = {}
    options         (1,1) struct = struct()
end

    combined_paths = {};

    % Data types to sum by default. Restricted to physically-additive maps:
    % absolute-temperature maps (heating / heating_end) are excluded because
    % summing them would stack the ~37 C baselines. Pass any of the other
    % PRESTUS units (MI, heating, heating_end, heatrise, CEM43, CEM43_iso) via
    % options.sequential_combine_units to override.
    DEFAULT_UNITS = {'intensity', 'pressure', 'heatrise_end', ...
                     'CEM43_end', 'CEM43_iso_end'};
    if isfield(options, 'sequential_combine_units') && ...
            ~isempty(options.sequential_combine_units)
        units = cellstr(options.sequential_combine_units);
    else
        units = DEFAULT_UNITS;
    end

    n_runs = numel(run_params_list);
    if n_runs < 2
        fprintf('combine_sequential_niftis: fewer than 2 runs — nothing to combine.\n');
        return
    end

    % Affix registry is authoritative when it covers every run; otherwise fall
    % back to each run's io.output_affix.
    has_affix_registry = numel(run_affixes) == n_runs;

    % Resolve io dirs so dir_nii_* / dir_output are populated even for params
    % built before the pipeline ran (mirrors generate_sequential_report).
    for ri = 1:n_runs
        run_params_list{ri} = resolve_io_dirs(run_params_list{ri});
    end

    base_p     = run_params_list{1};
    subject_id = base_p.subject_id;
    medium     = base_p.simulation.medium;

    % NIfTIs are only produced for layered / phantom media.
    if ~contains(medium, {'layered', 'phantom'})
        fprintf('combine_sequential_niftis: medium "%s" produces no NIfTIs — skipping.\n', medium);
        return
    end

    % Coordinate spaces: T1w always; MNI only for layered runs that saved MNI.
    spaces = struct('tag', {'T1w'}, 'dir_field', {'dir_nii_T1w'});
    if strcmp(medium, 'layered') && should_save_output(base_p.io, 'save_MNI')
        spaces(end+1) = struct('tag', 'MNI', 'dir_field', 'dir_nii_MNI');
    end

    % Output folder (outside nii/ so cleanup keeps it and the glob ignores it).
    combined_dir = fullfile(base_p.io.dir_output, 'combined');
    if ~exist(combined_dir, 'dir'); mkdir(combined_dir); end

    for si = 1:numel(spaces)
        space_tag = spaces(si).tag;
        dir_field = spaces(si).dir_field;

        for ui = 1:numel(units)
            unit = units{ui};

            % Collect the per-run files that actually exist on disk, one per
            % submitted affix (runs 1..N).
            run_files = {};
            for ri = 1:n_runs
                p = run_params_list{ri};
                if ~isfield(p.io, dir_field) || isempty(p.io.(dir_field))
                    continue
                end
                if has_affix_registry
                    o_affix = run_affixes{ri}.output_affix;
                else
                    o_affix = p.io.output_affix;
                end
                f = fullfile(p.io.(dir_field), ...
                    sprintf('sub-%03d_%s_%s%s_%s.nii.gz', ...
                        subject_id, medium, space_tag, o_affix, unit));
                if isfile(f)
                    run_files{end+1} = f; %#ok<AGROW>
                end
            end

            if numel(run_files) < 2
                % Not a genuine combination (0 or 1 run has this map).
                continue
            end

            % Voxelwise sum.
            vol_sum = [];
            hdr_ref = [];
            n_added = 0;
            for fi = 1:numel(run_files)
                try
                    vol = single(niftiread(run_files{fi}));
                catch ME
                    warning('combine_sequential_niftis:read', ...
                        'Could not read %s: %s', run_files{fi}, ME.message);
                    continue
                end
                if isempty(vol_sum)
                    vol_sum = vol;
                    hdr_ref = niftiinfo(run_files{fi});
                    n_added = 1;
                elseif isequal(size(vol), size(vol_sum))
                    vol_sum = vol_sum + vol;
                    n_added = n_added + 1;
                else
                    warning('combine_sequential_niftis:size', ...
                        'Grid mismatch, skipping %s (size %s vs %s).', ...
                        run_files{fi}, mat2str(size(vol)), mat2str(size(vol_sum)));
                end
            end

            if isempty(vol_sum) || n_added < 2
                continue
            end

            out_path = fullfile(combined_dir, ...
                sprintf('sub-%03d_%s_%s_seqcombined_%s.nii.gz', ...
                    subject_id, medium, space_tag, unit));
            try
                hdr_ref.Filename = out_path;
                hdr_ref.Datatype = 'single';
                niftiwrite(vol_sum, strrep(out_path, '.nii.gz', '.nii'), ...
                    hdr_ref, 'Compressed', true);
                combined_paths{end+1} = out_path; %#ok<AGROW>
                fprintf('Combined %d runs -> %s\n', n_added, out_path);
            catch ME
                warning('combine_sequential_niftis:write', ...
                    'Could not write %s: %s', out_path, ME.message);
            end
        end
    end

    if isempty(combined_paths)
        fprintf('combine_sequential_niftis: no combinable NIfTI sets found.\n');
    end
end
