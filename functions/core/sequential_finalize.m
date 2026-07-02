function sequential_finalize(all_run_params, run_affixes, options)
% SEQUENTIAL_FINALIZE  Post-chain outputs for a sequential simulation run
%
% Generates the multi-run summary report, the voxelwise-summed combined
% NIfTIs, and (optionally) cleans up per-run intermediates, once every run in
% a sequential chain has completed.
%
% This is invoked from prestus_pipeline STAGE 12 of the LAST run in the chain
% (not from the run that dispatches it).  Running it there guarantees that
% every config 1..N has already written its NIfTIs — critical under slurm,
% where each follow-up run is submitted as a separate non-blocking job, so
% finalising from the dispatching run would execute before the final run's
% outputs exist.
%
% Use as:
%   sequential_finalize(all_run_params, run_affixes, options)
%
% Input:
%   all_run_params - (1xN) cell array of PRESTUS parameters structs, one per
%                    run (base run first, then each sequential follow-up).
%   run_affixes    - (1xN) cell array of structs, each with fields
%                    output_affix and thermal_cache_affix; the affixes
%                    submitted for runs 1..N.
%   options        - pipeline options struct.  Recognised fields:
%                    .sequential_combine_niftis        (default true)
%                    .sequential_combine_units         (see combine_sequential_niftis)
%                    .sequential_cleanup_intermediate  (default false)
%
% See also: SEQUENTIAL_PIPELINE, GENERATE_SEQUENTIAL_REPORT,
%           COMBINE_SEQUENTIAL_NIFTIS

arguments
    all_run_params (1,:) cell
    run_affixes    (1,:) cell   = {}
    options        (1,1) struct = struct()
end

    if numel(all_run_params) <= 1
        return
    end

    % ---- multi-run summary report ----
    report_ok = false;
    try
        if numel(run_affixes) == numel(all_run_params)
            seq_labels = cellfun(@(a) label_from_affix(a.output_affix), run_affixes, ...
                'UniformOutput', false);
        else
            seq_labels = {};
        end
        generate_sequential_report(all_run_params, seq_labels, run_affixes);
        report_ok = true;
    catch ME_rep
        warning('prestus_pipeline:sequentialReport', ...
            'Sequential report generation failed: %s', ME_rep.message);
    end

    % ---- combine per-run NIfTIs into voxelwise-summed maps ----
    % Runs before cleanup so the per-run source NIfTIs still exist; combined
    % maps are written to <dir_output>/combined and survive cleanup.
    % Enabled by default; set options.sequential_combine_niftis = false to skip.
    do_combine = ~isfield(options, 'sequential_combine_niftis') || ...
                 options.sequential_combine_niftis;
    if do_combine
        try
            combine_sequential_niftis(all_run_params, run_affixes, options);
        catch ME_comb
            warning('prestus_pipeline:sequentialCombine', ...
                'Sequential NIfTI combination failed: %s', ME_comb.message);
        end
    end

    % ---- optional per-run NIfTI / image cleanup ----
    % Only runs after a successful report so integrated outputs exist first.
    % Cache (including heating timeseries .mat) is always retained.
    do_cleanup = report_ok && ...
                 isfield(options, 'sequential_cleanup_intermediate') && ...
                 options.sequential_cleanup_intermediate;
    if do_cleanup
        for ri = 1:numel(all_run_params)
            p = all_run_params{ri};
            if isfield(p.io, 'dir_output')
                base = p.io.dir_output;
            else
                continue
            end
            for subdir = {fullfile(base, 'nii'), fullfile(base, 'img')}
                d = subdir{1};
                if isfolder(d)
                    fprintf('Removing intermediate outputs: %s\n', d);
                    rmdir(d, 's');
                end
            end
        end
    end
end

function lbl = label_from_affix(affix)
    if isempty(affix)
        lbl = 'Base';
    else
        lbl = strtrim(regexprep(affix, '^_', ''));
        if isempty(lbl); lbl = affix; end
    end
end
