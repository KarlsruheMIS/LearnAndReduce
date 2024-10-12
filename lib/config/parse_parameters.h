/*****************************************************************************
 * parse_parameters.h
 * Purpose: Parse command line parameters.
 *
 *****************************************************************************/

#pragma once
#include <omp.h>

#include "configuration_mis.h"
#include "configuration_reduction.h"
#include "argtable3.h"
#include "string.h"

class BaseArguments
{
public:
    BaseArguments(int argc, char **argv) : argc(argc), argv(argv)
    {
        progname = argv[0];

        // Common parameters
        help = arg_lit0(NULL, "help", "Print help.");
        user_seed = arg_int0(NULL, "seed", NULL, "Seed to use for the PRNG.");
        verbose = arg_lit0(NULL, "verbose", "Print detailed information.");
        filename = arg_strn(NULL, NULL, "FILE", 1, 1, "Path to graph file.");
        output = arg_str0(NULL, "output", NULL, "Path to store resulting independent set.");
        time_limit = arg_dbl0(NULL, "time_limit", NULL, "Total time limit in s. Default 1000s.");
        console_log = arg_lit0(NULL, "console_log", "Stream the log into the console");
        disable_checks = arg_lit0(NULL, "disable_checks", "Disable sortedness check during I/O.");
        weight_source = arg_str0(NULL, "weight_source", NULL, "Choose how the weights are assigned. Can be either: file (default), hybrid, uniform, geometric, unit.");
        end = arg_end(100);
    }

    bool checkHelp(void *argtable[])
    {
        if (help->count > 0)
        {
            printf("Usage: %s", progname);
            arg_print_syntax(stdout, argtable, "\n");
            arg_print_glossary(stdout, argtable, "  %-40s %s\n");
            return true;
        }
        return false;
    }
    void checkErrors(int nerrors)
    {
        if (nerrors > 0)
        {
            arg_print_errors(stdout, end, progname);
            printf("Try '%s --help' for more information.\n", progname);
            exit(EXIT_FAILURE);
        }
    }

    int setConfig(Config &config);
    void parseParameters(Config &config);

protected:
    int argc;
    char **argv;
    const char *progname;
    struct arg_lit *help;
    struct arg_int *user_seed;
    struct arg_str *filename;
    struct arg_str *output;
    struct arg_dbl *time_limit;
    struct arg_lit *console_log;
    struct arg_lit *verbose;
    struct arg_lit *disable_checks;
    struct arg_str *weight_source;
    struct arg_end *end;
};

class ReductionArguments : public BaseArguments
{
public:
    ReductionArguments(int argc, char **argv) : BaseArguments(argc, argv)
    {

        // Unique parameters
        reduction_style = arg_str0(NULL, "reduction_style", NULL, "Choose the type of reductions appropriate for the input graph. Can be either: full, normal, dense, test1, test2.");
        // heuristic_style             = arg_str0(NULL, "heuristic_style", NULL, "Configuration to use in the heuristic reductions. ([single, multiple_very_safe, multiple_safe, all (if no gnn_filter this is the same as multiple_safe)]). Default: multiple_safe.");
        reduction_time_limit = arg_dbl0(NULL, "reduction_time_limit", NULL, "Time limit for reduction in s. Default equal to overall time limit.");
        reduction_config = arg_str0(NULL, "reduction_config", NULL, "Configuration to use. ([cyclicFast, cyclicStrong, kamis, mmwis, all_reductions_cyclicFast, all_reductions_CyclicStrong, extended_cyclicFast, all_decreasing, fast, very_fast]). Default: decreasing (not using increasing reductions).");
        kernel_filename = arg_str0(NULL, "kernel", NULL, "Path to store resulting kernel.");
        disable_reduction = arg_lit0(NULL, "disable_reduction", "Don't perforn any reductions.");
        print_reduction_info = arg_lit0(NULL, "print_reduction_info", "Print detailed information about each reduction");
        disable_early_termination = arg_lit0(NULL, "disable_early_termination", "Disable early termination of solving subgraphs in reductions.");
        gnn_filter = arg_str0(NULL, "gnn_filter", NULL, "Choose gnn filter strategy. Can be either: never (default), initial, initial_tight, always");

        // single reduction parameters
        disable_neighborhood = arg_lit0(NULL, "disable_neighborhood", "Disable neighborhood reduction.");
        disable_fold1 = arg_lit0(NULL, "disable_fold1", "Disable fold1 reduction.");
        disable_fold2 = arg_lit0(NULL, "disable_fold2", "Disable fold2 reduction.");
        disable_triangle_min = arg_lit0(NULL, "disable_triangle_min", "Disable triangle min reduction.");
        disable_triangle_mid = arg_lit0(NULL, "disable_triangle_mid", "Disable triangle mid reduction.");
        disable_v_shape_min = arg_lit0(NULL, "disable_v_shape_min", "Disable v_shape_min reduction.");
        disable_v_shape_max = arg_lit0(NULL, "disable_v_shape_max", "Disable v_shape_max reduction.");
        disable_v_shape_mid = arg_lit0(NULL, "disable_v_shape_mid", "Disable v_shape_mid reduction.");
        disable_domination = arg_lit0(NULL, "disable_domination", "Disable domination reduction.");
        disable_extended_domination = arg_lit0(NULL, "disable_extended_domination", "Disable extended domination reduction.");
        disable_extended_domination_reverse = arg_lit0(NULL, "disable_extended_domination_reverse", "Disable extended domination reverse reduction.");
        disable_extended_twin = arg_lit0(NULL, "disable_extended_twin", "Disable extended twin reduction.");
        disable_basic_se = arg_lit0(NULL, "disable_basic_se", "Disable basic single edge reduction.");
        disable_extended_se = arg_lit0(NULL, "disable_extended_se", "Disable extended single edge reduction.");
        disable_clique_neighborhood = arg_lit0(NULL, "disable_clique_neighborhood", "Disable neighborhood clique reduction.");
        disable_clique_neighborhood_fast = arg_lit0(NULL, "disable_clique_neighborhood_fast", "Disable neighborhood clique fast reduction.");
        disable_generalized_fold = arg_lit0(NULL, "disable_generalized_fold", "Disable generalized fold reduction.");
        disable_critical_set = arg_lit0(NULL, "disable_critical_set", "Disable critical set reduction.");
        disable_clique = arg_lit0(NULL, "disable_clique", "Disable clique reduction.");
        disable_twin = arg_lit0(NULL, "disable_twin", "Disable twin reduction.");
        disable_heavy_set = arg_lit0(NULL, "disable_heavy_set", "Disable heavy set reduction.");
        disable_heavy_set3 = arg_lit0(NULL, "disable_heavy_set3", "Disable heavy set3 reduction.");
        disable_cut_vertex = arg_lit0(NULL, "disable_cut_vertex", "Disable cut vertex reduction.");
        disable_unconfined = arg_lit0(NULL, "disable_unconfined", "Disable unconfined reduction.");
        disable_funnel = arg_lit0(NULL, "disable_funnel", "Disable funnel reduction.");
        disable_funnel_fold = arg_lit0(NULL, "disable_funnel_fold", "Disable funnel fold reduction.");
        // disable_high_degree         = arg_lit0(NULL, "disable_high_degree", "Disable high degree reduction.");
        // disable_bound               = arg_lit0(NULL, "disable_boundary", "Disable boundary reduction.");
        // disable_heuristic_include   = arg_lit0(NULL, "disable_heuristic_include", "Disable heuristic include reduction.");
        // disable_heursitic_exclude   = arg_lit0(NULL, "disable_heuristic_exclude", "Disable heuristic exclude reduction.");
        disable_decreasing_struction = arg_lit0(NULL, "disable_decreasing_struction", "Disable decreasing struction.");
        disable_plateau_struction = arg_lit0(NULL, "disable_plateau_struction", "Disable plateau struction.");
        subgraph_node_limit = arg_int0(NULL, "subgraph_node_limit", NULL, "Choose maximum number of nodes in subgraph.");

        // for struction and branch_and_reduce
        set_limit = arg_int0(NULL, "set_limit", NULL, "Choose maximum number of new vertices allowed to be created during struction application");
        struction_degree = arg_int0(NULL, "struction_degree", NULL, "Choose maximum degree of vertex to perform struction on it.");
        struction_type = arg_str0(NULL, "struction_type", NULL, "Choose used struction type. Can be either: original, modified, extended (default), extended_reduced or none");
        key_type = arg_str0(NULL, "key_type", NULL, "Choose used vertex selection strategy. Can be either: random, degree, increase (default), approximate_increase");
        key_reinsert_factor = arg_dbl0(NULL, "key_reinsert_factor", NULL, "Choose reinsert factor beta for approximate increase selection. default 2");
        global_blow_up_percent = arg_dbl0(NULL, "global_blow_up_factor", NULL, "Choose global blow up factor alpha.");
        phase_blow_up_factor = arg_dbl0(NULL, "phase_blow_up_factor", NULL, "Choose percentual phase blow up factor gamma.");
        phase_blow_ups = arg_int0(NULL, "phase_blow_ups", NULL, "Choose fixed blow ups per phase Z.");
        max_unimproving_phases = arg_int0(NULL, "max_unimproving_phases", NULL, "Choose maximum unimproving blow up phases X.");
        backtrack_style = arg_str0(NULL, "backtrack_style", NULL, "Choose backtrack strategy. Can be either: immediate_tie_breaking, immediate_exclude (default), end_min_kernel or no_backtrack");
        disable_blow_up = arg_lit0(NULL, "disable_blow_up", "Disable cyclic blow up algorithm.");
        plain_struction = arg_lit0(NULL, "plain_struction", "Only use struction to reduce graph.");
    }

    int setConfig(ReductionConfig &config);
    void parseParameters(ReductionConfig &config);

protected:
    struct arg_str *reduction_style;
    struct arg_dbl *reduction_time_limit;
    struct arg_lit *disable_early_termination;
    struct arg_str *gnn_filter;
    struct arg_str *reduction_config;
    // struct arg_str * heuristic_style;
    struct arg_str *kernel_filename;
    struct arg_lit *print_reduction_info;
    struct arg_lit *disable_reduction;
    struct arg_lit *disable_neighborhood;
    struct arg_lit *disable_fold1;
    struct arg_lit *disable_fold2;
    struct arg_lit *disable_triangle_min;
    struct arg_lit *disable_triangle_mid;
    struct arg_lit *disable_v_shape_min;
    struct arg_lit *disable_v_shape_max;
    struct arg_lit *disable_v_shape_mid;
    struct arg_lit *disable_domination;
    struct arg_lit *disable_extended_domination;
    struct arg_lit *disable_extended_domination_reverse;
    struct arg_lit *disable_extended_twin;
    struct arg_lit *disable_basic_se;
    struct arg_lit *disable_extended_se;
    struct arg_lit *disable_clique_neighborhood;
    struct arg_lit *disable_clique_neighborhood_fast;
    struct arg_lit *disable_generalized_fold;
    struct arg_lit *disable_critical_set;
    struct arg_lit *disable_clique;
    struct arg_lit *disable_twin;
    struct arg_lit *disable_heavy_set;
    struct arg_lit *disable_heavy_set3;
    struct arg_lit *disable_cut_vertex;
    struct arg_lit *disable_unconfined;
    struct arg_lit *disable_funnel;
    struct arg_lit *disable_funnel_fold;
    // struct arg_lit * disable_high_degree;
    // struct arg_lit * disable_bound;
    struct arg_lit *disable_decreasing_struction;
    struct arg_lit *disable_plateau_struction;
    // struct arg_lit * use_heuristic_reductions;
    // struct arg_lit * disable_heuristic_include;
    // struct arg_lit * disable_heursitic_exclude;
    struct arg_int *subgraph_node_limit;
    struct arg_int *struction_degree;
    struct arg_int *set_limit;
    struct arg_str *struction_type;
    struct arg_str *key_type;
    struct arg_dbl *key_reinsert_factor;
    struct arg_dbl *global_blow_up_percent;
    struct arg_dbl *phase_blow_up_factor;
    struct arg_int *phase_blow_ups;
    struct arg_int *max_unimproving_phases;
    struct arg_str *backtrack_style;
    struct arg_lit *plain_struction;
    struct arg_lit *disable_blow_up;
};

int BaseArguments::setConfig(Config &config)
{
    // Common argtable
    void *argtable[] = {
        help,
        user_seed,
        filename,
        output,
        time_limit,
        console_log,
        disable_checks,
        weight_source,
        nullptr};

    // Parse the arguments
    int nerrors = arg_parse(argc, argv, argtable);
    checkErrors(nerrors);
    if (checkHelp(argtable))
    {
        arg_freetable(argtable, sizeof(argtable) / sizeof(argtable[0]));
        exit(EXIT_SUCCESS);
    }
    parseParameters(config);
    arg_freetable(argtable, sizeof(argtable) / sizeof(argtable[0]));
    return 0;
}
void BaseArguments::parseParameters(Config &config)
{
    // Choose standard configuration
    configuration_mis cfg;
    cfg.standard(config);

    if (user_seed->count > 0)
    {
        config.seed = user_seed->ival[0];
    }

    if (filename->count > 0)
    {
        config.graph_filename = filename->sval[0];
    }

    if (output->count > 0)
    {
        config.output_filename = output->sval[0];
    }

    if (time_limit->count > 0)
    {
        config.time_limit = time_limit->dval[0];
    }

    if (console_log->count > 0)
    {
        config.console_log = true;
        config.print_log = false;
    }
    else
    {
        config.print_log = true;
    }

    if (verbose->count > 0)
    {
        config.verbose = true;
    }
    else
    {
        config.verbose = false;
        config.console_log = false;
        config.print_log = false;
    }

    if (disable_checks->count > 0)
    {
        config.check_sorted = false;
    }

    if (weight_source->count > 0)
    {
        config.setWeightSource(weight_source->sval[0]);
    }
}

int ReductionArguments::setConfig(ReductionConfig &config)
{
    // Common argtable
    void *argtable[] = {
        help,
        user_seed,
        filename,
        output,
        time_limit,
        console_log,
        verbose,
        disable_checks,
        weight_source,
        print_reduction_info,
        disable_early_termination,
        gnn_filter,
        // use_heuristic_reductions,
        kernel_filename,
        reduction_config,
        // heuristic_style,
        reduction_time_limit,
        reduction_style,
        disable_reduction,
        disable_neighborhood,
        disable_fold1,
        disable_fold2,
        disable_triangle_min,
        disable_triangle_mid,
        disable_v_shape_min,
        disable_v_shape_max,
        disable_v_shape_mid,
        disable_domination,
        disable_extended_domination,
        disable_extended_domination_reverse,
        disable_extended_twin,
        disable_basic_se,
        disable_extended_se,
        disable_clique_neighborhood,
        disable_clique_neighborhood_fast,
        disable_generalized_fold,
        disable_critical_set,
        disable_clique,
        disable_twin,
        disable_heavy_set,
        disable_heavy_set3,
        disable_cut_vertex,
        disable_unconfined,
        disable_funnel,
        disable_funnel_fold,
        // disable_high_degree,
        // disable_bound,
        disable_decreasing_struction,
        disable_plateau_struction,
        // disable_heuristic_include,
        // disable_heursitic_exclude,
        subgraph_node_limit,
        struction_degree,
        set_limit,
        struction_type,
        key_type,
        key_reinsert_factor,
        global_blow_up_percent,
        phase_blow_up_factor,
        phase_blow_ups,
        max_unimproving_phases,
        backtrack_style,
        plain_struction,
        disable_blow_up,
        end};

    // Parse the arguments
    int nerrors = arg_parse(argc, argv, argtable);
    checkErrors(nerrors);
    if (checkHelp(argtable))
    {
        arg_freetable(argtable, sizeof(argtable) / sizeof(argtable[0]));
        exit(EXIT_SUCCESS);
    }
    parseParameters(config);
    arg_freetable(argtable, sizeof(argtable) / sizeof(argtable[0]));
    return 0;
}
void ReductionArguments::parseParameters(ReductionConfig &config)
{

    // Choose standard configuration
    configuration_reduction cfg;
    cfg.all_decreasing(config);
    // Parse the arguments
    if (reduction_config->count > 0)
    {
        if (!strcmp(reduction_config->sval[0], "cyclicFast"))
        {
            cfg.original_cyclicFast(config);
            config.reduction_config_name = "cyclicFast";
        }
        else if (!strcmp(reduction_config->sval[0], "fast"))
        {
            cfg.fast(config);
            config.reduction_config_name = "fast";
        }
        else if (!strcmp(reduction_config->sval[0], "very_fast"))
        {
            cfg.very_fast(config);
            config.reduction_config_name = "very_fast";
        }
        else if (!strcmp(reduction_config->sval[0], "fast_reductions_cyclicFast"))
        {
            cfg.fast_reductions_cyclicFast(config);
            config.reduction_config_name = "fast_reductions_cyclicFast";
        }
        else if (!strcmp(reduction_config->sval[0], "all_reductions_cyclicFast"))
        {
            cfg.all_reductions_cyclicFast(config);
            config.reduction_config_name = "all_reductions_cyclicFast";
        }
        else if (!strcmp(reduction_config->sval[0], "extended_cyclicFast"))
        {
            cfg.extended_cyclicFast(config);
            config.reduction_config_name = "extended_cyclicFast";
        }
        else if (!strcmp(reduction_config->sval[0], "cyclicStrong"))
        {
            cfg.original_cyclicStrong(config);
            config.reduction_config_name = "cyclicStrong";
        }
        else if (!strcmp(reduction_config->sval[0], "all_reductions_cyclicStrong"))
        {
            cfg.all_reductions_cyclicStrong(config);
            config.reduction_config_name = "all_reductions_cyclicStrong";
        }
        else if (!strcmp(reduction_config->sval[0], "kamis"))
        {
            cfg.original_kamis(config);
            config.reduction_config_name = "kamis";
        }
        else if (!strcmp(reduction_config->sval[0], "mmwis"))
        {
            cfg.mmwis(config);
            config.reduction_config_name = "mmwis";
        }
        else
        {
            cfg.all_decreasing(config);
            config.reduction_config_name = "all_decreasing";
        }
    }

    // if (heuristic_style->count > 0) {
    //     if (!strcmp(heuristic_style->sval[0], "single"))
    //     {
    //         config.heuristic_style_name = "single";
    //         config.heuristic_style      = ReductionConfig::Heuristic_Style::single;
    //     }
    //     else if (!strcmp(heuristic_style->sval[0], "all"))
    //     {
    //         config.heuristic_style_name = "all";
    //         config.heuristic_style      = ReductionConfig::Heuristic_Style::all;
    //     }
    //     else if (!strcmp(heuristic_style->sval[0], "multiple_very_safe"))
    //     {
    //         config.heuristic_style_name = "multiple_very_safe";
    //         config.heuristic_style      = ReductionConfig::Heuristic_Style::multiple_very_safe;
    //     }
    //     else { // default
    //         config.heuristic_style_name = "multiple_safe";
    //         config.heuristic_style      = ReductionConfig::Heuristic_Style::multiple_safe;
    //     }
    // }

    BaseArguments::parseParameters(config);

    if (reduction_style->count > 0)
    {
        config.setReductionStyle(reduction_style->sval[0]);
    }
    if (reduction_time_limit->count > 0)
    {
        config.reduction_time_limit = reduction_time_limit->dval[0];
    }
    else
    {
        config.reduction_time_limit = config.time_limit;
    }
    if (disable_early_termination->count > 0)
    {
        config.disable_early_termination = true;
    }
    else
    {
        config.disable_early_termination = false;
    }
    if (gnn_filter->count > 0)
    {
        config.setGNNFilterStyle(gnn_filter->sval[0]);
    }
    if (print_reduction_info->count > 0)
    {
        config.print_reduction_info = true;
    }
    if (disable_neighborhood->count > 0)
    {
        config.disable_neighborhood = true;
    }
    if (disable_fold1->count > 0)
    {
        config.disable_fold1 = true;
    }
    if (disable_fold2->count > 0)
    {
        config.disable_fold2 = true;
    }
    if (disable_triangle_min->count > 0)
    {
        config.disable_triangle_min = true;
    }
    if (disable_triangle_mid->count > 0)
    {
        config.disable_triangle_mid = true;
    }
    if (disable_v_shape_min->count > 0)
    {
        config.disable_v_shape_min = true;
    }
    if (disable_v_shape_max->count > 0)
    {
        config.disable_v_shape_max = true;
    }
    if (disable_v_shape_mid->count > 0)
    {
        config.disable_v_shape_mid = true;
    }
    if (disable_domination->count > 0)
    {
        config.disable_domination = true;
    }
    if (disable_extended_domination->count > 0)
    {
        config.disable_extended_domination = true;
    }
    if (disable_extended_domination_reverse->count > 0)
    {
        config.disable_extended_domination_reverse = true;
    }
    if (disable_extended_twin->count > 0)
    {
        config.disable_extended_twin = true;
    }
    if (disable_basic_se->count > 0)
    {
        config.disable_basic_se = true;
    }
    if (disable_extended_se->count > 0)
    {
        config.disable_extended_se = true;
    }
    if (disable_clique_neighborhood->count > 0)
    {
        config.disable_clique_neighborhood = true;
    }
    if (disable_clique_neighborhood_fast->count > 0)
    {
        config.disable_clique_neighborhood_fast = true;
    }
    if (disable_generalized_fold->count > 0)
    {
        config.disable_generalized_fold = true;
    }
    if (disable_critical_set->count > 0)
    {
        config.disable_critical_set = true;
    }
    if (disable_clique->count > 0)
    {
        config.disable_clique = true;
    }
    if (disable_twin->count > 0)
    {
        config.disable_twin = true;
    }
    if (disable_heavy_set->count > 0)
    {
        config.disable_heavy_set = true;
    }
    if (disable_heavy_set3->count > 0)
    {
        config.disable_heavy_set3 = true;
    }
    if (disable_cut_vertex->count > 0)
    {
        config.disable_cut_vertex = true;
    }
    if (disable_unconfined->count > 0)
    {
        config.disable_unconfined = true;
    }
    if (disable_funnel_fold->count > 0)
    {
        config.disable_funnel_fold = true;
    }
    if (disable_funnel->count > 0)
    {
        config.disable_funnel = true;
    }
    // if (disable_high_degree->count > 0) {
    //     config.disable_high_degree = true;
    // }
    // if (disable_bound->count > 0) {
    //     config.disable_bound_reduction = true;
    // }
    // if (use_heuristic_reductions->count > 0) {
    //     config.disable_heuristic_exclude = false;
    //     config.disable_heuristic_include = false;
    //     if (config.heuristic_style == ReductionConfig::Heuristic_Style::none)
    //     {
    //         config.heuristic_style      = ReductionConfig::Heuristic_Style::single;
    //         config.heuristic_style_name = "single";
    //     }
    // }
    // if (disable_heuristic_include->count > 0) {
    //     config.disable_heuristic_include = true;
    // }
    // if (disable_heursitic_exclude->count > 0) {
    //     config.disable_heuristic_exclude = true;
    // }
    if (subgraph_node_limit->count > 0)
    {
        config.subgraph_node_limit = subgraph_node_limit->ival[0];
    }

    if (kernel_filename->count > 0)
    {
        config.kernel_filename = kernel_filename->sval[0];
        config.write_kernel = true;
    }
    else
    {
        config.write_kernel = false;
    }

    if (disable_decreasing_struction->count > 0)
    {
        config.disable_decreasing_struction = true;
    }
    if (disable_plateau_struction->count > 0)
    {
        config.disable_plateau_struction = true;
    }
    if (struction_degree->count > 0)
    {
        config.struction_degree = struction_degree->ival[0];
    }
    if (set_limit->count > 0)
    {
        config.set_limit = set_limit->ival[0];
    }
    if (struction_type->count > 0)
    {
        config.setStructionType(struction_type->sval[0]);
    }
    if (key_type->count > 0)
    {
        config.setKeyType(key_type->sval[0]);
    }
    if (key_reinsert_factor->count > 0)
    {
        config.key_reinsert_factor = key_reinsert_factor->dval[0];
    }
    if (global_blow_up_percent->count > 0)
    {
        config.global_blow_up_factor = global_blow_up_percent->dval[0];
    }
    if (phase_blow_up_factor->count > 0)
    {
        config.phase_blow_up_factor = phase_blow_up_factor->dval[0];
    }
    if (phase_blow_ups->count > 0)
    {
        config.phase_blow_ups = phase_blow_ups->ival[0];
    }
    if (max_unimproving_phases->count > 0)
    {
        config.max_unimproving_phases = max_unimproving_phases->ival[0];
    }
    if (backtrack_style->count > 0)
    {
        config.setBacktrackType(backtrack_style->sval[0]);
    }
    if (plain_struction->count)
    {
        config.plain_struction = true;
    }
    if (disable_blow_up->count > 0)
    {
        config.disable_blow_up = true;
    }
}
