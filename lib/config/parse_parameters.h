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
        filename = arg_strn(NULL, NULL, "FILE", 1, 1, "Path to graph file.");
        user_seed = arg_int0(NULL, "seed", NULL, "Seed to use for the PRNG.");
        time_limit = arg_dbl0(NULL, "time_limit", NULL, "Total time limit in s. Default 1000s.");
        weight_source = arg_str0(NULL, "weight_source", NULL, "Choose how the weights are assigned. Can be either: file (default), hybrid, uniform, geometric, unit.");
        verbose = arg_lit0(NULL, "verbose", "Print detailed information.");
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
        if (nerrors > 0 && help->count == 0)
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
    struct arg_str *filename;
    struct arg_lit *help;
    struct arg_int *user_seed;
    struct arg_dbl *time_limit;
    struct arg_lit *verbose;
    struct arg_str *weight_source;
    struct arg_end *end;
};

class ReductionArguments : public BaseArguments
{
public:
    ReductionArguments(int argc, char **argv) : BaseArguments(argc, argv)
    {

        // Unique parameters
        print_reduction_info = arg_lit0(NULL, "print_reduction_info", "Print detailed information about each reduction.");
        reduction_config = arg_str0(NULL, "reduction_config", NULL, "What Reductions to use [no_gnn_reductions, all_reductions (default)].");
        kernel_filename = arg_str0(NULL, "kernel", NULL, "Path to store resulting kernel.");
        gnn_filter = arg_str0(NULL, "gnn_filter", NULL, "Choose gnn filter strategy. Can be either: never, initial, initial_tight (default), always.");
        cyclicFast = arg_lit0(NULL, "cyclicFast", "Set struction configuration cyclicFast.");
        cyclicStrong = arg_lit0(NULL, "cyclicStrong", "Set struction configuration cyclicStrong.");

        // single reduction parameters
        // disable_neighborhood = arg_lit0(NULL, "disable_neighborhood", "Disable neighborhood reduction.");
        // disable_degree_1 = arg_lit0(NULL, "disable_degree_1", "Disable degree_1 reduction.");
        // disable_degree_2 = arg_lit0(NULL, "disable_degree_2", "Disable degree_2 reduction.");
        // disable_triangle_min = arg_lit0(NULL, "disable_triangle_min", "Disable triangle min reduction.");
        // disable_triangle_mid = arg_lit0(NULL, "disable_triangle_mid", "Disable triangle mid reduction.");
        // disable_v_shape_min = arg_lit0(NULL, "disable_v_shape_min", "Disable v_shape_min reduction.");
        // disable_v_shape_max = arg_lit0(NULL, "disable_v_shape_max", "Disable v_shape_max reduction.");
        // disable_v_shape_mid = arg_lit0(NULL, "disable_v_shape_mid", "Disable v_shape_mid reduction.");
        // disable_domination = arg_lit0(NULL, "disable_domination", "Disable domination reduction.");
        // disable_extended_domination = arg_lit0(NULL, "disable_extended_domination", "Disable extended domination reduction.");
        // disable_extended_domination_reverse = arg_lit0(NULL, "disable_extended_domination_reverse", "Disable extended domination reverse reduction.");
        // disable_extended_twin = arg_lit0(NULL, "disable_extended_twin", "Disable extended twin reduction.");
        // disable_single_edge = arg_lit0(NULL, "disable_single_edge", "Disable basic single edge reduction.");
        // disable_extended_single_edge = arg_lit0(NULL, "disable_extended_single_edge", "Disable extended single edge reduction.");
        // disable_clique_neighborhood = arg_lit0(NULL, "disable_clique_neighborhood", "Disable neighborhood clique reduction.");
        // disable_clique_neighborhood_fast = arg_lit0(NULL, "disable_clique_neighborhood_fast", "Disable neighborhood clique fast reduction.");
        // disable_generalized_fold = arg_lit0(NULL, "disable_generalized_fold", "Disable generalized fold reduction.");
        // disable_critical_set = arg_lit0(NULL, "disable_critical_set", "Disable critical set reduction.");
        // disable_simplicial_vertex = arg_lit0(NULL, "disable_simplicial_vertex", "Disable clique reduction.");
        // disable_twin = arg_lit0(NULL, "disable_twin", "Disable twin reduction.");
        // disable_heavy_set = arg_lit0(NULL, "disable_heavy_set", "Disable heavy set reduction.");
        // disable_heavy_set3 = arg_lit0(NULL, "disable_heavy_set3", "Disable heavy set3 reduction.");
        // disable_cut_vertex = arg_lit0(NULL, "disable_cut_vertex", "Disable cut vertex reduction.");
        // disable_unconfined = arg_lit0(NULL, "disable_unconfined", "Disable unconfined reduction.");
        // disable_funnel = arg_lit0(NULL, "disable_funnel", "Disable funnel reduction.");
        // disable_struction_decrease = arg_lit0(NULL, "disable_struction_decrease", "Disable decreasing struction.");
        // disable_struction_plateau = arg_lit0(NULL, "disable_struction_plateau", "Disable plateau struction.");
    }

    int setConfig(ReductionConfig &config);
    void parseParameters(ReductionConfig &config);

protected:
    struct arg_lit *cyclicFast;
    struct arg_lit *cyclicStrong;
    struct arg_str *gnn_filter;
    struct arg_str *reduction_config;
    struct arg_str *kernel_filename;
    struct arg_lit *print_reduction_info;
    // struct arg_lit *disable_neighborhood;
    // struct arg_lit *disable_degree_1;
    // struct arg_lit *disable_degree_2;
    // struct arg_lit *disable_triangle_min;
    // struct arg_lit *disable_triangle_mid;
    // struct arg_lit *disable_v_shape_min;
    // struct arg_lit *disable_v_shape_max;
    // struct arg_lit *disable_v_shape_mid;
    // struct arg_lit *disable_domination;
    // struct arg_lit *disable_extended_domination;
    // struct arg_lit *disable_extended_domination_reverse;
    // struct arg_lit *disable_extended_twin;
    // struct arg_lit *disable_single_edge;
    // struct arg_lit *disable_extended_single_edge;
    // struct arg_lit *disable_clique_neighborhood;
    // struct arg_lit *disable_clique_neighborhood_fast;
    // struct arg_lit *disable_generalized_fold;
    // struct arg_lit *disable_critical_set;
    // struct arg_lit *disable_simplicial_vertex;
    // struct arg_lit *disable_twin;
    // struct arg_lit *disable_heavy_set;
    // struct arg_lit *disable_heavy_set3;
    // struct arg_lit *disable_cut_vertex;
    // struct arg_lit *disable_unconfined;
    // struct arg_lit *disable_funnel;
    // struct arg_lit *disable_struction_decrease;
    // struct arg_lit *disable_struction_plateau;
};

int BaseArguments::setConfig(Config &config)
{
    // Common argtable
    void *argtable[] = {
        help,
        filename,
        user_seed,
        time_limit,
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

    if (time_limit->count > 0)
    {
        config.time_limit = time_limit->dval[0];
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
        filename,
        user_seed,
        time_limit,
        verbose,
        cyclicFast,
        cyclicStrong,
        weight_source,
        print_reduction_info,
        gnn_filter,
        kernel_filename,
        reduction_config,
        // disable_neighborhood,
        // disable_degree_1,
        // disable_degree_2,
        // disable_triangle_min,
        // disable_triangle_mid,
        // disable_v_shape_min,
        // disable_v_shape_max,
        // disable_v_shape_mid,
        // disable_domination,
        // disable_extended_domination,
        // disable_extended_domination_reverse,
        // disable_extended_twin,
        // disable_single_edge,
        // disable_extended_single_edge,
        // disable_clique_neighborhood,
        // disable_clique_neighborhood_fast,
        // disable_generalized_fold,
        // disable_critical_set,
        // disable_simplicial_vertex,
        // disable_twin,
        // disable_heavy_set,
        // disable_heavy_set3,
        // disable_cut_vertex,
        // disable_unconfined,
        // disable_funnel,
        // disable_struction_decrease,
        // disable_struction_plateau,
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
    cfg.standard(config);
    // Parse the arguments
    if (cyclicFast->count > 0)
    {
        cfg.cyclicFast(config);
        config.struction_config_name = "cyclicFast";
    }
    if (cyclicStrong->count > 0)
    {
        cfg.cyclicStrong(config);
        config.struction_config_name = "cyclicStrong";
    }
    if (reduction_config->count > 0)
    {
        if (!strcmp(reduction_config->sval[0], "no_gnn_reductions"))
        {
            cfg.no_gnn_reductions(config);
            config.reduction_config_name = "no_gnn_reductions";
        }
        else if (!strcmp(reduction_config->sval[0], "all_reductions"))
        {
            cfg.all_reductions(config);
            config.reduction_config_name = "all_reductions";
        }
        else
        {
            cfg.all_reductions(config);
            config.reduction_config_name = "all_reductions";
        }
    }

    BaseArguments::parseParameters(config);

    if (gnn_filter->count > 0)
    {
        config.setGNNFilterStyle(gnn_filter->sval[0]);
    }
    if (print_reduction_info->count > 0)
    {
        config.verbose = true;
        config.print_reduction_info = true;
    }
    // if (disable_neighborhood->count > 0)
    // {
    //     config.disable_neighborhood = true;
    // }
    // if (disable_degree_1->count > 0)
    // {
    //     config.disable_fold1 = true;
    // }
    // if (disable_degree_2->count > 0)
    // {
    //     config.disable_fold2 = true;
    // }
    // if (disable_triangle_min->count > 0)
    // {
    //     config.disable_triangle_min = true;
    // }
    // if (disable_triangle_mid->count > 0)
    // {
    //     config.disable_triangle_mid = true;
    // }
    // if (disable_v_shape_min->count > 0)
    // {
    //     config.disable_v_shape_min = true;
    // }
    // if (disable_v_shape_max->count > 0)
    // {
    //     config.disable_v_shape_max = true;
    // }
    // if (disable_v_shape_mid->count > 0)
    // {
    //     config.disable_v_shape_mid = true;
    // }
    // if (disable_domination->count > 0)
    // {
    //     config.disable_domination = true;
    // }
    // if (disable_extended_domination->count > 0)
    // {
    //     config.disable_extended_domination = true;
    // }
    // if (disable_extended_domination_reverse->count > 0)
    // {
    //     config.disable_extended_domination_reverse = true;
    // }
    // if (disable_extended_twin->count > 0)
    // {
    //     config.disable_extended_twin = true;
    // }
    // if (disable_single_edge->count > 0)
    // {
    //     config.disable_single_edge = true;
    // }
    // if (disable_extended_single_edge->count > 0)
    // {
    //     config.disable_extended_single_edge = true;
    // }
    // if (disable_clique_neighborhood->count > 0)
    // {
    //     config.disable_clique_neighborhood = true;
    // }
    // if (disable_clique_neighborhood_fast->count > 0)
    // {
    //     config.disable_clique_neighborhood_fast = true;
    // }
    // if (disable_generalized_fold->count > 0)
    // {
    //     config.disable_generalized_fold = true;
    // }
    // if (disable_critical_set->count > 0)
    // {
    //     config.disable_critical_set = true;
    // }
    // if (disable_simplicial_vertex->count > 0)
    // {
    //     config.disable_simplicial_vertex = true;
    // }
    // if (disable_twin->count > 0)
    // {
    //     config.disable_twin = true;
    // }
    // if (disable_heavy_set->count > 0)
    // {
    //     config.disable_heavy_set = true;
    // }
    // if (disable_heavy_set3->count > 0)
    // {
    //     config.disable_heavy_set3 = true;
    // }
    // if (disable_cut_vertex->count > 0)
    // {
    //     config.disable_cut_vertex = true;
    // }
    // if (disable_unconfined->count > 0)
    // {
    //     config.disable_unconfined = true;
    // }
    // if (disable_funnel->count > 0)
    // {
    //     config.disable_funnel = true;
    // }
    // if (disable_struction_decrease->count > 0)
    // {
    //     config.disable_struction_decrease = true;
    // }
    // if (disable_struction_plateau->count > 0)
    // {
    //     config.disable_struction_plateau = true;
    // }

    if (kernel_filename->count > 0)
    {
        config.kernel_filename = kernel_filename->sval[0];
        config.write_kernel = true;
    }
    else
    {
        config.write_kernel = false;
    }

}
