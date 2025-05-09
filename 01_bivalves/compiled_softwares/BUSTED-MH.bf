RequireVersion ("2.5.23");


LoadFunctionLibrary("libv3/all-terms.bf"); // must be loaded before CF3x4
LoadFunctionLibrary("libv3/UtilityFunctions.bf");
LoadFunctionLibrary("libv3/IOFunctions.bf");
LoadFunctionLibrary("libv3/tasks/estimators.bf");
LoadFunctionLibrary("libv3/tasks/alignments.bf");
LoadFunctionLibrary("libv3/models/codon.bf");
LoadFunctionLibrary("libv3/tasks/trees.bf");
LoadFunctionLibrary("libv3/tasks/genetic_code.bf");
LoadFunctionLibrary("SelectionAnalyses/modules/io_functions.ibf");
LoadFunctionLibrary("SelectionAnalyses/modules/selection_lib.ibf");
LoadFunctionLibrary("libv3/models/codon/BS_REL.bf");
LoadFunctionLibrary("libv3/convenience/math.bf");
LoadFunctionLibrary("libv3/models/rate_variation.bf");
LoadFunctionLibrary("libv3/models/codon/MG_REV.bf");
LoadFunctionLibrary("libv3/models/codon/MG_REV_MH.bf");
LoadFunctionLibrary("libv3/models/codon/MG_REV_TRIP.bf");


utility.SetEnvVariable ("NORMALIZE_SEQUENCE_NAMES", TRUE);
utility.SetEnvVariable ("LF_SMOOTHING_SCALER", 0.1);


busted.analysis_description = {
                               terms.io.info :
"BUSTED-SMSH (branch-site unrestricted statistical test of episodic diversification) uses a random effects branch-site model fitted
jointly to all or a subset of tree branches in order to test for alignment-wide evidence of episodic diversifying selection and.
Assuming there is evidence of positive selection (i.e. there is an omega > 1),  BUSTED will also perform a quick evidence-ratio
style analysis to explore which individual sites may have been subject to selection. v2.0 adds support for synonymous rate variation,
and relaxes the test statistic to 0.5 (chi^2_0 + chi^2_2). Version 2.1 adds a grid search for the initial starting point.
Version 2.2 changes the grid search to LHC, and adds an initial search phase to use adaptive Nedler-Mead.
",
                               terms.io.version : "2.2",
                               terms.io.reference : "*Gene-wide identification of episodic selection*, Mol Biol Evol. 32(5):1365-71",
                               terms.io.authors : "Sergei L Kosakovsky Pond, Sadie Wisotsky",
                               terms.io.contact : "spond@temple.edu",
                               terms.io.requirements : "in-frame codon alignment and a phylogenetic tree (optionally annotated with {})"
                              };

io.DisplayAnalysisBanner (busted.analysis_description);

busted.FG = "Test";
busted.BG = "Background";
busted.SRV = "Synonymous site-to-site rates";
busted.background = "background";
busted.unconstrained = "unconstrained";
busted.constrained = "constrained";
busted.optimized_null = "optimized null";
busted.MG94 = terms.json.mg94xrev_sep_rates;
busted.MG94x2 = "MG94 with double instantaneous substitutions";
busted.MG94x3 = "MG94 with double and triple instantaneous substitutions";


busted.json.background = busted.background;
busted.json.site_logl  = "Site Log Likelihood";
busted.json.evidence_ratios  = "Evidence Ratios";
busted.json.srv_posteriors  = "Synonymous site-posteriors";
busted.rate_classes = 3;
busted.synonymous_rate_classes = 3;
busted.initial_grid.N = 250;
busted.delta.parameter = "busted.delta";
busted.psi.parameter = "busted.psi";
busted.psi_islands.parameter = "busted.psi_islands";


busted.json    = { terms.json.analysis: busted.analysis_description,
                   terms.json.input: {},
                   busted.json.background: {},
                   terms.json.fits : {},
                   terms.json.timers : {},
                   busted.json.site_logl : {},
                   busted.json.evidence_ratios: {},
                   busted.json.site_logl : {}
                  };


busted.display_orders = {terms.original_name: -1,
                         terms.json.nucleotide_gtr: 0,
                         busted.MG94: 1,
                         busted.double: 2,
                         busted.triple: 3,
                         busted.unconstrained: 4,
                         busted.constrained: 5
                        };


selection.io.startTimer (busted.json [terms.json.timers], "Overall", 0);

KeywordArgument ("code",      "Which genetic code should be used", "Universal");
    /**
        keyword, description (for inline documentation and help messages), default value
    */
KeywordArgument ("alignment", "An in-frame codon alignment in one of the formats supported by HyPhy");
    /**
        keyword, description (for inline documentation and help messages), no default value,
        meaning that it will be required
    */

KeywordArgument ("tree",      "A phylogenetic tree (optionally annotated with {})", null, "Please select a tree file for the data:");
    /** the use of null as the default argument means that the default expectation is for the
        argument to be missing, i.e. the tree is expected to be in the file
        the fourth, optional argument, can match this keyword with the dialog prompt / choice list title,
        meaning that it can only be consumed when this dialog prompt / choice list is invoked
        This allows handling some branching logic conditionals
    */

KeywordArgument ("branches",  "Branches to test", "All");
KeywordArgument ("srv", "Include synonymous rate variation in the model", "Yes");
KeywordArgument ("rates", "The number omega rate classes to include in the model [1-10, default 3]", busted.rate_classes);


namespace busted {
    LoadFunctionLibrary ("SelectionAnalyses/modules/shared-load-file.bf");
    load_file ("busted");
}


busted.do_srv = io.SelectAnOption ({"Yes" : "Allow synonymous substitution rates to vary from site to site (but not from branch to branch)",
                                    "No"  : "Synonymous substitution rates are constant across sites. This is the 'classic' behavior, i.e. the original published test"},
                                    "Synonymous rate variation"
                                    ) == "Yes";

        //setting multihit rate as default for now
        /*
busted.multihit = io.SelectAnOption ({"1" : "Use standard MG94 model.",
                                      "2"  : "Use MG94 with instantaneous double hits allowed.",
                                      "3" : "Use MG94 with instantaneous double and triple hits allowed."},
                                      "Multiple hit rate"
                                      ) == "3";
                                      */


busted.rate_classes = io.PromptUser ("The number omega rate classes to include in the model", busted.rate_classes, 1, 10, TRUE);

if (busted.do_srv) {
    KeywordArgument ("syn-rates", "The number alpha rate classes to include in the model [1-10, default 3]", busted.synonymous_rate_classes);
    busted.synonymous_rate_classes = io.PromptUser ("The number omega rate classes to include in the model", busted.synonymous_rate_classes, 1, 10, TRUE);
}

KeywordArgument ("grid-size", "The number of points in the initial distributional guess for likelihood fitting", 250);
busted.initial_grid.N = io.PromptUser ("The number of points in the initial distributional guess for likelihood fitting", 250, 1, 10000, TRUE);
KeywordArgument ("starting-points", "The number of initial random guesses to seed rate values optimization", 1);
busted.N.initial_guesses = io.PromptUser ("The number of initial random guesses to 'seed' rate values optimization", 1, 1, busted.initial_grid.N$10, TRUE);

KeywordArgument ("output", "Write the resulting JSON to this file (default is to save to the same path as the alignment file + 'BUSTED.json')", busted.codon_data_info [terms.json.json]);
busted.codon_data_info [terms.json.json] = io.PromptUserForFilePath ("Save the resulting JSON file to");

io.ReportProgressMessageMD('BUSTED',  'selector', 'Branches to test for selection in the BUSTED analysis');

utility.ForEachPair (busted.selected_branches, "_partition_", "_selection_",
    "_selection_ = utility.Filter (_selection_, '_value_', '_value_ == terms.tree_attributes.test');
     io.ReportProgressMessageMD('BUSTED',  'selector', '* Selected ' + Abs(_selection_) + ' branches to test in the BUSTED analysis: \\\`' + Join (', ',utility.Keys(_selection_)) + '\\\`')");

// check to see if there are any background branches
busted.has_background = FALSE;

utility.ForEachPair (busted.selected_branches, "_partition_", "_selection_",
    "_selection_ = utility.Filter (_selection_, '_value_', '_value_ != terms.tree_attributes.test');
     if (utility.Array1D (_selection_)) { busted.has_background = TRUE;} ");
busted.json[busted.json.background] =  busted.has_background;

selection.io.startTimer (busted.json [terms.json.timers], "Preliminary model fitting", 1);

namespace busted {
    doGTR ("busted");
}


estimators.fixSubsetOfEstimates(busted.gtr_results, busted.gtr_results[terms.global]);


io.ReportProgressMessageMD ("BUSTED", "codon-refit", "Improving branch lengths, nucleotide substitution biases, and global dN/dS ratios under a full codon model");
///this is where busted fits MG94, so I think we want to add similar stuff for double and MH

busted.final_partitioned_mg_results = estimators.FitMGREV (busted.filter_names, busted.trees, busted.codon_data_info [terms.code], {
    terms.run_options.model_type: terms.local,
    terms.run_options.partitioned_omega: busted.selected_branches,
}, busted.partitioned_mg_results);


io.ReportProgressMessageMD("BUSTED", "codon-refit", "* " + selection.io.report_fit (busted.final_partitioned_mg_results, 0, busted.codon_data_info[terms.data.sample_size]));
busted.global_dnds = selection.io.extract_global_MLE_re (busted.final_partitioned_mg_results, "^" + terms.parameters.omega_ratio);
utility.ForEach (busted.global_dnds, "_value_", 'io.ReportProgressMessageMD ("BUSTED", "codon-refit", "* " + _value_[terms.description] + " = " + Format (_value_[terms.fit.MLE],8,4));');


utility.Extend (busted.final_partitioned_mg_results[terms.global],
                {
                    terms.parameters.multiple_hit_rate : { utility.getGlobalValue ("terms.fit.MLE") : 0.05, terms.fix : FALSE}

                });


busted.two_hit_results = busted.run_model_fit (busted.MG94x2, "models.codon.MG_REV_MH.ModelDescription", busted.final_partitioned_mg_results, "busted.init_delta");


utility.Extend (busted.two_hit_results[terms.global],
                {
                    terms.parameters.triple_hit_rate : { utility.getGlobalValue ("terms.fit.MLE") : 0.05, terms.fix : FALSE},
                    terms.parameters.triple_hit_rate_syn : { utility.getGlobalValue ("terms.fit.MLE") : 0.05, terms.fix : FALSE}

                });

busted.three_hit_results = busted.run_model_fit (busted.MG94x3, "models.codon.MG_REV_TRIP.ModelDescription", busted.two_hit_results, "busted.init_triple");



selection.io.stopTimer (busted.json [terms.json.timers], "Preliminary model fitting");

//Store MG94 to JSON
selection.io.json_store_lf_withEFV (busted.json,
                            busted.MG94,
                            busted.three_hit_results[terms.fit.log_likelihood],
                            busted.three_hit_results[terms.parameters],
                            busted.sample_size,
                            utility.ArrayToDict (utility.Map (busted.global_dnds, "_value_", "{'key': _value_[terms.description], 'value' : Eval({{_value_ [terms.fit.MLE],1}})}")),
                            (busted.three_hit_results[terms.efv_estimate])["VALUEINDEXORDER"][0],
                            busted.display_orders[busted.MG94]);

utility.ForEachPair (busted.filter_specification, "_key_", "_value_",
    'selection.io.json_store_branch_attribute(busted.json, busted.MG94, terms.branch_length, busted.display_orders[busted.MG94],
                                             _key_,
                                             selection.io.extract_branch_info((busted.three_hit_results[terms.branch_length])[_key_], "selection.io.branch.length"));');



busted.model_generator = "models.codon.BS_REL.ModelDescription.MH";

lfunction models.codon.BS_REL.ModelDescription.MH (type, code, components) {
    def = models.codon.BS_REL.ModelDescription (type, code, components);
    def [utility.getGlobalValue("terms.model.defineQ")] = "models.codon.BS_REL._DefineQ.MH";
    return def;
}

if (busted.do_srv) {
    lfunction busted.model.with.GDD (type, code, rates) {
        def = models.codon.BS_REL.ModelDescription.MH (type, code, rates);
        def [utility.getGlobalValue("terms.model.rate_variation")] = rate_variation.types.GDD.factory ({utility.getGlobalValue("terms.rate_variation.bins") : utility.getGlobalValue("busted.synonymous_rate_classes")});
        return def;
    }
    busted.model_generator = "busted.model.with.GDD";
}

busted.test.bsrel_model =  model.generic.DefineMixtureModel(busted.model_generator,
        "busted.test", {
            "0": parameters.Quote(terms.global),
            "1": busted.codon_data_info[terms.code],
            "2": parameters.Quote (busted.rate_classes) // the number of rate classes
        },
        busted.filter_names,
        None);

busted.background.bsrel_model =  model.generic.DefineMixtureModel(busted.model_generator,
        "busted.background", {
            "0": parameters.Quote(terms.global),
            "1": busted.codon_data_info[terms.code],
            "2": parameters.Quote (busted.rate_classes) // the number of rate classes
        },
        busted.filter_names,
        None);


models.BindGlobalParameters ({"0" : busted.test.bsrel_model, "1" : busted.background.bsrel_model}, terms.nucleotideRate("[ACGT]","[ACGT]"));


if (busted.do_srv) {
    models.BindGlobalParameters ({"0" : busted.test.bsrel_model, "1" : busted.background.bsrel_model}, "GDD rate category");
    models.BindGlobalParameters ({"0" : busted.test.bsrel_model, "1" : busted.background.bsrel_model}, utility.getGlobalValue("terms.mixture.mixture_aux_weight") + " for GDD category ");
}

busted.distribution = models.codon.BS_REL.ExtractMixtureDistribution(busted.test.bsrel_model);

// set up parameter constraints

for (busted.i = 1; busted.i < busted.rate_classes; busted.i += 1) {
    parameters.SetRange (model.generic.GetGlobalParameter (busted.test.bsrel_model , terms.AddCategory (terms.parameters.omega_ratio,busted.i)), terms.range01);
    parameters.SetRange (model.generic.GetGlobalParameter (busted.background.bsrel_model , terms.AddCategory (terms.parameters.omega_ratio,busted.i)), terms.range01);
}


parameters.SetRange (model.generic.GetGlobalParameter (busted.test.bsrel_model , terms.AddCategory (terms.parameters.omega_ratio,busted.rate_classes)), terms.range_gte1);

/* create an initial grid to initialize the optimization */
busted.initial_grid = {};

/* if populated, use this as a baseline to generate the distributions from */
busted.initial_grid_presets = {"0" : 0.1};
busted.initial_ranges = {};

busted.init_grid_setup (busted.distribution);

/** setup parameter optimization groups */

PARAMETER_GROUPING = {};
PARAMETER_GROUPING + busted.distribution["rates"];
PARAMETER_GROUPING + busted.distribution["weights"];

busted.initial_ranges [((busted.test.bsrel_model[terms.parameters])[terms.global])[terms.parameters.multiple_hit_rate]] = {
                    terms.lower_bound : 0,
                    terms.upper_bound : 1
                };
busted.initial_ranges [((busted.test.bsrel_model[terms.parameters])[terms.global])[terms.parameters.triple_hit_rate]] = {
                    terms.lower_bound : 0,
                    terms.upper_bound : 1
                };
busted.initial_ranges [((busted.test.bsrel_model[terms.parameters])[terms.global])[terms.parameters.triple_hit_rate_syn]] = {
                    terms.lower_bound : 0,
                    terms.upper_bound : 1
                };

if (busted.has_background) {
    busted.model_object_map = { "busted.background" : busted.background.bsrel_model,
                                "busted.test" :       busted.test.bsrel_model };
    busted.background_distribution = models.codon.BS_REL.ExtractMixtureDistribution(busted.background.bsrel_model);
    busted.init_grid_setup (busted.background_distribution);

    busted.initial_ranges [((busted.background.bsrel_model[terms.parameters])[terms.global])[terms.parameters.multiple_hit_rate]] = {
                        terms.lower_bound : 0,
                        terms.upper_bound : 1
                    };
    busted.initial_ranges [((busted.background.bsrel_model[terms.parameters])[terms.global])[terms.parameters.triple_hit_rate]] = {
                        terms.lower_bound : 0,
                        terms.upper_bound : 1
                    };
    busted.initial_ranges [((busted.background.bsrel_model[terms.parameters])[terms.global])[terms.parameters.triple_hit_rate_syn]] = {
                        terms.lower_bound : 0,
                        terms.upper_bound : 1
                    };
    PARAMETER_GROUPING = {};
    PARAMETER_GROUPING + busted.background_distribution["rates"];
    PARAMETER_GROUPING + busted.background_distribution["weights"];

} else {
    busted.model_object_map = { "busted.test" :       busted.test.bsrel_model };
}

if (busted.do_srv)  {
    busted.srv_rate_regex  = "GDD rate category [0-9]+";
    busted.srv_weight_regex = "Mixture auxiliary weight for GDD category [0-9]+";
    busted.srv_distribution = regexp.PartitionByRegularExpressions(utility.Keys ((busted.test.bsrel_model[terms.parameters])[terms.global]), {"0" : busted.srv_rate_regex, "1" : busted.srv_weight_regex});


    busted.srv_distribution = {
        'rates' : utility.UniqueValues (utility.Map (busted.srv_distribution [busted.srv_rate_regex ]  , "_value_", '((busted.test.bsrel_model[terms.parameters])[terms.global])[_value_]')),
        'weights' : utility.UniqueValues (utility.Map (busted.srv_distribution [busted.srv_weight_regex ]  , "_value_", '((busted.test.bsrel_model[terms.parameters])[terms.global])[_value_]'))
    };

    PARAMETER_GROUPING + busted.srv_distribution["rates"];
    PARAMETER_GROUPING + busted.srv_distribution["weights"];

    busted.init_grid_setup (busted.srv_distribution);

}

busted.initial.test_mean = ((busted.global_dnds)["0"])["MLE"];
//busted.initial.test_mean    = ((selection.io.extract_global_MLE_re (busted.three_hit_results, "^" + terms.parameters.omega_ratio + ".+test.+"))["0"])[terms.fit.MLE];
busted.initial_grid         = estimators.LHC (busted.initial_ranges,busted.initial_grid.N);

busted.initial_grid = utility.Map (busted.initial_grid, "_v_",
    'busted._renormalize (_v_, "busted.distribution", busted.initial.test_mean)'
);

if (busted.has_background) { //GDD rate category
    busted.initial.background_mean    = ((selection.io.extract_global_MLE_re (busted.three_hit_results, "^" + terms.parameters.omega_ratio + ".+background.+"))["0"])[terms.fit.MLE];
    busted.initial_grid = utility.Map (busted.initial_grid, "_v_",
        'busted._renormalize (_v_, "busted.background_distribution", busted.initial.background_mean)'
    );
}


busted.model_map = {};

for (busted.partition_index = 0; busted.partition_index < busted.partition_count; busted.partition_index += 1) {
    selection.io.json_store_branch_attribute(busted.json, terms.original_name, terms.json.node_label, 0,
                                             busted.partition_index,
                                             busted.name_mapping);

    busted.model_map + { "busted.test" : utility.Filter (busted.selected_branches[busted.partition_index], '_value_', '_value_ == terms.tree_attributes.test'),
    		       	 	            "busted.background" : utility.Filter (busted.selected_branches[busted.partition_index], '_value_', '_value_ != terms.tree_attributes.test')};
}

utility.SetEnvVariable ("ASSUME_REVERSIBLE_MODELS", TRUE);

selection.io.startTimer (busted.json [terms.json.timers], "Unconstrained BUSTED model fitting", 2);

io.ReportProgressMessageMD ("BUSTED", "main", "Performing the full (dN/dS > 1 allowed) branch-site model fit");

/**
    perform the initial fit using a constrained branch length model and
    a low precision Nedler-Mead algorithm pass to find good initial values
    for the rate distribution parameters
*/


busted.nm.precision = -0.00025*busted.three_hit_results[terms.fit.log_likelihood];

//VERBOSITY_LEVEL = 10;

parameters.DeclareGlobalWithRanges ("busted.bl.scaler", 1, 0, 1000);
busted.grid_search.results =  estimators.FitLF (busted.filter_names, busted.trees, busted.model_map, busted.three_hit_results, busted.model_object_map, {
    "retain-lf-object": TRUE,
    terms.run_options.proportional_branch_length_scaler :
                                            {"0" : "busted.bl.scaler"},

    terms.run_options.optimization_settings :
        {
            "OPTIMIZATION_METHOD" : "nedler-mead",
            "MAXIMUM_OPTIMIZATION_ITERATIONS" : 500,
            "OPTIMIZATION_PRECISION" : busted.nm.precision
        } ,

    terms.search_grid     : busted.initial_grid,
    terms.search_restarts : busted.N.initial_guesses
});

busted.full_model =  estimators.FitLF (busted.filter_names, busted.trees, busted.model_map, busted.grid_search.results, busted.model_object_map, {
        "retain-lf-object": TRUE,
        terms.run_options.optimization_settings :
            {
                "OPTIMIZATION_METHOD" : "hybrid"
            }

    });


KeywordArgument ("save-fit", "Save BUSTED model fit to this file (default is not to save)", "/dev/null");
io.SpoolLFToPath(busted.full_model[terms.likelihood_function], io.PromptUserForFilePath ("Save BUSTED model fit to this file ['/dev/null' to skip]"));

io.ReportProgressMessageMD("BUSTED", "main", "* " + selection.io.report_fit (busted.full_model, 9, busted.codon_data_info[terms.data.sample_size]));
busted.global_dnds = selection.io.extract_global_MLE_re (busted.full_model, "(" + terms.parameters.multiple_hit_rate  + "|" + terms.parameters.triple_hit_rate + "|" + terms.parameters.triple_hit_rate + ")");
utility.ForEach (busted.global_dnds, "_value_", 'io.ReportProgressMessageMD ("BUSTED", "main", "* " + _value_[terms.description] + " = " + Format (_value_[terms.fit.MLE],8,4));');

io.ReportProgressMessageMD("BUSTED", "main", "* For *test* branches, the following rate distribution for branch-site combinations was inferred");

selection.io.stopTimer (busted.json [terms.json.timers], "Unconstrained BUSTED model fitting");

busted.inferred_test_distribution = parameters.GetStickBreakingDistribution (busted.distribution) % 0;
selection.io.report_dnds (busted.inferred_test_distribution);


busted.distribution_for_json = {busted.FG : utility.Map (utility.Range (busted.rate_classes, 0, 1),
                                                         "_index_",
                                                         "{terms.json.omega_ratio : busted.inferred_test_distribution [_index_][0],
                                                           terms.json.proportion : busted.inferred_test_distribution [_index_][1]}")
                                };


if (busted.has_background) {
    io.ReportProgressMessageMD("BUSTED", "main", "* For *background* branches, the following rate distribution for branch-site combinations was inferred");
    busted.inferred_background_distribution = parameters.GetStickBreakingDistribution (busted.background_distribution) % 0;
    selection.io.report_dnds (busted.inferred_background_distribution);
    busted.distribution_for_json [busted.BG] = utility.Map (utility.Range (busted.rate_classes, 0, 1),
                                                         "_index_",
                                                         "{terms.json.omega_ratio : busted.inferred_background_distribution [_index_][0],
                                                           terms.json.proportion : busted.inferred_background_distribution [_index_][1]}");
}

if (busted.do_srv) {
    busted.srv_info = Transpose((rate_variation.extract_category_information(busted.test.bsrel_model))["VALUEINDEXORDER"][0])%0;
    io.ReportProgressMessageMD("BUSTED", "main", "* The following rate distribution for site-to-site **synonymous** rate variation was inferred");
    selection.io.report_distribution (busted.srv_info);

    busted.distribution_for_json [busted.SRV] = (utility.Map (utility.Range (busted.synonymous_rate_classes, 0, 1),
                                                             "_index_",
                                                             "{terms.json.rate :busted.srv_info [_index_][0],
                                                               terms.json.proprtion : busted.srv_info [_index_][1]}"));

    ConstructCategoryMatrix (busted.cmx, ^(busted.full_model[terms.likelihood_function]));
    ConstructCategoryMatrix (busted.cmx_weights, ^(busted.full_model[terms.likelihood_function]), WEIGHTS);
    busted.cmx_weighted         = busted.cmx_weights $ busted.cmx;
    busted.column_weights       = {1, Rows (busted.cmx_weights)}["1"] * busted.cmx_weighted;
    busted.column_weights       = busted.column_weights["1/_MATRIX_ELEMENT_VALUE_"];
    (busted.json [busted.json.srv_posteriors]) =  busted.cmx_weighted $ busted.column_weights;

}



selection.io.json_store_lf (busted.json,
                            "Unconstrained model",
                            busted.full_model[terms.fit.log_likelihood],
                            busted.full_model[terms.parameters] + 9 , // +9 comes from CF3x4
                            busted.codon_data_info[terms.data.sample_size],
                            busted.distribution_for_json,
                            busted.display_orders[busted.unconstrained]);

utility.ForEach (busted.global_dnds, "_value_", '((busted.json[terms.json.fits])["Unconstrained model"])[_value_[terms.description]]  = _value_[terms.fit.MLE]');



busted.run_test = busted.inferred_test_distribution [busted.rate_classes-1][0] > 1 && busted.inferred_test_distribution [busted.rate_classes-1][1] > 0;

utility.ForEachPair (busted.filter_specification, "_key_", "_value_",
    'selection.io.json_store_branch_attribute(busted.json, busted.unconstrained, terms.branch_length, 0,
                                             _key_,
                                             selection.io.extract_branch_info((busted.full_model[terms.branch_length])[_key_], "selection.io.branch.length"));');


(busted.json [busted.json.site_logl])[busted.unconstrained] = busted.ComputeSiteLikelihoods (busted.full_model[terms.likelihood_function]);

if (!busted.run_test) {
    io.ReportProgressMessageMD ("BUSTED", "Results", "No evidence for episodic diversifying positive selection under the unconstrained model, skipping constrained model fitting");
    busted.json [terms.json.test_results] = busted.ComputeLRT (0, 0);
} else {
    selection.io.startTimer (busted.json [terms.json.timers], "Constrained BUSTED model fitting", 3);


    io.ReportProgressMessageMD ("BUSTED", "test", "Performing the constrained (dN/dS > 1 not allowed) model fit");
    parameters.SetConstraint (model.generic.GetGlobalParameter (busted.test.bsrel_model , terms.AddCategory (terms.parameters.omega_ratio,busted.rate_classes)), terms.parameters.one, terms.global);
    (busted.json [busted.json.site_logl])[busted.constrained] = busted.ComputeSiteLikelihoods (busted.full_model[terms.likelihood_function]);
    busted.null_results = estimators.FitExistingLF (busted.full_model[terms.likelihood_function], busted.model_object_map);
    (busted.json [busted.json.site_logl])[busted.optimized_null] = busted.ComputeSiteLikelihoods (busted.full_model[terms.likelihood_function]);
    io.ReportProgressMessageMD ("BUSTED", "test", "* " + selection.io.report_fit (busted.null_results, 9, busted.codon_data_info[terms.data.sample_size]));
    busted.LRT = busted.ComputeLRT (busted.full_model[terms.fit.log_likelihood], busted.null_results[terms.fit.log_likelihood]);
    busted.json [terms.json.test_results] = busted.LRT;


    selection.io.stopTimer (busted.json [terms.json.timers], "Constrained BUSTED model fitting");

    utility.ForEachPair (busted.filter_specification, "_key_", "_value_",
        'selection.io.json_store_branch_attribute(busted.json, busted.constrained, terms.branch_length, 1,
                                                 _key_,
                                                 selection.io.extract_branch_info((busted.null_results[terms.branch_length])[_key_], "selection.io.branch.length"));');

    busted.global_dnds = selection.io.extract_global_MLE_re (busted.null_results, "(" + terms.parameters.multiple_hit_rate  + "|" + terms.parameters.triple_hit_rate + "|" + terms.parameters.triple_hit_rate + ")");
    utility.ForEac (busted.global_dnds, "_value_", 'io.ReportProgressMessageMD ("BUSTED", "test", "* " + _value_[terms.description] + " = " + Format (_value_[terms.fit.MLE],8,4));');

    io.ReportProgressMessageMD("BUSTED", "test", "* For *test* branches under the null (no dN/dS > 1 model), the following rate distribution for branch-site combinations was inferred");

    busted.inferred_test_distribution = parameters.GetStickBreakingDistribution (busted.distribution) % 0;
    selection.io.report_dnds (parameters.GetStickBreakingDistribution (busted.distribution) % 0);

    busted.distribution_for_json = {busted.FG : utility.Map (utility.Range (busted.rate_classes, 0, 1),
                                                         "_index_",
                                                         "{terms.json.omega_ratio : busted.inferred_test_distribution [_index_][0],
                                                           terms.json.proportion : busted.inferred_test_distribution [_index_][1]}")};

    if (busted.has_background) {
        busted.inferred_background_distribution = parameters.GetStickBreakingDistribution (busted.background_distribution) % 0;
        //selection.io.report_dnds (busted.inferred_background_distribution);
        busted.distribution_for_json [busted.BG] = utility.Map (utility.Range (busted.rate_classes, 0, 1),
                                                             "_index_",
                                                             "{terms.json.omega_ratio : busted.inferred_background_distribution [_index_][0],
                                                               terms.json.proportion : busted.inferred_background_distribution [_index_][1]}");
    }

    if (busted.do_srv) {
        busted.srv_info = Transpose((rate_variation.extract_category_information(busted.test.bsrel_model))["VALUEINDEXORDER"][0])%0;
        io.ReportProgressMessageMD("BUSTED", "main", "* The following rate distribution for site-to-site **synonymous** rate variation was inferred");
        selection.io.report_distribution (busted.srv_info);

         busted.distribution_for_json [busted.SRV] = (utility.Map (utility.Range (busted.synonymous_rate_classes, 0, 1),
                                                                 "_index_",
                                                                 "{terms.json.rate :busted.srv_info [_index_][0],
                                                                   terms.json.proportion : busted.srv_info [_index_][1]}"));

    }

    selection.io.json_store_lf (busted.json,
                            "Constrained model",
                            busted.null_results[terms.fit.log_likelihood],
                            busted.null_results[terms.parameters] + 9 , // +9 comes from CF3x4
                            busted.codon_data_info[terms.data.sample_size],
                            busted.distribution_for_json,
                            busted.display_orders[busted.constrained]);

    utility.ForEach (busted.global_dnds, "_value_", '((busted.json[terms.json.fits])["Constrained model"])[_value_[terms.description]]  = _value_[terms.fit.MLE]');

    (busted.json [busted.json.evidence_ratios])[busted.constrained] = busted.EvidenceRatios ( (busted.json [busted.json.site_logl])[busted.unconstrained],  (busted.json [busted.json.site_logl])[busted.constrained]);
    (busted.json [busted.json.evidence_ratios ])[busted.optimized_null] = busted.EvidenceRatios ( (busted.json [busted.json.site_logl])[busted.unconstrained],  (busted.json [busted.json.site_logl])[busted.optimized_null]);
}


console.log ("----\n## Branch-site unrestricted statistical test of episodic diversification [BUSTED]");
console.log ( "Likelihood ratio test for episodic diversifying positive selection, **p = " + Format ((busted.json [terms.json.test_results])[terms.p_value], 8, 4) + "**.");

selection.io.stopTimer (busted.json [terms.json.timers], "Overall");

io.SpoolJSON (busted.json, busted.codon_data_info [terms.json.json]);

return busted.json;

/// HELPERS

//-----------------------------------------------------------------------------
lfunction busted.ComputeSiteLikelihoods (id) {
    ConstructCategoryMatrix (sl, ^id, SITE_LOG_LIKELIHOODS);
   return sl;
}
//------------------------------------------------------------------------------

function busted.ComputeLRT (ha, h0) {
    return {terms.LRT : 2*(ha-h0),
            terms.p_value : 0.5*(1-CChi2 (2*(ha-h0),2))};
}

//------------------------------------------------------------------------------
lfunction busted.EvidenceRatios (ha, h0) {
    return ha["Exp(_MATRIX_ELEMENT_VALUE_-h0[_MATRIX_ELEMENT_COLUMN_])"];
}
//------------------------------------------------------------------------------

lfunction busted.DistributionGuess (mean) {
    /*guess = {{Random (0,0.25),Random (0.5,0.7)}
             {Random (0.25,0.6),Random (0.1,0.2)}
             {Random (2,20),Random (0.01, 0.1)}};
    */

    guess = {{0,0.7}{0.25,0.2}{10,0.1}};

    norm = + guess[-1][1];
    guess_mean = 1/(+(guess [-1][0] $ guess [-1][1]))/norm;
    return guess["_MATRIX_ELEMENT_VALUE_*(guess_mean*(_MATRIX_ELEMENT_COLUMN_==0)+(_MATRIX_ELEMENT_COLUMN_==1)*(1/norm))"];
}

//------------------------------------------------------------------------------

// renormalize the grid to have the same mean as the initial omega value


lfunction busted._renormalize (v, distro, mean) {
    parameters.SetValues (v);
    m = parameters.GetStickBreakingDistribution (^distro);
    d = Rows (m);
    m = +(m[-1][0] $ m[-1][1]); // current mean
    for (i = 0; i < d; i+=1) {
        (v[((^"busted.distribution")["rates"])[i]])[^"terms.fit.MLE"] = (v[((^"busted.distribution")["rates"])[i]])[^"terms.fit.MLE"] / m * mean;
    }
    return v;

}

//------------------------------------------------------------------------------

function busted.init_grid_setup (omega_distro) {
    utility.ForEachPair (omega_distro[terms.parameters.rates], "_index_", "_name_",
        '
            if (_index_[0] < busted.rate_classes - 1) { // not the last rate
                busted.initial_grid  [_name_] = {
                    {
                        0.01, 0.1, 0.25, 0.75
                    }
                }["_MATRIX_ELEMENT_VALUE_^(busted.rate_classes-_index_[0]-1)"];
                busted.initial_grid_presets [_name_] = 0;
                busted.initial_ranges [_name_] = {
                    terms.lower_bound : 0,
                    terms.upper_bound : 1
                };
            }  else {
                busted.initial_grid  [_name_] = {
                    {
                        1, 1.5, 2, 4, 10
                    }
                };
                busted.initial_ranges [_name_] = {
                    terms.lower_bound : 1,
                    terms.upper_bound : 10
                };
                busted.initial_grid_presets [_name_] = 2;
            }
        '
    );


    utility.ForEachPair (omega_distro[terms.parameters.weights], "_index_", "_name_",
        '
            busted.initial_grid  [_name_] = {
                {
                    0.2, 0.5, 0.7, 0.8, 0.9
                }
            };
            busted.initial_grid_presets [_name_] = 3;
            busted.initial_ranges [_name_] = {
                terms.lower_bound : 0,
                terms.upper_bound : 1
            };
        '
    );

}

function busted.init_triple (lf_id, components, data_filter, tree, model_map, initial_values, model_objects) {
    // define the shared global delta (one for all branches, could extend to be partition specific)
    // should be a single model here
    busted.init_delta (lf_id, components, data_filter, tree, model_map, initial_values, model_objects);
    parameters.DeclareGlobalWithRanges (busted.psi.parameter, 0.05, 0, 1000);
    model.generic.AddGlobal (model_objects[utility.Keys(model_objects)[0]],busted.psi.parameter,terms.parameters.triple_hit_rate);
    model.generic.AddGlobal (model_objects[utility.Keys(model_objects)[0]],busted.psi_islands.parameter,terms.parameters.triple_hit_rate_syn);
    estimators.TraverseLocalParameters (lf_id, model_objects, "busted.set.psi");
    estimators.TraerseLocalParameters (lf_id, model_objects, "busted.set.psi_islands");
    return 0;
}

function busted.init_delta  (lf_id, components, data_filter, tree, model_map, initial_values, model_objects) {
    // define the shared global delta (one for all branches, could extend to be partition specific)
    // should be a single model here
    parameters.DeclareGlobalWithRanges (busted.delta.parameter, 0.05, 0, 1000);
    model.generic.AddGlobal (model_objects[utility.Keys(model_objects)[0]],busted.delta.parameter,terms.parameters.multiple_hit_rate);
    estimators.TraverseLocalParameters (lf_id, model_objects, "busted.set.delta");
    return 0;
}

//------------------------------------------------------------------------------

lfunction busted.set.delta (tree_name, node_name, model_description, unused) {
    if (utility.Has (model_description [utility.getGlobalValue ("terms.local")], utility.getGlobalValue ("terms.parameters.multiple_hit_rate"), "String")) {
        k = (model_description [utility.getGlobalValue ("terms.local")])[utility.getGlobalValue ("terms.parameters.multiple_hit_rate")];
        //console.log (k);
        parameters.SetConstraint (tree_name + "." + node_name + "." + k, utility.getGlobalValue ("busted.delta.parameter"), "");
        return tree_name + "." + node_name + "." + k;
    }
    return "";
}

lfunction busted.set.psi (tree_name, node_name, model_description, unused) {
    if (utility.Has (model_description [utility.getGlobalValue ("terms.local")], utility.getGlobalValue ("terms.parameters.triple_hit_rate"), "String")) {
        k = (model_description [utility.getGlobalValue ("terms.local")])[utility.getGlobalValue ("terms.parameters.triple_hit_rate")];
        //console.log (k);
        parameters.SetConstraint (tree_name + "." + node_name + "." + k, utility.getGlobalValue ("busted.psi.parameter"), "");
        return tree_name + "." + node_name + "." + k;
    }
    return "";
}

lfunction busted.set.psi_islands (tree_name, node_name, model_description, unused) {
    if (utility.Has (model_description [utility.getGlobalValue ("terms.local")], utility.getGlobalValue ("terms.parameters.triple_hit_rate_syn"), "String")) {
        k = (model_description [utility.getGlobalValue ("terms.local")])[utility.getGlobalValue ("terms.parameters.triple_hit_rate_syn")];
        //console.log (k);
        parameters.SetConstraint (tree_name + "." + node_name + "." + k, utility.getGlobalValue ("busted.psi_islands.parameter"), "");
        return tree_name + "." + node_name + "." + k;
    }
    return "";
}

//------------------------------------------------------------------------------


function busted.run_model_fit (model_name, model_generator, initial_values, rate_initializer) {
    io.ReportProgressMessageMD ("busted", model_name, "Fitting `model_name`");
    selection.io.startTimer (busted.json [terms.json.timers], model_name, busted.display_order [model_name]);


    busted.results =  estimators.FitCodonModel (busted.filter_names, busted.trees, model_generator, busted.codon_data_info [utility.getGlobalValue("terms.code")],
     {
        utility.getGlobalValue("terms.run_options.retain_lf_object"): TRUE,
        utility.getGlobalValue("terms.run_options.model_type"): utility.getGlobalValue("terms.local"),
        //utility.getGlobalValue("terms.run_options.proportional_branch_length_scaler"): scaler_variables,
        utility.getGlobalValue("terms.run_options.partitioned_omega"): busted.selected_branches,
        utility.getGlobalValue("terms.run_options.apply_user_constraints"): rate_initializer

    },initial_values);


    //ConstructCategoryMatrix (busted.run_model_fit.sl, ^(busted.results[terms.likelihood_function]), SITE_LOG_LIKELIHOODS);
    //(busted.json [busted.json.site_logl])[model_name] = busted.run_model_fit.sl;

    io.ReportProgressMessageMD("busted", model_name, "* " + selection.io.report_fit (busted.results, 0, busted.codon_data_info[terms.data.sample_size]));
    busted.global_dnds = selection.io.extract_global_MLE_re (busted.results, "^(" + terms.parameters.omega_ratio + "|" + terms.parameters.multiple_hit_rate  + "|" + terms.parameters.triple_hit_rate + ")");

    utility.ForEach (busted.global_dnds, "_value_", 'io.ReportProgressMessageMD ("busted", model_name, "* " + _value_[terms.description] + " = " + Format (_value_[terms.fit.MLE],8,4));');

    selection.io.json_store_lf (busted.json,
                                model_name,
                                busted.results[terms.fit.log_likelihood],
                                busted.results[terms.parameters],
                                busted.sample_size,
                                utility.Map (busted.results[terms.global], "_value_", '_value_ [terms.fit.MLE]'),
                                busted.display_orders[model_name]);


    utility.ForEachPair (busted.filter_specification, "_key_", "_value_",
        'selection.io.json_store_branch_attribute(busted.json, model_name, terms.branch_length, busted.display_orders[model_name],
                                                 _key_,
                                                 selection.io.extract_branch_info((busted.results[terms.branch_length])[_key_], "selection.io.branch.length"));');

    selection.io.stopTimer (busted.json [terms.json.timers], model_name);
    return busted.results;
}

/**
 * @name models.codon.BS_REL.BS_REL._DefineQ
 * @param {Dict} mg_rev
 * @param {String} namespace
 * @returns {Dict} updated model
 */

lfunction models.codon.BS_REL._DefineQ.MH (bs_rel, namespace) {


    rate_matrices = {};

    bs_rel [utility.getGlobalValue("terms.model.q_ij")] = &rate_generator;
    bs_rel [utility.getGlobalValue("terms.mixture.mixture_components")] = {};

    _aux = parameters.GenerateSequentialNames (namespace + ".bsrel_mixture_aux", bs_rel[utility.getGlobalValue("terms.model.components")] - 1, "_");
    _wts = parameters.helper.stick_breaking (_aux, None);


    for (component = 1; component <= bs_rel[utility.getGlobalValue("terms.model.components")]; component += 1) {
       key = "component_" + component;
       ExecuteCommands ("
        function rate_generator (fromChar, toChar, namespace, model_type, model) {
               return models.codon.MG_REV_TRIP._GenerateRate_generic (fromChar, toChar, namespace, model_type, model[utility.getGlobalValue('terms.translation_table')],
                'alpha', utility.getGlobalValue('terms.parameters.synonymous_rate'),
                'beta_`component`', terms.AddCategory (utility.getGlobalValue('terms.parameters.nonsynonymous_rate'), component),
                'omega`component`', terms.AddCategory (utility.getGlobalValue('terms.parameters.omega_ratio'), component),
                'delta', utility.getGlobalValue('terms.parameters.multiple_hit_rate'),
                'psi', utility.getGlobalValue('terms.parameters.triple_hit_rate'),
                'psi_islands', utility.getGlobalValue('terms.parameters.triple_hit_rate_syn')
               );
            }"
       );

       if ( component < bs_rel[utility.getGlobalValue("terms.model.components")]) {
            model.generic.AddGlobal ( bs_rel, _aux[component-1], terms.AddCategory (utility.getGlobalValue("terms.mixture.mixture_aux_weight"), component ));
            parameters.DeclareGlobalWithRanges (_aux[component-1], 0.5, 0, 1);
       }
       models.codon.generic.DefineQMatrix(bs_rel, namespace);
       rate_matrices [key] = bs_rel[utility.getGlobalValue("terms.model.rate_matrix")];
       (bs_rel [^'terms.mixture.mixture_components'])[key] = _wts [component-1];
    }

    bs_rel[utility.getGlobalValue("terms.model.rate_matrix")] = rate_matrices;
    parmeters.SetConstraint(((bs_rel[utility.getGlobalValue("terms.parameters")])[utility.getGlobalValue("terms.global")])[terms.nucleotideRate("A", "G")], "1", "");

    return bs_rel;
}
