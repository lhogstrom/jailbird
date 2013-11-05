import test_modules.pcla_svm_classifier as psc

wkdir = '/xchip/cogs/projects/pharm_class/svm_pcla_classifier_NOV1'

reload(psc)
pso = psc.svm_pcla(out=wkdir)
pso.set_classes()
pso.classification_by_cell()
pso.classification_across_cell(loo_type='by_cp',max_signatures_per_cp=3,groups_to_model=pso.test_groups)
rnkptMed_file = '/xchip/cogs/projects/pharm_class/TTd_Oct29/PCL_group_rankpt_medians.txt'
pso.test_classes_incrementally(rnkpt_med_file=rnkptMed_file,n_test_max=40)