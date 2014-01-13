import test_modules.pcla_svm_classifier as psc
import os

wkdir = '/xchip/cogs/projects/pharm_class/svm_pcla_classifier_Jan11'
if not os.path.exists(wkdir):
    os.mkdir(wkdir)

reload(psc)
pso = psc.svm_pcla(out=wkdir)
# pso.set_classes()
pso.set_clique_n69()
# pso.classification_by_cell()
pso.load_expression_data(max_signatures_per_cp=3,groups_to_model=None,keep_by_cell_line=True)
pso.classification_across_cell(loo_type='by_cp',n_procs=20)
pso.model_accuracy_description()
pso.write_results()
# rnkptMed_file = '/xchip/cogs/projects/pharm_class/TTd_Oct29/PCL_group_rankpt_medians.txt'
# pso.test_classes_incrementally(rnkpt_med_file=rnkptMed_file,n_test_max=40)