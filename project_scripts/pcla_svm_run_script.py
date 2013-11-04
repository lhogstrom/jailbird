import test_modules.pcla_svm_classifier as psc

wkdir = '/xchip/cogs/projects/pharm_class/svm_pcla_classifier_NOV1'

reload(psc)
pso = psc.svm_pcla(out=wkdir)
pso.set_classes()
pso.classification_by_cell()