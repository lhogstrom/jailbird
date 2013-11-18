   source("~/CGP2013/CNMF/CNMF.4.R")

   MSIG.Preprocess.Dataset(
      input.ds            = "~/CGP2013/TA/TA_PATHWAYS/TAref4pathways_OE_L1000_NMF_TA.OE003_A375_72H_20131012.gct",
      output.ds           = "~/CGP2013/TA/TA_PATHWAYS/TAref4pathways_OE_L1000_NMF_TA.OE003_A375_72H_20131012.NORM.gct",
      normalization       = 6)   # replace values with rank/total genes.

   MSIG.Preprocess.Dataset(
      input.ds            = "~/CGP2013/TA/TA_PATHWAYS/TAref4pathways_OE_L1000_NMF_TA.OE003_PC3_72H_20131012.gct",
      output.ds           = "~/CGP2013/TA/TA_PATHWAYS/TAref4pathways_OE_L1000_NMF_TA.OE003_PC3_72H_20131012.NORM.gct",
      normalization       = 6)   # replace values with rank/total genes.

   consensusNMF.2(
      input.ds            = "~/CGP2013/TA/TA_PATHWAYS/TAref4pathways_OE_L1000_NMF_TA.OE003_A375_72H_20131012.NORM.gct",
      k.init              = 2,
      k.final             = 11,
      num.clusterings     = 30,
      maxniter            = 1000,
      error.function      = "divergence",
      rseed               = 1898234,
      directory           = "~/CGP2013/TA/TA_PATHWAYS/",
      stopconv            = 20,
      stopfreq            = 10,
      non.interactive.run = F,
      doc.string          = "CNMF_TAref4pathways_OE_A375.v4")


#----------------------------------------


   source("~/CGP2013/CNMF/CNMF.4.R")

   MSIG.Preprocess.Dataset(
      input.ds            = "~/CGP2013/TA/CC/MCF7_top_intra_connecting_compound_classes_n79x978.gct",
      output.ds           = "~/CGP2013/TA/CC/MCF7_top_intra_connecting_compound_classes_n79x978.NORM.gct",
      normalization       = 6)   # replace values with rank/total genes.

   MSIG.Preprocess.Dataset(
      input.ds            = "~/CGP2013/TA/CC/PC3_top_intra_connecting_compound_classes_n83x978.gct",
      output.ds           = "~/CGP2013/TA/CC/PC3_top_intra_connecting_compound_classes_n83x978.NORM.gct",
      normalization       = 6)   # replace values with rank/total genes.

   consensusNMF.2(
      input.ds            = "~/CGP2013/TA/CC/PC3_top_intra_connecting_compound_classes_n83x978.NORM.gct",
      k.init              = 20,
      k.final             = 25,
      num.clusterings     = 30,
      maxniter            = 1000,
      error.function      = "divergence",
      rseed               = 1898234,
      directory           = "~/CGP2013/TA/CC/",
      stopconv            = 20,
      stopfreq            = 10,
      non.interactive.run = F,
      doc.string          = "CNMF_CC_PC3.v2")

#----------------------------------

   source("~/CGP2013/CNMF/CNMF.4.R")

   MSIG.Preprocess.Dataset(
      input.ds            = "~/CGP2012/signatures/L1000/signatures.QNORM_collapsed_to_symbols.gct",
      output.ds           = "~/CGP2012/signatures/L1000/signatures.QNORM_collapsed_to_symbols.NORM.gct",
      normalization       = 6)   # replace values with rank/total genes.

   consensusNMF.2(
      input.ds            = "~/CGP2012/signatures/L1000/signatures.QNORM_collapsed_to_symbols.NORM.gct",
      k.init              = 2,
      k.final             = 20,
      num.clusterings     = 30,
      maxniter            = 1000,
      error.function      = "divergence",
      rseed               = 1898234,
      directory           = "~/CGP2013/TA/Oncogenic/",
      stopconv            = 20,
      stopfreq            = 10,
      non.interactive.run = F,
      doc.string          = "CNMF_Oncogenic.v1")
