##### Annotate components
#
#  1. Project the CCLE in the space of components (generate an H matrix)
#  2. Make a heatmap of the projection
#  3. Annotate each components by matching genomic features to it

   source("~/CGP2013/TA/DISSECTOR_lib.v3.R")

   DISSECTOR_project_dataset.v1(     # Project a dataset in the space defined by a W matrix
      input_dataset           = "~/CGP2013/CCLE/rnaseq.v3.gct",  # Input dataset (GCT)
      input_normalization     = "rank",     # Normalization for the input dataset: "rank"
      normalize_after_match   = T,        # Normalize input dataset after matching with rows of W
      input_W_dataset         = "~/CGP2013/TA/TA_PATHWAYS/CCLE_TAref4pathways_k7_4/TAref4pathways_A375.W.v1.gct",  # Input W matrix (GCT)
      W_normalization         = "none",         # Normalization for W                                            
      output_H_dataset        = "~/CGP2013/TA/TA_PATHWAYS/CCLE_TAref4pathways_k7_4/rnaseq.v3_TAref4pathways_A375.H_proj.v1.gct",   # Output dataset H (GCT)
      output_W_dataset        = "~/CGP2013/TA/TA_PATHWAYS/CCLE_TAref4pathways_k7_4/TAref4pathways_A375.W.v1.gct")  # Output dataset normalized W (GCT)

   source("~/CGP2013/TA/DISSECTOR_lib.v3.R")

   DISSECTOR_make_heatmap_of_matrix.v1(
      input_dataset            = "~/CGP2013/TA/TA_PATHWAYS/CCLE_TAref4pathways_k7_4/rnaseq.v3_TAref4pathways_A375.H_proj.v1.gct", # Input dataset (GCT).
      annot_file               = c("~/CGP2013/signatures/components3/datasets/CCLE_Kim_Tamayo_96AffyCEL_FAZES_Oct12.v2.txt","CCLE_name", "Site.Primary", F),
      transpose_data           = F,           # Transpose input matrix
      append_annot             = T,           # Append annotation to column names
      sorting_method           = "HC",        # Sorting method for cols inside phenotype: MDS (Multi_dimensinal Scaling) or HC (Hiererachical Clustering)  
      output_plot_landscape    = "~/CGP2013/TA/TA_PATHWAYS/CCLE_TAref4pathways_k7_4/rnaseq.v3_TAref4pathways_A375.H_proj.v1_LPLOT.v1.pdf",  # Output (PDF) file
      output_plot_portrait     = "~/CGP2013/TA/TA_PATHWAYS/CCLE_TAref4pathways_k7_4/rnaseq.v3_TAref4pathways_A375.H_proj.v1_PPLOT.v1.pdf")  # Output (PDF) file

   source("~/CGP2013/TA/DISSECTOR_lib.v3.R")

   DISSECTOR_annotate_components.v1(  #  Match genomic features to each components profile
      input_dataset         = "~/CGP2013/TA/TA_PATHWAYS/CCLE_TAref4pathways_k7_4/rnaseq.v3_TAref4pathways_A375.H_proj.v1.gct", # Input dataset (GCT).
      directory             = "~/CGP2013/TA/TA_PATHWAYS/CCLE_TAref4pathways_k7_4/",                            
      identifier            = "Comp_annot.v1",                                                            # string or prefix to identify this analysis
      feature.type.files    =  list("ACHILLES"     = "~/CGP2013/Distiller/Achilles_v2.4.1.rnai.Gs.gct",
                                    "RPPA"         = "~/CGP2013/Distiller/RPPA.dat.gct"),
                                  #  "MUT_CNA"      = "~/CGP2013/Distiller/RNAseqMUTs_CNA_20130729.gct",  
                                  #  "EXP_PATHWAYS" = "~/CGP2013/Distiller/CCLE_MSigDB_plus_oncogenic.PATHWAYS.v2.gct",
                                  #  "EXP_GENES"    = "~/CGP2013/CCLE/rnaseq.v3.gct"),
       feature.directions   = c(0, 1, 1, 1, 1),                                         
       n.markers            = 25,                         # Number of top hits shown in the heatmaps
       n.perm               = 2,                                         
       char.scaling         = 0.575,                   # Character scaling for heatmaps
       locs.table.file      = "~/CGP2013/Distiller/hgnc_downloads.txt",
       log.table.file       = "~/CGP2013/TA/TA_PATHWAYS/CCLE_TAref4pathways_k7_4/annot_comp.log.txt",                            
       min.thres            = 10,
       character.scaling    = 0.65,
       phen.table           = NULL,
       phenotypes           = NULL)  

