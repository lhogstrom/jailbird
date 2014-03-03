##### Annotate components
#
#  1. Project the CCLE in the space of components (generate an H matrix)
#  2. Make a heatmap of the projection
#  3. Annotate each components by matching genomic features to it
   # source("/xchip/cogs/hogstrom/analysis/pablos_NMF_analysis/TA/CNMF.4.R")
   source("/xchip/cogs/projects/NMF/Comp_Annot/FS.library.v8.6.R")
   source("/xchip/cogs/projects/NMF/Comp_Annot/DISSECTOR_lib.v3.R")

   # path1 <- "/xchip/cogs/projects/NMF/MCF7_7_PCLs_w_DMSO"
   # prefix1 <- "MCF7_top_intra_connecting_compound_classes_n129x978"
   # outpath <- "/xchip/cogs/projects/NMF/MCF7_comp_annot_to_CCLE_space2"
   prefix1 <- "PC3_top_intra_connecting_compound_classes_n134x978"
   outpath <- "/xchip/cogs/projects/NMF/PC3_comp_annot_to_CCLE_space"

   DISSECTOR_project_dataset.v1(     # Project a dataset in the space defined by a W matrix
      # input_dataset           = paste(path1,"/",prefix1,".gct",sep=""),  # Input dataset (GCT)
      input_dataset           = "/xchip/cogs/projects/NMF/Comp_Annot/rnaseq.v3.gct",  # Input dataset (GCT)      
      input_normalization     = "rank",     # Normalization for the input dataset: "rank"
      normalize_after_match   = T,        # Normalize input dataset after matching with rows of W
      input_W_dataset         = paste(outpath,"/",prefix1,".W.k9.gct",sep=""),  # Input W matrix (GCT)
      W_normalization         = "none",         # Normalization for W                                            
      output_H_dataset        = paste(outpath,"/",prefix1,".H_proj.v1.gct",sep=""),   # Output dataset H (GCT)
      output_W_dataset        = paste(outpath,"/",prefix1,".W.v1.gct",sep=""))   # Output dataset normalized W (GCT) - 932 tissue types x components

   source("/xchip/cogs/projects/NMF/Comp_Annot/DISSECTOR_lib.v3.R")

   DISSECTOR_make_heatmap_of_matrix.v1(
      input_dataset            = paste(outpath,"/",prefix1,".H_proj.v1.gct",sep=""), # Input dataset (GCT).
      # annot_file               = c(paste(path1,"/MCF7_top_intra_connecting_compound_classes.v3.txt",sep=""),"sig_id", "annot", F), # Phenotype annotation file (TXT, optional) in format c(file, name_column, annot_column, use_prefix)
      annot_file               = c("/xchip/cogs/projects/NMF/Comp_Annot/CCLE_Kim_Tamayo_96AffyCEL_FAZES_Oct12.v2.txt","CCLE_name", "Site.Primary", F),
      transpose_data           = F,           # Transpose input matrix
      append_annot             = T,           # Append annotation to column names
      sorting_method           = "HC",        # Sorting method for cols inside phenotype: MDS (Multi_dimensinal Scaling) or HC (Hiererachical Clustering)  
      output_plot_landscape    = paste(outpath,"/PC3.H_proj.v1_LPLOT.v1.pdf",sep=""),  # Output (PDF) file
      output_plot_portrait     = paste(outpath,"/PC3.H_proj.v1_PPLOT.v1.pdf",sep=""))    # Output (PDF) file

   source("/xchip/cogs/projects/NMF/Comp_Annot/DISSECTOR_lib.v3.R")

   DISSECTOR_annotate_components.v1(  #  Match genomic features to each components profile
      input_dataset         = paste(outpath,"/",prefix1,".H_proj.v1.gct",sep=""), # Input dataset (GCT).
      directory             = paste(outpath,"/c_annot",sep=""),
      identifier            = "Comp_annot.v1",                                                            # string or prefix to identify this analysis
      feature.type.files    =  list("ACHILLES"    = "/xchip/cogs/projects/NMF/Comp_Annot/Achilles_v2.4.1.rnai.Gs.gct",
                                   "RPPA"         = "/xchip/cogs/projects/NMF/Comp_Annot/RPPA.dat.gct", #phospho protein
                                   "MUT_CNA"      = "/xchip/cogs/projects/NMF/Comp_Annot/RNAseqMUTs_CNA_20130729.gct",  
                                   "EXP_PATHWAYS" = "/xchip/cogs/projects/NMF/Comp_Annot/CCLE_MSigDB_plus_oncogenic.PATHWAYS.v2.gct",
                                   "EXP_GENES"    = "/xchip/cogs/projects/NMF/Comp_Annot/rnaseq.v3.gct"),
       feature.dir          = c(0, 1, 1, 1, 1),   # 0=negative, 1=positive
       # feature.dir          = c(0, 1), # 0=negative, 1=positive
       n.markers            = 25,                         # Number of top hits shown in the heatmaps
       n.perm               = 2,                                         
       char.scaling         = 0.575,                   # Character scaling for heatmaps
       locs.table.file      = "/xchip/cogs/projects/NMF/Comp_Annot/hgnc_downloads.txt",
       # log.table.file       = "~/CGP2013/TA/TA_PATHWAYS/CCLE_TAref4pathways_k7_4/annot_comp.log.txt",                            
       log.table.file       = paste(outpath,"/annot_comp.log.txt",sep=""),
       min.thres            = 10,
       character.scaling    = 0.65,
       phen.table           = NULL,
       phenotypes           = NULL)



