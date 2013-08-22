%% run summly tool query_tool result
resDir = '/xchip/cogs/projects/target_id/CTD2_14June2013/drug_KD_connection/MCF7/sig_query_out';
outDir = '/xchip/cogs/hogstrom/analysis/summly/24July2013';

sig_summly_tool(resDir, '--out', outDir)
% sig_summly_tool('/query_tool/folder/', '--out', '/outpath/')

%test for self connections
resDir = '/xchip/cogs/hogstrom/analysis/summly/trichostatin-a/jul24/my_analysis.sig_query_tool.2013072414472132/';
outDir = '/xchip/cogs/hogstrom/analysis/summly/trichostatin-a';
sig_summly_tool(resDir, '--out', outDir)

%% test outside signatures
metric = 'wtcs'
fup = '/xchip/cogs/hogstrom/analysis/summly/cflix_pannel/cflix_up_sm.gmt'
fdn = '/xchip/cogs/hogstrom/analysis/summly/cflix_pannel/cflix_dn_sm.gmt'
sig_query_tool('uptag', fup, 'dntag', fdn,'metric','wtcs','row_space','full')

resDir = '/xchip/cogs/hogstrom/analysis/summly/cflix_pannel_minus_phenox/jul26/my_analysis.sig_query_tool.2013072615304993';
outDir = '/xchip/cogs/hogstrom/analysis/summly/cflix_pannel/';
sig_summly_tool(resDir, '--out', outDir)


%% test known drug-gene connections
resDir = '/xchip/cogs/hogstrom/analysis/summly/bioAs/BRD-K81418486/sig_query';
outDir = '/xchip/cogs/hogstrom/analysis/summly/bioAs/BRD-K81418486/summly_result_2';
sig_summly_tool(resDir, '--out', outDir)

%% combined perts - more than one pert into summly 
% resDir = '/xchip/cogs/hogstrom/analysis/summly/cp_class/sig_query_no_7371';
resDir = '/xchip/cogs/hogstrom/analysis/summly/cp_class/sig_query';
outDir = '/xchip/cogs/hogstrom/analysis/summly/cp_class/sig_query';
sig_summly_tool(resDir, '--out', outDir)

%% avicin
% resDir = '/xchip/cogs/hogstrom/analysis/summly/avicin/sig_query_10um';
% outDir = '/xchip/cogs/hogstrom/analysis/summly/avicin';
resDir = '/xchip/cogs/projects/avicins/gold_10um/sig_query_10um';
outDir = '/xchip/cogs/projects/avicins/gold_10um/sig_query_10um';
sig_summly_tool(resDir, '--out', outDir)

%% cdt2
resDir = '/xchip/cogs/hogstrom/analysis/summly/cp_class/ctd2_sig_query';
outDir = '/xchip/cogs/hogstrom/analysis/summly/cp_class/ctd2_sig_query';
sig_summly_tool(resDir, '--out', outDir)

%% tp53 STITCH cps
resDir = '/xchip/cogs/hogstrom/analysis/summly/TP53_dgo/tp53_stitch_query';
outDir = '/xchip/cogs/hogstrom/analysis/summly/TP53_dgo/tp53_stitch_query';
sig_summly_tool(resDir, '--out', outDir)
