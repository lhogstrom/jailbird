%% test known drug-gene connections
% brd = 'BRD-K81418486'
baseDir = '/xchip/cogs/hogstrom/analysis/summly/bioAs';
resDir = fullfile(baseDir,brd,'sig_query');
outDir = fullfile(baseDir,brd,'summly_result');
mkdir(outDir)
sig_summly_tool(resDir, '--out', outDir)

