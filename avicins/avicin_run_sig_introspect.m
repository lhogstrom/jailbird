% run sig_introspect

sigGrp = '/xchip/cogs/projects/avicins//pert_explorer_avicins_Oct8/hydroxyl.grp';
out = '/xchip/cogs/projects/avicins//pert_explorer_avicins_Oct8';
% run introspect using wtcs and the default rank space and query group 
sig_introspect_tool('sig_id', sigGrp, 'out', out, 'metric', 'wtcs')

% compute self-ranks over a custom space 
% sig_introspect_tool('sig_id', 'sig_ids.grp', 'out', 'out_path', 'rank_space', 'rank_space_mcf7.grp', 'metric','wtcs')