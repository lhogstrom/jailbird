function matched_inds = get_cell_inds(cell_array,name_to_match,use_strcmpi)

if nargin < 3
   use_strcmpi = 0;
end

if use_strcmpi
    matched = cellfun(@(x) strcmpi(x,name_to_match),cell_array);
    matched_inds = find(matched == 1);
else
    matched = cellfun(@(x) regexpi(x,name_to_match),cell_array,'UniformOutput',0);
    notEmpties = ~cellfun(@isempty,matched);
    matched_inds = find(notEmpties==1);
end
