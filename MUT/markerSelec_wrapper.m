function snr_val = markerSelec_SNR(gct,out,varargin)
% Convert gctx to gct format

    pnames = {'class_field1','class1_name1','class1_name2','cls','edge_size','flip_contrast','std_fixflow'};
    %edge length
    dflts = {{}, {}, {}, {}, 50,false,false};

    args = parse_args(pnames, dflts, varargin{:});

pnames = {'in', 'out', 'rpt'};
dflts = {pwd, pwd, 'my_analysis'};
args = parse_args(pnames, dflts, varargin{:});
print_args(mfilename, 1, args);

[fn, fp] = find_file(args.in);
isgctx = ~cellfun(@isempty, regexp(fn,'.gctx$'));
fp = fp(isgctx);
fn = fn(isgctx);
nds = length(fp);
out = mktoolfolder(args.out, mfilename, 'rpt', args.rpt);
print_args(mfilename, fullfile(out, sprintf('%s_params.txt', mfilename)), args)
for ii=1:nds
    dbg (1, '%d / %d %s', ii, nds, fn{ii});
    ds = parse_gctx(fp{ii});
    [~, f, ~] = fileparts(fn{ii});
    mkgct(fullfile(out, f), ds); 
end

end