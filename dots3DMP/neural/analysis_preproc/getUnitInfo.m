
function unitInfo = getUnitInfo(thisDir, keepMU)

% import ch/depth, contamPct etc.

filename = fullfile(thisDir,'cluster_info.tsv');
fid = fopen(filename);
hdr  = textscan(fid, '%s%s%s%s%s%s%s%s%s%s%s',1);
main = textscan(fid,'%d%f%f%s%f%d%d%f%s%d%d');
fclose(fid);

if keepMU
    removethese = strcmp(main{9},'noise');
else
    removethese = ~strcmp(main{9},'good');
end

for f=1:length(main)
    main{f}(removethese) = [];
    unitInfo.(hdr{f}{:}) = main{f};
end


