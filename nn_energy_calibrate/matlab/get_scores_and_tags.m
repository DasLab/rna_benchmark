function [data, tags] = get_scores_and_tags( scorefile  );
[data, tags] = load_score_data( scorefile  );
for n = 1:length(tags); 
  name = strsplit( tags{n}, '/' );
  tags{n} = strrep( name(length(name)), '_helix.pdb', '' );
  %tags{n} = strrep( basename( tags{n} ), '_helix.pdb', '' );
end;
