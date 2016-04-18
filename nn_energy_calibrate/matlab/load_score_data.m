function [data, tags] = load_score_data( file )

fid = fopen( file  );

found_score_labels = 0;

count = 0;

while ~feof( fid ) & ~found_score_labels
  l = fgetl( fid );
  
  if length(l) < 2 % problem!
    data = {};
    return;
  end

  cols = strsplit( l );
  %cols = split_string( l );
  
  if length( cols ) > 0 & strcmp( cols{1}, 'SCORE:' )
    if (~found_score_labels )
      found_score_labels = 1;
      score_labels = { cols{2: (length(cols)-1)} };
    end
  end
end

keyword = 'SCORE';
scorenums = [1:length( score_labels)] + 1;

command = ['grep ',keyword,' ',file,' | ',...
	   'awk ''{print $',num2str(scorenums(1))];
for i = 2: length( scorenums)
  command = [command, ',$',num2str(scorenums(i))];
end
command = [command,...
	   '}'' | grep -v inp | ',...
	   'grep -v R | grep -v H | grep -v score | grep -v pdb | grep -v descr | ', ...
	    ' grep -v total | grep -v S_ > data.txt'];

system(command);
scores = load('data.txt');

data.scores = scores;
data.score_labels = score_labels;

% get the tags 
% they are in column = length( scorenums ) +1
% grep -v description

tag_command = ['grep ',keyword,' ',file,' | ',...
		'awk ''{print $',num2str(length(scorenums)+2),'}'' | grep -v description > data_tags.txt'];
system(tag_command);

%tags = {}
ftagid = fopen( 'data_tags.txt'  );
tags = { fgetl( ftagid ) };
while ~feof( ftagid )
  tag = fgetl( ftagid );
  tags = { tags{1:length(tags)} tag };
end

fclose( ftagid );


fclose( fid );
