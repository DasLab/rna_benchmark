function out_color = fade_color(  start_color, fade_level );

out_color =  start_color + fade_level * ( [1 1 1] - start_color );

out_color = min( out_color, 1.0 );