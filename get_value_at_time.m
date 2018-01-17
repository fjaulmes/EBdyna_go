function val_time = get_value_at_time(valarray,ts,timeline)

val_time=interp1(timeline,valarray,ts);

return 
end