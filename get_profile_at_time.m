function prof_time = get_profile_at_time(profarray,ts,timeline)

prof_time=interp1(timeline,profarray,ts);

return 
end

