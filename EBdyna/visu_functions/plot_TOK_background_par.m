if strcmp(par.MACHINE_NAME,'MAST-U')
    plot_MASTU_TOK_background;
elseif strcmp(par.MACHINE_NAME,'COMPASS-U')
    plot_CU_TOK_background;
elseif strcmp(par.MACHINE_NAME,'COMPASS')
    plot_COMPASS_vessel;
end