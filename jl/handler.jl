include("mcflux.jl")

skipsecs = 30;
in_dir = "data-kivalokh/";
out_dir = "jl-kivalokh/";
fig_dir = "figdir/";

co2outs = Output[];
ch4outs = Output[];

process_directory(in_dir,out_dir,fig_dir,co2outs,ch4outs,skipsecs,true,true);
