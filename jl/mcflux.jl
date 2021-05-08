## 2021-01-20
## Jani Anttila

using Dates;
using JSON;
using CSV;
using DataFrames;
using DataFramesMeta;
using GLM;
using Plots;
pyplot()

const dform0 = DateFormat("d.m.y");

const A = pi * (0.315/2.0)^2; ## kammion pohjapinta-ala
const V = A * 0.305;          ## kammion tilavuus
const Va = A * (0.305 + 0.2); ## lisäkauluksellisen kammion tilavuus
const Mco2 = 44.0095;         ## CO2 moolimassa
const Mch4 = 16.0425;         ## CH4 moolimassa
const Vm = 0.0224;

""" Calculates the flux value """
function getflux(area::Float64,volume::Float64,
                 temp::Float64,molmass::Float64)::Float64
    vm = 0.0224;
    zerot = 273.15;
    convf = 3600.0; # [ug (gas) m-2 h-1]
    return (molmass/vm)*(zerot/(zerot + temp))*(volume/area)*convf
end

"""
Calculates residual sum-of-squares variation from a yv~xv linear fit
with coefficients c0, c1. Skip seconds can be provided with start_idx.
"""
function residual_variation(c0::Float64,c1::Float64,start_idx::Int64,
                            xv::Array{Int64,1},yv::Array{Float64,1})::Float64
    resid_sq_sum = 0.0;
    for i in start_idx:length(xv)
        resid = (c0 + c1 * xv[i]) - yv[i];
        resid_sq_sum += (resid * resid);
    end
    return sqrt(resid_sq_sum)
end

mutable struct Output
    lin_flux_s0::Array{Float64}   ## linear flux without skipping
    lin_flux_skip::Array{Float64} ## linear flux with skipping
    residual::Array{Float64}      ## total residual variation (w skips)
    start_pp::Array{Float64}      ## start ppm/ppb (w skips)
    end_pp::Array{Float64}        ## end ppm/ppb
    icept::Array{Float64}         ## intercept (w skips)
    lcoef::Array{Float64}         ## linear coefficient (w skips)
end

""" Creates an empty output object """
function make_output(n::Integer)::Output
    o = Output(zeros(n),zeros(n),zeros(n),zeros(n),
               zeros(n),zeros(n),zeros(n));
    return o;
end

"""
Reads mdf: metadata file (i.e. field form) and ddf: data file from licor, 
and fills Output objects for co2o: carbon dioxide and ch4o: methane. 
Seconds to be skipped from the beginning are supplied as skips.
"""
function fill_out_obj(mdf::DataFrame,ddf::DataFrame,
                      co2o::Output,ch4o::Output,
                      skips::Int64)
    for r = 1:nrow(mdf)
        ##println(r)
        day = mdf[r,"Date"];
        start_time = mdf[r,"Start time"];
        end_time = mdf[r,"End time"];
        ## sometimes times are read in as Strings -> try to infer the Time objects
        if typeof(start_time) == String; start_time = Time(start_time); end
        if typeof(end_time) == String; end_time = Time(end_time); end
        ## find the time-series window from the data file
        subs = filter(row -> ((row.DATE == day) &&
                              (row.TIME >= start_time) &&
                              (row.TIME <= end_time)), ddf);
        ## set :secs column to start seconds from zero
        mintime = subs[1,:SECONDS];
        transform!(subs, :SECONDS => (x -> x .- mintime) => :secs);
        ## get duration in rows/seconds
        full_duration = nrow(subs);
        duration = full_duration - skips + 1;
        ## set Output obj start and end values: ppm for co2 and ppb for ch4
        co2o.start_pp[r] = subs[skips,:CO2];
        co2o.end_pp[r] = subs[full_duration,:CO2];
        ch4o.start_pp[r] = subs[skips,:CH4];
        ch4o.end_pp[r] = subs[full_duration,:CH4];
        ## get mean temp during measurement
        T_mean = (mdf[r,"Chamber T_start, C"] + mdf[r,"Chamber T_end, C"]) / 2.0;
        ## fit linear regression coefs for no skips case
        co2lms0 = lm(@formula(CO2 ~ secs),subs[1:full_duration,:]);
        ch4lms0 = lm(@formula(CH4 ~ secs),subs[1:full_duration,:]);
        co2cs0 = coef(co2lms0); ch4cs0 = coef(ch4lms0);
        ## fit linear regression coefs for skips case
        co2lm = lm(@formula(CO2 ~ secs),subs[skips:full_duration,:]);
        ch4lm = lm(@formula(CH4 ~ secs),subs[skips:full_duration,:]);
        co2c = coef(co2lm); ch4c = coef(ch4lm);
        ## store the flux values
        ## also: handle the very annoying additional collar case
        ## todo: this could be much neater...
        if !ismissing(mdf[r,"Notice!"]) && !isnothing(findfirst("LISÄKAULUS",mdf[r,"Notice!"]))
            co2o.lin_flux_s0[r] = round((co2cs0[2] *
                                         getflux(A,Va,T_mean,Mco2)) / 1000,
                                        digits=3);
            ch4o.lin_flux_s0[r] = round((ch4cs0[2] *
                                         getflux(A,Va,T_mean,Mch4)),
                                        digits=3);
            co2o.lin_flux_skip[r] = round((co2c[2] *
                                           getflux(A,Va,T_mean,Mco2)) / 1000,digits=3);
            ch4o.lin_flux_skip[r] = round((ch4c[2] *
                                           getflux(A,Va,T_mean,Mch4)),digits=3);
            ##println("using lisäkaulus at row ",r)
        else
            ## [mg CO2 m-2 h-1]
            co2o.lin_flux_s0[r] = round((co2cs0[2] *
                                         getflux(A,V,T_mean,Mco2)) / 1000,
                                        digits=3);
            ch4o.lin_flux_s0[r] = round((ch4cs0[2] *
                                         getflux(A,V,T_mean,Mch4)),
                                        digits=3);
            co2o.lin_flux_skip[r] = round((co2c[2] *
                                           getflux(A,V,T_mean,Mco2)) / 1000,digits=3);
            ch4o.lin_flux_skip[r] = round((ch4c[2] *
                                           getflux(A,V,T_mean,Mch4)),digits=3);
        end
        ## calculate residual variation for the skips case
        co2rv = residual_variation(co2c[1],co2c[2],skips,subs[:,:secs],subs[:,:CO2]);
        ch4rv = residual_variation(ch4c[1],ch4c[2],skips,subs[:,:secs],subs[:,:CH4]);
        co2o.residual[r] = round(co2rv,digits=3);
        ch4o.residual[r] = round(ch4rv,digits=3);
        ## store the (skips) linear coefficients for plotting later
        co2o.icept[r] = co2c[1]; co2o.lcoef[r] = co2c[2];
        ch4o.icept[r] = ch4c[1]; ch4o.lcoef[r] = ch4c[2];
    end
    return nothing;
end

function plot_out_obj(mdf::DataFrame,ddf::DataFrame,modnum::String,basename::String,
                      fdir::String,co2o::Output,ch4o::Output,skips::Int64)
    for mr = 1:nrow(mdf)
        day = mdf[mr,"Date"];
        start_time = mdf[mr,"Start time"]
        end_time = mdf[mr,"End time"]
        if typeof(start_time) == String; start_time = Time(start_time); end
        if typeof(end_time) == String; end_time = Time(end_time); end
        dateobj = mdf[mr,"Date"]
        sitename = mdf[mr,"Site"]
        plotname = mdf[mr,"Plot"]
        co2titlestring = string("co2-",dateobj,"-",sitename,"-",plotname,"-",modnum,
                                "-",start_time)
        ch4titlestring = string("ch4-",dateobj,"-",sitename,"-",plotname,"-",modnum,
                                "-",start_time)
        psub = filter(row -> ((row.DATE == day) &&
                              (row.TIME >= start_time) &&
                              (row.TIME <= end_time)), ddf);
        mintime = psub[1,:SECONDS];
        transform!(psub, :SECONDS => (x -> x .- mintime) => :secs);
        full_duration = nrow(psub); duration = full_duration - skipsecs + 1;
        init = psub[1:skipsecs,:]; fitp = psub[skipsecs:full_duration,:];
        co2icept = co2o.icept[mr]; co2lcoef = co2o.lcoef[mr];
        ch4icept = ch4o.icept[mr]; ch4lcoef = ch4o.lcoef[mr];
        lxv = skipsecs:full_duration;
        co2yxv = map(x -> co2icept + x*co2lcoef, lxv);
        ch4yxv = map(x -> ch4icept + x*ch4lcoef, lxv);
        figname = string(fdir,basename,"-",lpad(mr,4,"0"),".pdf")
        p1 = plot(init[:,:secs], init[:,:CO2], seriestype = :scatter, color = :gray,
                  marker = (:hex, 2, 0.8, Plots.stroke(1, :gray)), label = "skipped",
                  title = co2titlestring, show = false,
                  annotation = (0,maximum(fitp[:,:CO2])-2,
                                Plots.text(string("RMSE = ",co2o.residual[mr]),
                                           12,:black,:left)));
        plot!(p1,fitp[:,:secs], fitp[:,:CO2], seriestype = :scatter, color = :red,
              marker = (:hex, 3, 0.8, Plots.stroke(1, :gray)), label = "fitted",
              show = false);
        plot!(p1,lxv, co2yxv, color = :blue, label = "fit", show = false);
        ylabel!("ppm"); xlabel!("seconds");
        p2 = plot(init[:,:secs], init[:,:CH4], seriestype = :scatter, color = :gray,
                  marker = (:hex, 2, 0.8, Plots.stroke(1, :gray)), label = "skipped",
                  title = ch4titlestring, show = false,
                  annotation = (0,minimum(fitp[:,:CH4])+1,
                                Plots.text(string("RMSE = ",ch4o.residual[mr]),
                                           12,:black,:left)));
        plot!(p2,fitp[:,:secs], fitp[:,:CH4], seriestype = :scatter, color = :red,
              marker = (:hex, 3, 0.8, Plots.stroke(1, :gray)), label = "fitted",
              show = false);
        plot!(p2,lxv, ch4yxv, color = :blue, label = "fit", show = false);
        ylabel!("ppb"); xlabel!("seconds");
        plot(p1, p2, layout = (1, 2), size = (1280,480), legend = false, show = false);
        savefig(figname);
    end
end


"""
Reads all mfiles: metadata (i.e. field form) files and dfiles: data files 
in directory. Fills arrays of co2 and ch4 Output structs. Writes corresponding
csv ouput files with bnames: filenames without extensions, to outdirectory.
The metadata file extension should be .csv and the datafile extension .data.
"""
function process_directory(directory::String,outdirectory::String,figdir::String,
                           co2os::Array{Output},ch4os::Array{Output},
                           skipsecs::Int64,writecsv::Bool,plotflux::Bool)
    dd_contents = readdir(directory);
    basenames = String[];
    modelnums = String[];
    metafiles = String[];
    datafiles = String[];
    for i in 1:length(dd_contents)
        fn = dd_contents[i];
        extension = SubString(fn,length(fn)-3,length(fn));
        if(extension == ".csv")
            push!(metafiles,fn)
            basename = SubString(fn,1,length(fn)-4);
            modelnum = SubString(basename,length(basename)-4,length(basename));
            push!(basenames,basename);
            push!(modelnums,modelnum);
            push!(datafiles,string(basename,".data"));
        end
    end
    for f in 1:length(metafiles)
        println("file ",metafiles[f]);
        metaname = metafiles[f];
        dataname = datafiles[f];
        metapath = string(directory,metaname);
        datapath = string(directory,dataname);
        metadf = CSV.File(metapath) |> DataFrame;
        datadf = CSV.File(datapath, header=6, datarow=8) |> DataFrame;
        ## field forms often contain empty rows -> remove rows without 'Start time'
        dropmissing!(metadf,"Start time");
        num_meas = nrow(metadf);
        ## sometimes dates are read wrong and it may help to re-read them
        metadf.Date = Date.(metadf.Date,Dates.DateFormat("dd.mm.yyyy"));
        if(typeof(datadf.DATE) != Array{Date,1})
            datadf.DATE = Date.(datadf.DATE,Dates.DateFormat("dd.mm.yyyy"));
        end
        co2o = make_output(num_meas);
        ch4o = make_output(num_meas);
        fill_out_obj(metadf,datadf,co2o,ch4o,skipsecs);
        push!(co2os,co2o);
        push!(ch4os,ch4o);
        if plotflux
            plot_out_obj(metadf,datadf,modelnums[f],basenames[f],figdir,co2o,ch4o,skipsecs);
        end
        if writecsv
            metadf.co2_lin_flux_s0   = co2o.lin_flux_s0;
            metadf.co2_lin_flux_skip = co2o.lin_flux_skip;
            metadf.co2_residual      = co2o.residual;
            metadf.co2_start_ppm     = co2o.start_pp;
            metadf.co2_end_ppm       = co2o.end_pp;
            metadf.ch4_lin_flux_s0   = ch4o.lin_flux_s0;
            metadf.ch4_lin_flux_skip = ch4o.lin_flux_skip;
            metadf.ch4_residual      = ch4o.residual;
            metadf.ch4_start_ppb     = ch4o.start_pp;
            metadf.ch4_end_ppb       = ch4o.end_pp;
            CSV.write(string(outdirectory,basenames[f],"-flux.csv"),metadf);
        end
    end
    return nothing
end
