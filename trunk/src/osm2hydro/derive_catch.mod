# this PCRaster script derives:
# - a local drain direction map, 
# - a strahler order map, a reduced strahler order above a threshold,
# - a height above nearest drainage map,
# - a catchment map for each strahler threshold
# - a river stretch Id for each strahler subcatchment
# 
binding
 # Basin parameters and local drain direction map
dem          = $1; # maps_out\dem_burn.map
dem_src      = $2; # maps_out\dem_src.map;
Ldd          = $3;
drain        = $4;
# number of streamorder, used to delineate rivers
no_orders = 3;


timer
 1 1 1; # not relevant

dynamic

# Ldd = lddcreate(dem, 1e31,1e31,1e31,1e31);
strahler = streamorder(Ldd);
stream_thres = ordinal(mapmaximum(scalar(strahler))-no_orders); # Anything above this threshold is delineated as a subcatchment
strahler_thres = if(strahler ge stream_thres, strahler);
rivers_bool = boolean(strahler_thres);
Ldd_mask = if(strahler_thres ge 1, Ldd);
drain    = accuflux(Ldd_mask, 1);

# Parameters
alf = 50; # Alf ranges from 5 to > 60. 5 for hardrock. large values for sediments
Qavg = 170; # Yearly average Q at outlet (can be found online (Wikipedia))
Qmax = 6000; # Q at bankfull
n = 0.1; # Mannings roughness coefficient. Later to be extracted from land use maps?
degreelength = 120000; # number of meters of a one degree distance

# Create river width map (From Jaap Schellekens)
upstr = catchmenttotal(1, Ldd);
slope_avg = if(strahler_thres ge 1, slope(windowaverage(dem_src,celllength() * 16.0)))/degreelength;
Qavg_scale = upstr/mapmaximum(upstr) * Qavg;
Qmax_scale = upstr/mapmaximum(upstr) * Qmax;
W = (alf * (alf + 2.0)**(0.6666666667))**(0.375) * Qavg_scale**(0.375) * (max(slope_avg, 0.001))**(-0.1875) * n **(0.375);
RiverWidth = if(strahler_thres ge 1, W);

# Create embankment map
RiverEmb = if(strahler_thres ge 1, windowaverage(dem_src,celllength() * 16.0));

# Create depth map with chezy (replacing R by H assuming low H/W)
H = ((Qmax_scale * n)/(max(slope_avg, 0.001)**0.5 * W))**(3/5);
RiverDepth = if(strahler_thres ge 1, H);

# now make a point map where downstream value is not equal to value itself
upstream = upstream(Ldd_mask, scalar(rivers_bool));
downstream = downstream(Ldd_mask, upstream);
confluences = boolean(if(downstream eq 2, 1));
boundaries = boolean(if(scalar(Ldd_mask) eq 5, 1));
catch_points = nominal(uniqueid(cover(confluences, boundaries)));
catchments   = nominal(subcatchment(Ldd, catch_points));

report $5\rivers.map        =strahler_thres;
#report $5\ldd.map        =Ldd;
report $5\drain.map      = drain;
report $5\river_width.map      = RiverWidth;
report $5\river_depth.map      = RiverDepth;
report $5\river_embank.map      = RiverEmb;
report $5\slope.map      = slope_avg;
report $5\catch_point.map   = catch_points;
report $5\boundaries.map    = boundaries;


catchments2  = nominal(if(scalar(catchments) ge 1, catchments));
riversid     = if(strahler_thres ne 0, catchments2);
report $5\catchments.map    = catchments2;
report $5\riversid.map      = riversid;
