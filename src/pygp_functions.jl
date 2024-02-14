#######################################################
# Parsers
#######################################################


# User types in lat/long points, that IDs the closest polygon?
# 


module pygp_functions
__precompile__(false)  # will cause using / import to load it directly into the 
                       # current process  skip the precompile  caching. 
                       # This also thereby prevents the module from being 
                       # imported by any other precompiled module.
                       # https://docs.julialang.org/en/v1/manual/modules/

print("\n\nStarting modulen pygp_functions'...loading dependencies...\n")
using BenchmarkTools # for @time
using InvertedIndices # for Not
using LSODA           # for lsoda()
using Sundials        # for CVODE_BDF(linear_solver=:GMRES)
using DifferentialEquations
using Distributed
using Random					# for MersenneTwister()
using Dates						# for e.g. DateTime, Dates.now()
using PhyloNetworks
using CSV
#using Plots						# for plot
using DataFrames          # for DataFrame()
using Query               # for Querying dataframes!
using BioGeoJulia.TrUtils # for flat2() (similar to unlist)
using BioGeoJulia.StateSpace
using BioGeoJulia.TreePass
using BioGeoJulia.SSEs

# (1) List all function names here:
export locations_used, land_last_touch, distance_given, distance_interp, oldest_polygon, pygptxt_read, pygpcsv_read


# DataValue{ANY} comes from the fact that it's querying through 



#######################################################
# Temporary file to store functions under development
#
# Start with:
# 
# Setup:

"""
cd("/GitHub/BioGeoJulia.jl/notes/")
includen pygp_functions.jl")
"""
#######################################################


#######################################################
# (2) write the functions here
#######################################################

"""
    locations_used(df)

Provides a String Vector of all Land Names used by GPlates printer within the pygplates dataframe 
ONLY USABLE FOR DATAFRAMES CREATED USING PYGPLATES OUTPUT CREATED BY WALLIS BLAND. (SEE GPLATES PRINTER FOR PYTHON DOCUMENTATION)

* `df` - dataframe created by pygplates output created by Wallis Bland
contains variables: 'Reconstruction_Time_Ma', 'Land1', 'Land2', 'Closest_Distance_Ma',
'Point1_Lat', 'Point1_Lon', 'Point2_Lat', 'Point2_Lon'

# Examples
```julia-repl

julia> test_df = DataFrame(Reconstruction_Time_Ma=[0, 10, 20, 30, 40, 50, 60, 0, 10, 20, 30, 40, 50, 60], Land1=["India", "India", "India", "India", "India", "India", "India", "India", "India", "India", "India", "India", "India", "India"], Land2=["Asia", "Asia", "Asia", "Asia", "Asia", "Asia", "Asia", "Madagascar", "Madagascar", "Madagascar", "Madagascar", "Madagascar", "Madagascar","Madagascar"], Closest_Distance_km=[0, 1000, 2000, 3500, 4100, 4900, "NA", 4900, 4100, 3500, 2000, 1000, 0, "NA"], Point1_Lat=[1,2,3,4,5,6,"NA",1,2,3,4,5,6,"NA"], Point1_Lon=[1,2,3,4,5,6,"NA",1,2,3,4,5,6,"NA"], Point2_Lat=[1,2,3,4,5,6,7,1,2,3,4,5,6,"NA"], Point2_Lon=[1,2,3,4,5,6,7,1,2,3,4,5,6,"NA"])

14×8 DataFrame
 Row │ Reconstruction_Time_Ma  Land1   Land2       Closest_Distance_km  Point1_Lat  Point1_Lon  Point2_Lat  Point2_Lon 
     │ Int64                   String  String      Any                  Any         Any         Any         Any        
─────┼─────────────────────────────────────────────────────────────────────────────────────────────────────────────────
   1 │                      0  India   Asia        0                    1           1           1           1
   2 │                     10  India   Asia        1000                 2           2           2           2
   3 │                     20  India   Asia        2000                 3           3           3           3
   4 │                     30  India   Asia        3500                 4           4           4           4
   5 │                     40  India   Asia        4100                 5           5           5           5
   6 │                     50  India   Asia        4900                 6           6           6           6
   7 │                     60  India   Asia        NA                   NA          NA          7           7
   8 │                      0  India   Madagascar  4900                 1           1           1           1
   9 │                     10  India   Madagascar  4100                 2           2           2           2
  10 │                     20  India   Madagascar  3500                 3           3           3           3
  11 │                     30  India   Madagascar  2000                 4           4           4           4
  12 │                     40  India   Madagascar  1000                 5           5           5           5
  13 │                     50  India   Madagascar  0                    6           6           6           6
  14 │                     60  India   Madagascar  NA                   NA          NA          NA          NA

julia> locations_used(test_df)
3-element Vector{String}:
   "India"
   "Asia"
   "Madagascar"

julia>   df = DataFrame(name=["John", "Sally", "Roger", "Alice"], age=[54., 34., 79., 41.], children=[0, 2, 4, 5], adult=["Alice", "Bob", "Charlie", "John"])

julia> locations_used(df)
ERROR: ArgumentError: column name :Land1 not found in the data frame

# POTENTIAL ERRORS:

* Land1 or Land2 not found 
    Not using the correct dataframe, or have edited the original PyGplates output headers

    Be sure all variable headers align with original given headers in GPlates Printer code. 
    Other headings will not work with all functions in this function set

    julia> df = DataFrame(name=["John", "Sally", "Roger", "Alice"], age=[54., 34., 79., 41.], children=[0, 2, 4, 5], adult=["Alice", "Bob", "Charlie", "John"])
    julia> locations_used(df)
    ERROR: ArgumentError: column name :Land1 not found in the data frame


```

"""


function locations_used(df)
    lands = append!(df.Land1, df.Land2)
    return unique(lands)
end


"""
    distance_given(df, time, land1, land2)

Provides the distance at a given time point, provided that the time point requested is
along the interval originally outputted from pygplates into the dataframe. Allows the 
user to quickly access distance already recorded within their dataframe. 

For interpolated distances at time points not previously recorded, please use 'distance_interp'

* `df` - dataframe created by pygplates output created by Wallis Bland
contains variables: 'Reconstruction_Time_Ma', 'Land1', 'Land2', 'Closest_Distance_Ma',
'Point1_Lat', 'Point1_Lon', 'Point2_Lat', 'Point2_Lon'

* `time` - timepoint for distance requested. Within distance_given, this timepoint should be
within the original dataframe.

* `land1` - first landmass chosen for distance comparison. Ensure this is within "string" format.

* `land2` - second landmass chosed for distance comparison. Ensure this is within "string" format

NOTE: Lands 1 and 2 do not need to line up with Land1 or 2 within the dataframe in this function, as it will search both!

# Examples
```julia-repl

julia> test_df = DataFrame(Reconstruction_Time_Ma=[0, 10, 20, 30, 40, 50, 60, 0, 10, 20, 30, 40, 50, 60], Land1=["India", "India", "India", "India", "India", "India", "India", "India", "India", "India", "India", "India", "India", "India"], Land2=["Asia", "Asia", "Asia", "Asia", "Asia", "Asia", "Asia", "Madagascar", "Madagascar", "Madagascar", "Madagascar", "Madagascar", "Madagascar","Madagascar"], Closest_Distance_km=[0, 1000, 2000, 3500, 4100, 4900, "NA", 4900, 4100, 3500, 2000, 1000, 0, "NA"], Point1_Lat=[1,2,3,4,5,6,"NA",1,2,3,4,5,6,"NA"], Point1_Lon=[1,2,3,4,5,6,"NA",1,2,3,4,5,6,"NA"], Point2_Lat=[1,2,3,4,5,6,7,1,2,3,4,5,6,"NA"], Point2_Lon=[1,2,3,4,5,6,7,1,2,3,4,5,6,"NA"])

14×8 DataFrame
 Row │ Reconstruction_Time_Ma  Land1   Land2       Closest_Distance_km  Point1_Lat  Point1_Lon  Point2_Lat  Point2_Lon 
     │ Int64                   String  String      Any                  Any         Any         Any         Any        
─────┼─────────────────────────────────────────────────────────────────────────────────────────────────────────────────
   1 │                      0  India   Asia        0                    1           1           1           1
   2 │                     10  India   Asia        1000                 2           2           2           2
   3 │                     20  India   Asia        2000                 3           3           3           3
   4 │                     30  India   Asia        3500                 4           4           4           4
   5 │                     40  India   Asia        4100                 5           5           5           5
   6 │                     50  India   Asia        4900                 6           6           6           6
   7 │                     60  India   Asia        NA                   NA          NA          7           7
   8 │                      0  India   Madagascar  4900                 1           1           1           1
   9 │                     10  India   Madagascar  4100                 2           2           2           2
  10 │                     20  India   Madagascar  3500                 3           3           3           3
  11 │                     30  India   Madagascar  2000                 4           4           4           4
  12 │                     40  India   Madagascar  1000                 5           5           5           5
  13 │                     50  India   Madagascar  0                    6           6           6           6
  14 │                     60  India   Madagascar  NA                   NA          NA          NA          NA


julia> distance_given(test_df, 40, "India", "Asia")
4100

julia> distance_given(test_df, 40, "Asia", "India")
4100


# POTENTIAL ERRORS:

* TTIMESTAMP NOT WITHIN TABLE
    Variables potentially identified, but time stamp not within the written intervals of the dataframe. Use function pygp_functions.distance_interp() instead.

    julia> distance_given(test_df, 44, "Asia", "India")
    ERROR: STOP ERROR in pygp_functions.distance_given().
     Pairing: 'Asia' and 'India' has not been found within given dataframe at this timestamp.
     Please check spelling of compared landmasses.
     If using function distance_given(), please check dataframe timestamps.
 	    if time intervals taken from original pygplates do not have your needed timestamp, please us function pygp_functions.distance_interp() 
     If using function distance_interp(), please ensure time intervals within dataframe go PAST the requested timestamp. Thank you! 



* INCORRECT VARIABLES / MISSPELLINGS
    In this case the variable was either spelt wrong or for some other reason Variable was never printed within the given dataframe.

    julia> distance_given(test_df, 40, "Asia", "Indai")
    ERROR: STOP ERROR in pygp_functions.distance_given().
     Pairing: 'Asia' and 'Indai' has not been found within given dataframe at this timestamp.
     Please check spelling of compared landmasses.
     If using function distance_given(), please check dataframe timestamps.
        if time intervals taken from original pygplates do not have your needed timestamp, please us function pygp_functions.distance_interp() 
     If using function distance_interp(), please ensure time intervals within dataframe go PAST the requested timestamp. Thank you! 



* INCORRECT PAIRING
    Lands exist within the dataframe, but were never compared distance wise 
    If original pygplates code was run correctly, this SHOULD NOT OCCUR. More likely need to check spelling

    julia> distance_given(test_df, 40, "Asia", "Madagascar")
    ERROR: STOP ERROR in pygp_functions.distance_given().
     Pairing: 'Asia' and 'Madagascar' has not been found within given dataframe at this timestamp.
     Please check spelling of compared landmasses.
     If using function distance_given(), please check dataframe timestamps.
        if time intervals taken from original pygplates do not have your needed timestamp, please us function pygp_functions.distance_interp() 
     If using function distance_interp(), please ensure time intervals within dataframe go PAST the requested timestamp. Thank you! 



* BOTH LANDS MISSING POLYGONS
    Variables identified, but polygons for original gplates gpml file for BOTH landmasses were nvever identified at that timestamp

    julia> distance_given(test_df, 60, "India", "Madagascar")
    ERROR: STOP ERROR in pygp_functions.distance_given(). Both landmasses, 'India' and 'Madagascar', 
     do not have polygons at this timestamp (time=60).
     Please use land_begin to find the earliest timestamp with polygons present. If land_begin returns an earlier timestamp, please check gplates file.



* SINGLE LAND MISSING POLYGONS
    Variables identified, but polygons from original gplates gpml file were never identified at that timestamp.
    Note below, whether India is land 1 or 2, it is noted that that is the one missing polygons


    julia> distance_given(test_df, 60, "India", "Asia") 
    ERROR: STOP ERROR in pygp_functions.distance_given(). The landmass, 'India', 
     does not have a polygon at this timestamp (time=60).
     Please use land_begin to find the earliest timestamp with a polygon present. If land_begin returns an earlier timestamp, please check gplates file.

    
    julia> distance_given(test_df, 60, "Asia", "India")
    ERROR: STOP ERROR in pygp_functions.distance_given(). The landmass, 'India', 
     does not have a polygon at this timestamp (time=60).
     Please use land_begin to find the earliest timestamp with a polygon present. If land_begin returns an earlier timestamp, please check gplates file.
```
"""

# Hey Nick! A question here! mshould land1  land2 be input as strings when someone uses the fuction?
# i.e: distance_given(df, 40, "India", "Asia") ??


# or should I enter them as "land1", "land2" within this definition of the function?




function distance_given(df, time, land1, land2)
    
    q_first = @from i in test_df begin
        @where i.Reconstruction_Time_Ma == time && ((i.Land1 == land1 && i.Land2 == land2) || (i.Land2 == land1 && i.Land1 == land2))
        @select i.Closest_Distance_km
        @collect
    end

    q = @from i in q_first begin
        @where i.hasvalue == true
        @select i.value
        @collect
    end

    if length(q) == 0
        error("STOP ERROR in pygp_functions.distance_given().\n Pairing: '", land1, "' and '", land2, "' has not been found within given dataframe at this timestamp.\n Please check spelling of compared landmasses.\n If using function distance_given(), please check dataframe timestamps.\n \t if time intervals taken from original pygplates do not have your needed timestamp, please us function pygp_functions.distance_interp() \n If using function distance_interp(), please ensure time intervals within dataframe go PAST the requested timestamp. Thank you! ")
    end


    if q[1] == "NA"
        q1 = @from i in df begin
            @where i.Reconstruction_Time_Ma == time && i.Land1 == land1 && i.Land2 == land2
            @select i.Point1_Lat
            @collect
        end

        q2 = @from i in df begin
            @where i.Reconstruction_Time_Ma == time && i.Land1 == land1 && i.Land2 == land2
            @select i.Point2_Lat
            @collect
        end
        
        if length(q1) == 0
            q3 = @from i in df begin
                @where i.Reconstruction_Time_Ma == time && i.Land1 == land2 && i.Land2 == land1
                @select i.Point1_Lat
                @collect
            end

            q4 = @from i in df begin
                @where i.Reconstruction_Time_Ma == time && i.Land2 == land1 && i.Land1 == land2 
                @select i.Point2_Lat
                @collect
            end

            if q4[1]=="NA" && q3[1]!="NA"
                error("STOP ERROR in pygp_functions.distance_given(). The landmass, '", land1, "',\n does not have a polygon at this timestamp (time=", time, ").\n Please use land_begin to find the earliest timestamp with a polygon present. If land_begin returns an earlier timestamp, please check gplates file.")
            end # end error check

            if q3[1]=="NA" && q4[1]!="NA"
                error("STOP ERROR in pygp_functions.distance_given(). The landmass, '", land2, "',\n does not have a polygon at this timestamp (time=", time, ").\n Please use land_begin to find the earliest timestamp with a polygon present. If land_begin returns an earlier timestamp, please check gplates file.")
            end # end error check

            if q3[1]=="NA" && q4[1]=="NA"
                error("STOP ERROR in pygp_functions.distance_given(). Both landmasses, '", land1, "' and '", land2, "',\n do not have polygons at this timestamp\n (time=", time, ").\n Please use land_begin to find the earliest timestamp with polygons present. If land_begin returns an earlier timestamp, please check gplates file.")
            end
        end


        if q1[1]=="NA" && q2[1]!="NA"
            error("STOP ERROR in pygp_functions.distance_given(). The landmass, '", land1, "',\n does not have a polygon at this timestamp (time=", time, ").\n Please use land_begin to find the earliest timestamp with a polygon present. If land_begin returns an earlier timestamp, please check gplates file.")
        end # end error check

        if q2[1]=="NA" && q1[1]!="NA"
            error("STOP ERROR in pygp_functions.distance_given(). The landmass, '", land2, "',\n does not have a polygon at this timestamp (time=", time, ").\n Please use land_begin to find the earliest timestamp with a polygon present. If land_begin returns an earlier timestamp, please check gplates file.")
        end # end error check

        if q1[1]=="NA" && q2[1]=="NA"
            error("STOP ERROR in pygp_functions.distance_given(). Both landmasses, '", land1, "' and '", land2, "',\n do not have polygons at this timestamp \n (time=", time, ").\n Please use land_begin to find the earliest timestamp with polygons present. If land_begin returns an earlier timestamp, please check gplates file.")
        end

    end 


    distance_km = q[1]
    return distance_km

end


"""
    distance_interp(df, time, land1, land2)

Provides the distance at a given time point, to be used when timepoint requested is
not within the dataframe.

To pull up distance timestamps already recorded within the dataframe, a faster function 
would be distance_given

* `df` - dataframe created by pygplates output created by Wallis Bland
contains variables: 'Reconstruction_Time_Ma', 'Land1', 'Land2', 'Closest_Distance_Ma',
'Point1_Lat', 'Point1_Lon', 'Point2_Lat', 'Point2_Lon'

* `time` - timepoint for distance requested. May be any given time as long as the recorded
timespan within the dataframe is extends beyond requested timepoint.

* `land1` - first landmass chosen for distance comparison. Ensure this is within 'string' format

* `land2` - second landmass chosed for distance comparison. Ensure this is within 'string' format

NOTE: Lands 1 and 2 do not need to line up with Land1 or 2 within the dataframe in this function, as it will search both!

Process: Takes the two closest timestamps recorded in dataframe  produces a weighed average

# Examples
```julia-repl


julia> test_df = DataFrame(Reconstruction_Time_Ma=[0, 10, 20, 30, 40, 50, 60, 0, 10, 20, 30, 40, 50, 60], Land1=["India", "India", "India", "India", "India", "India", "India", "India", "India", "India", "India", "India", "India", "India"], Land2=["Asia", "Asia", "Asia", "Asia", "Asia", "Asia", "Asia", "Madagascar", "Madagascar", "Madagascar", "Madagascar", "Madagascar", "Madagascar","Madagascar"], Closest_Distance_km=[0, 1000, 2000, 3500, 4100, 4900, "NA", 4900, 4100, 3500, 2000, 1000, 0, "NA"], Point1_Lat=[1,2,3,4,5,6,"NA",1,2,3,4,5,6,"NA"], Point1_Lon=[1,2,3,4,5,6,"NA",1,2,3,4,5,6,"NA"], Point2_Lat=[1,2,3,4,5,6,7,1,2,3,4,5,6,"NA"], Point2_Lon=[1,2,3,4,5,6,7,1,2,3,4,5,6,"NA"])

14×8 DataFrame
 Row │ Reconstruction_Time_Ma  Land1   Land2       Closest_Distance_km  Point1_Lat  Point1_Lon  Point2_Lat  Point2_Lon 
     │ Int64                   String  String      Any                  Any         Any         Any         Any        
─────┼─────────────────────────────────────────────────────────────────────────────────────────────────────────────────
   1 │                      0  India   Asia        0                    1           1           1           1
   2 │                     10  India   Asia        1000                 2           2           2           2
   3 │                     20  India   Asia        2000                 3           3           3           3
   4 │                     30  India   Asia        3500                 4           4           4           4
   5 │                     40  India   Asia        4100                 5           5           5           5
   6 │                     50  India   Asia        4900                 6           6           6           6
   7 │                     60  India   Asia        NA                   NA          NA          7           7
   8 │                      0  India   Madagascar  4900                 1           1           1           1
   9 │                     10  India   Madagascar  4100                 2           2           2           2
  10 │                     20  India   Madagascar  3500                 3           3           3           3
  11 │                     30  India   Madagascar  2000                 4           4           4           4
  12 │                     40  India   Madagascar  1000                 5           5           5           5
  13 │                     50  India   Madagascar  0                    6           6           6           6
  14 │                     60  India   Madagascar  NA                   NA          NA          NA          NA

julia> distance_km = distance_interp(test_df, 44, "India", "Asia")

4420.0


# POTENTIAL ERRORS

* INCORRECT VARIABLES / MISSPELLINGS
    In this case the variable was either spelt wrong or for some other reason Variable was never printed within the given dataframe.

    julia> distance_interp(test_df, 44, "Asia", "Indai")
    ERROR: STOP ERROR in pygp_functions.distance_interp().
     Pairing: 'Asia' and 'Indai' has not been found within given dataframe at this timestamp.



* INCORRECT PAIRING
    Lands exist within the dataframe, but were never compared distance wise 
    If original pygplates code was run correctly, this SHOULD NOT OCCUR. More likely need to check spelling

    julia> distance_interp(test_df, 44, "Asia", "Madagascar")
    ERROR: STOP ERROR in pygp_functions.distance_interp().
     Pairing: 'Asia' and 'Madagascar' has not been found within given dataframe.



* TIME INVALID
    If the timestamp requested is beyond the reconstructed years (in example our further back year is 60, we requested 64)
    there is no higher time for the interpolator to use as a weighted average. 
    Return to pygplates python code to change the output years to cover the correct range.

    julia> distance_km = distance_interp(test_df, 64, "India", "Asia")
    ERROR: STOP ERROR in pygpfun_ctions.distance_interp().
     Timestamp: 64 is beyond the reconstructed time in the given dataframe



* NAs within Outputs
    These errors will appear as errors within the pygrplatesfunction.distance_given() code.
    When the selected higher and lower time are sent through the distance_given function to retrieve the distances for
    the weighted average, they will return an error stating which land has failed at which time stamp.

    This is most likely to occur if there are missing polygons within outputs, especially when polygons appear and disappear.
    Polygon location and 'closest distance' can be set to 0s before they appear. However that is up to the user to decide.

    julia> distance_interp(test_df, 56, "India", "Asia") 
    ERROR: STOP ERROR in pygp_functions.distance_given(). The landmass, 'India', 
     does not have a polygon at this timestamp (time=60).
     Please use land_begin to find the earliest timestamp with a polygon present. If land_begin returns an earlier timestamp, please check gplates file.
    

```
"""


# HEY SHOULD THIS SKIP ONES THAT HAVE NA DISTANCES?

function distance_interp(df, time, land1, land2)
    
    q = @from i in df begin
        @where ((i.Land1 == land1 && i.Land2 == land2) || (i.Land2 == land1 && i.Land1 == land2))
        @select i.Reconstruction_Time_Ma
        @collect DataFrame
    end

    List_of_times = @from i in q begin
        @where i.hasvalue == true
        @select i.value
        @collect
    end

    # should this actually be
    """
    List_of_times = @from i in df begin
        @where i.Closest_Distance_km != "NA" && ((i.Land1 == land1 && i.Land2 == land2) || (i.Land2 == land1 && i.Land1 == land2))
        @select i.Reconstruction_Time_Ma
        @collect
    end

    """

    if (length(List_of_times) == 0)
        error("STOP ERROR in pygp_functions.distance_interp().\n Pairing '", land1, "' and '", land2, "' has not been found within given dataframe")
    end

    if time - last(List_of_times) > 0
        error("STOP ERROR in pygp_functions.distance_interp().\n Timestamp: ", time, " (mya) is beyond the reconstructed time in the given dataframe")
    end

    i = 1
    global time_low, time_high
    for timecheck in List_of_times

        if timecheck == "NA"
            continue
        end

        if (time - timecheck <= 0)
            time_low = List_of_times[i-1]
            time_high = List_of_times[i]

            break
        end
        
        i = i+1
    end
    
    dist_high = distance_given(df, time_high, land1, land2)
    dist_low = distance_given(df, time_low, land1, land2)
    
    t2b_high = abs(time - time_high)
    t2b_low = abs(time - time_low)
    t2bsum = t2b_high + t2b_low
    
    weight_high = (t2bsum - t2b_high)/t2bsum
    weight_low = (t2bsum - t2b_low)/t2bsum
    
    weight_dist_high = weight_high * dist_high
    weight_dist_low = weight_low * dist_low
    
    distance_km =  weight_dist_high + weight_dist_low
    return distance_km

end

"""
    oldest_polygon(df, land)

Provides the last timestamp in which a polygon can be found. Some gpml files have polygons that disappear and reappear!

* `df` - dataframe created by pygplates output created by Wallis Bland
contains variables: 'Reconstruction_Time_Ma', 'Land1', 'Land2', 'Closest_Distance_Ma',
'Point1_Lat', 'Point1_Lon', 'Point2_Lat', 'Point2_Lon'

* `land` - name of landmass. Ensure this is within 'string' format


# Examples
```julia-repl


julia> test_df = DataFrame(Reconstruction_Time_Ma=[0, 10, 20, 30, 40, 50, 60, 0, 10, 20, 30, 40, 50, 60], Land1=["India", "India", "India", "India", "India", "India", "India", "India", "India", "India", "India", "India", "India", "India"], Land2=["Asia", "Asia", "Asia", "Asia", "Asia", "Asia", "Asia", "Madagascar", "Madagascar", "Madagascar", "Madagascar", "Madagascar", "Madagascar","Madagascar"], Closest_Distance_km=[0, 1000, 2000, 3500, 4100, 4900, "NA", 4900, 4100, 3500, 2000, 1000, 0, "NA"], Point1_Lat=[1,2,3,4,5,6,"NA",1,2,3,4,5,6,"NA"], Point1_Lon=[1,2,3,4,5,6,"NA",1,2,3,4,5,6,"NA"], Point2_Lat=[1,2,3,4,5,6,7,1,2,3,4,5,6,"NA"], Point2_Lon=[1,2,3,4,5,6,7,1,2,3,4,5,6,"NA"])

14×8 DataFrame
 Row │ Reconstruction_Time_Ma  Land1   Land2       Closest_Distance_km  Point1_Lat  Point1_Lon  Point2_Lat  Point2_Lon 
     │ Int64                   String  String      Any                  Any         Any         Any         Any        
─────┼─────────────────────────────────────────────────────────────────────────────────────────────────────────────────
   1 │                      0  India   Asia        0                    1           1           1           1
   2 │                     10  India   Asia        1000                 2           2           2           2
   3 │                     20  India   Asia        2000                 3           3           3           3
   4 │                     30  India   Asia        3500                 4           4           4           4
   5 │                     40  India   Asia        4100                 5           5           5           5
   6 │                     50  India   Asia        4900                 6           6           6           6
   7 │                     60  India   Asia        NA                   NA          NA          7           7
   8 │                      0  India   Madagascar  4900                 1           1           1           1
   9 │                     10  India   Madagascar  4100                 2           2           2           2
  10 │                     20  India   Madagascar  3500                 3           3           3           3
  11 │                     30  India   Madagascar  2000                 4           4           4           4
  12 │                     40  India   Madagascar  1000                 5           5           5           5
  13 │                     50  India   Madagascar  0                    6           6           6           6
  14 │                     60  India   Madagascar  NA                   NA          NA          NA          NA

  julia> oldest_polygon(test_df, "Asia")
  60

  julia> oldest_polygon(test_df, "India")
  50

  julia> oldest_polygon(test_df, "Indai")
  ERROR: STOP ERROR in pygp_functions.oldest_polygon.
   Landmass 'Indai' either does not exist in dataframe or does not contain any polygons.

```
"""


function oldest_polygon(df, land)
    
    q = @from i in df begin
        @where i.Point1_Lat != "NA" && i.Land1 == land
        @select i.Reconstruction_Time_Ma
        @collect DataFrame
    end

    List_of_times = @from i in q begin
        @where i.hasvalue == true
        @select i.value
        @collect
    end
    
    if length(List_of_times) == 0
        List_of_times = @from i in df begin
            @where i.Point2_Lat != "NA" && i.Land2 == land
            @select i.Reconstruction_Time_Ma
            @collect
        end
    end

    if length(List_of_times) == 0
        error("STOP ERROR in pygp_functions.oldest_polygon.\n Landmass '", land, "' either does not exist in dataframe or does not contain any polygons.")
    end

    oldest_time_stamp_with_polygon = last(List_of_times)
    return oldest_time_stamp_with_polygon

end

"""
    land_last_touch(df, land1, land2)

Provides the most recent timestamp in which the landmasses have touched according to your dataframe

* `df` - dataframe created by pygplates output created by Wallis Bland
contains variables: 'Reconstruction_Time_Ma', 'Land1', 'Land2', 'Closest_Distance_Ma',
'Point1_Lat', 'Point1_Lon', 'Point2_Lat', 'Point2_Lon'

* `land1` - first landmass chosen for distance comparison. Ensure this is within 'string' format

* `land2` - second landmass chosed for distance comparison. Ensure this is within 'string' format

NOTE: Lands 1 and 2 do not need to line up with Land1 or 2 within the dataframe in this function, as it will search both!


# Examples
```julia-repl


julia> test_df = DataFrame(Reconstruction_Time_Ma=[0, 10, 20, 30, 40, 50, 60, 0, 10, 20, 30, 40, 50, 60], Land1=["India", "India", "India", "India", "India", "India", "India", "India", "India", "India", "India", "India", "India", "India"], Land2=["Asia", "Asia", "Asia", "Asia", "Asia", "Asia", "Asia", "Madagascar", "Madagascar", "Madagascar", "Madagascar", "Madagascar", "Madagascar","Madagascar"], Closest_Distance_km=[0, 1000, 2000, 3500, 4100, 4900, "NA", 4900, 4100, 3500, 2000, 1000, 0, "NA"], Point1_Lat=[1,2,3,4,5,6,"NA",1,2,3,4,5,6,"NA"], Point1_Lon=[1,2,3,4,5,6,"NA",1,2,3,4,5,6,"NA"], Point2_Lat=[1,2,3,4,5,6,7,1,2,3,4,5,6,"NA"], Point2_Lon=[1,2,3,4,5,6,7,1,2,3,4,5,6,"NA"])

14×8 DataFrame
 Row │ Reconstruction_Time_Ma  Land1   Land2       Closest_Distance_km  Point1_Lat  Point1_Lon  Point2_Lat  Point2_Lon 
     │ Int64                   String  String      Any                  Any         Any         Any         Any        
─────┼─────────────────────────────────────────────────────────────────────────────────────────────────────────────────
   1 │                      0  India   Asia        0                    1           1           1           1
   2 │                     10  India   Asia        1000                 2           2           2           2
   3 │                     20  India   Asia        2000                 3           3           3           3
   4 │                     30  India   Asia        3500                 4           4           4           4
   5 │                     40  India   Asia        4100                 5           5           5           5
   6 │                     50  India   Asia        4900                 6           6           6           6
   7 │                     60  India   Asia        NA                   NA          NA          7           7
   8 │                      0  India   Madagascar  4900                 1           1           1           1
   9 │                     10  India   Madagascar  4100                 2           2           2           2
  10 │                     20  India   Madagascar  3500                 3           3           3           3
  11 │                     30  India   Madagascar  2000                 4           4           4           4
  12 │                     40  India   Madagascar  1000                 5           5           5           5
  13 │                     50  India   Madagascar  0                    6           6           6           6
  14 │                     60  India   Madagascar  NA                   NA          NA          NA          NA



```
"""


function land_last_touch(df, land1, land2)
    distance = 0    
    
    q = @from i in df begin
        @where i.Closest_Distance_km == 0 && ((i.Land1 == land1 && i.Land2 == land2) || (i.Land2 == land1 && i.Land1 == land2))
        @select i.Reconstruction_Time_Ma
        @collect DataFrame
    end
    List_of_times = @from i in q begin
        @where i.hasvalue == true
        @select i.value
        @collect
    end
    
    if List_of_times == 0
        times_by_lands = @from i in df begin
            @where (i.Land1 == land1 && i.Land2 == land2) || (i.Land2 == land1 && i.Land1 == land2)
            @select i.Reconstruction_Time_Ma
            @collect
        end

        if length(times_by_lands) == 0
            error("STOP ERROR in pygp_functions.distance_interp().\n Pairing '", land1, "' and '", land2, "' has not been found within given dataframe.")
        else
            error("STOP ERROR in pygp_functions.distance_interp().\n Pairing '", land1, "' and '", land2, "' do not touch, based on the time intervals reconstructed by gplates in your dataframe.")
        end 
    end

    Last = List_of_times[1]
    return Last
end




"""
    pygptxt_read(file)

Function to read in pygplates txt files (output from previous pygplates code by wallis bland). For csv files, please use pygpcsv_read()

* `file` - location of txt file. 
Contains dataframe created by pygplates output created by Wallis Bland
with variables: 'Reconstruction_Time_Ma', 'Land1', 'Land2', 'Closest_Distance_Ma',
'Point1_Lat', 'Point1_Lon', 'Point2_Lat', 'Point2_Lon'

# Examples
```julia-repl

pygplates_txt("Desktop/output.txt")

14×8 DataFrame
 Row │ Reconstruction_Time_Ma  Land1  Land2       Closest_Distance_km  Point1_Lat  Point1_Lon  Point2_Lat  Point2_Lon 
     │ Any                     Any    Any         Any                  Any         Any         Any         Any        
─────┼────────────────────────────────────────────────────────────────────────────────────────────────────────────────
   1 │ 0.0                     India  Asia        0.0                  1.0         1.0         1.0         1.0
   2 │ 10.0                    India  Asia        1000.0               2.0         2.0         2.0         2.0
   3 │ 20.0                    India  Asia        2000.0               3.0         3.0         3.0         3.0
   4 │ 30.0                    India  Asia        3500.0               4.0         4.0         4.0         4.0
   5 │ 40.0                    India  Asia        4100.0               5.0         5.0         5.0         5.0
   6 │ 50.0                    India  Asia        4900.0               6.0         6.0         6.0         6.0
   7 │ 60.0                    India  Asia        NA                   NA          NA          7.0         7.0
   8 │ 0.0                     India  Madagascar  4900.0               1.0         1.0         1.0         1.0
   9 │ 10.0                    India  Madagascar  4100.0               2.0         2.0         2.0         2.0
  10 │ 20.0                    India  Madagascar  3500.0               3.0         3.0         3.0         3.0
  11 │ 30.0                    India  Madagascar  2000.0               4.0         4.0         4.0         4.0
  12 │ 40.0                    India  Madagascar  1000.0               5.0         5.0         5.0         5.0
  13 │ 50.0                    India  Madagascar  0.0                  6.0         6.0         6.0         6.0
  14 │ 60.0                    India  Madagascar  NA                   NA          NA          NA          NA



```


"""


function pygptxt_read(file)

    fhandle = open(file)

    # Read the lines to an array
    lines = readlines(fhandle)
    close(fhandle)

    if length(lines) == 0
        error("Oops! looks like file: '", file, "' is empty! Please ensure this is the correct file, or that the pygplates python dataframe is correctly printing to file.")
    end

    # lines = readlines("Desktop/School/Rhacophoridae/Gplates/output.txt")

    # Go through each line and split by \t
    for i in 1:length(lines)
        line = lines[i]

        if length(line) == 0 # if line is blank. Skip
            continue
        end

        if i == 1
            parts = split(line, "\t") #split the column headings out! I dont think we actually use these? I tried to below, but can't get that to work..

            # Attempting to make the headers into vectors, so that this read TXT can be used anywhere!
            # HOW?

            #for i in length(parts)
            #    headers = Any[]
            #    push!(headers, Symbol(parts[1]))
            #    headers[1] = Any[]
            #end

            # create the vectors just the once. global allowsit to be accessed later outside of the loop.
            global Reconstruction_Time_Ma = Any[]
            global Land1 = Any[]
            global Land2 = Any[]
            global Closest_Distance_km = Any[]
            global Point1_Lat = Any[]
            global Point1_Lon = Any[]
            global Point2_Lat = Any[]
            global Point2_Lon = Any[]

        else
            datam = split(line, "\t") # split each line with data by \t 
            
            # input each into the corresponding array (bin). All info is in string format that the moment.
            push!(Reconstruction_Time_Ma, datam[1])
            push!(Land1, datam[2])
            push!(Land2, datam[3])
            push!(Closest_Distance_km, datam[4])
            push!(Point1_Lat, datam[5])
            push!(Point1_Lon, datam[6])
            push!(Point2_Lat, datam[7])
            push!(Point2_Lon, datam[8])
        end
    end


    # NOW have it go through CDK, P1Lat, P2Lat, P1Long and P2Lon to change non Nas into Floats!
    # this is the long hand version, shorthand version actually used: 

    """
    a = 1
    for line in Reconstruction_Time_Ma
        if line != "NA"
            Reconstruction_Time_Ma[a] = parse(Float64, line)
        end
        a = a + 1
    end
    a = 1
    for line in Closest_Distance_km
        if line != "NA"
            Closest_Distance_km[a] = parse(Float64, line)
        end
        a = a + 1
    end
    a = 1
    for line in Point1_Lat
        if line != "NA"
            Point1_Lat[a] = parse(Float64, line)
        end
        a = a + 1
    end
    a = 1
    for line in Point1_Lon
        if line != "NA"
            Point1_Lon[a] = parse(Float64, line)
        end
        a = a + 1
    end
    a = 1
    for line in Point2_Lat
        if line != "NA"
            Point2_Lat[a] = parse(Float64, line)
        end
        a = a + 1
    end
    a = 1
    for line in Point2_Lon
        if line != "NA"
            Point2_Lon[a] = parse(Float64, line)
        end
        a = a + 1
    end

    df_2 = DataFrame(Reconstruction_Time_Ma = Reconstruction_Time_Ma, Land1 = Land1, Land2 = Land2, Closest_Distance_km = Closest_Distance_km, Point1_Lat = Point1_Lat, Point1_Lon = Point1_Lon, Point2_Lat = Point2_Lat, Point2_Lon = Point2_Lon)
    return df
    """

    # Shorthand version has us loop through them in an array of arrays
    df_array = Any[Reconstruction_Time_Ma, Land1, Land2, Closest_Distance_km, Point1_Lat, Point1_Lon, Point2_Lat, Point2_Lon]

    groups_to_change =  [1,4,5,6,7,8] # we know which groups need to be changed to floats

    for group in df_array[groups_to_change] # loop through array of array
        a = 1 # start at line one of data again
        for line in group
            if line != "NA" # if data is anything but NA, change it to a float (only other string info should be in the LANDS, all info in these columns should be actual numbers!)
                group[a] = parse(Float64, line) 
            end 
            a = a + 1
        end
    end

    df = DataFrame(df_array, [:Reconstruction_Time_Ma, :Land1, :Land2, :Closest_Distance_km, :Point1_Lat, :Point1_Lon, :Point2_Lat, :Point2_Lon])
    return df
end


"""
    pygpcsv_read(file)

Function to read in pygplates csv files (output from previous pygplates code by wallis bland).

* `file` - location of txt file. 
Contains dataframe created by pygplates output created by Wallis Bland
with variables: 'Reconstruction_Time_Ma', 'Land1', 'Land2', 'Closest_Distance_Ma',
'Point1_Lat', 'Point1_Lon', 'Point2_Lat', 'Point2_Lon'

# Examples
```julia-repl

julia> pygpcsv_read("/Users/wbla447/Desktop/Files/Gplates/Nickgraphicoutput.csv")

126×8 DataFrame
 Row │ Reconstruction_Time_Ma  Land1        Land2          Closest_Distance_km ⋯
     │ Int64                   String15     String15       String31            ⋯
─────┼──────────────────────────────────────────────────────────────────────────
   1 │                      0  New Zealand  Antarctica     1347.0694100138985  ⋯
   2 │                      5  New Zealand  Antarctica     1306.5091361116174
   3 │                     10  New Zealand  Antarctica     1269.8588986531095
   4 │                     15  New Zealand  Antarctica     1276.630907584367
   5 │                     20  New Zealand  Antarctica     1271.728477569835   ⋯
   6 │                     25  New Zealand  Antarctica     1254.9240899936028
   7 │                     30  New Zealand  Antarctica     1212.368755489743
   8 │                     35  New Zealand  Antarctica     1162.6193444486703
  ⋮  │           ⋮                  ⋮             ⋮                 ⋮          ⋱
 120 │                     70  Australia    New Caledonia  NA                  ⋯
 121 │                     75  Australia    New Caledonia  NA
 122 │                     80  Australia    New Caledonia  NA
 123 │                     85  Australia    New Caledonia  NA
 124 │                     90  Australia    New Caledonia  NA                  ⋯
 125 │                     95  Australia    New Caledonia  NA
 126 │                    100  Australia    New Caledonia  NA
                                                  4 columns and 111 rows omitted



```


"""

function pygpcsv_read(file)

    csv_reader = CSV.read(file, DataFrame)
    df = csv_reader[!, 2:9]
    return df

end

end # ENDING Parsers