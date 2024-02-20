
# Run this before starting python
# Ensure that it directs to the location of your pygplates download
export PYTHONPATH=$PYTHONPATH:/Users/wbla447/Desktop/School/Gplates/pygplates_rev28_python38_MacOS64



'''
Dear User:

Only edit 'INPUT' and 'OUTPUT' sections

Thank you
Wallis :)
'''


# Input the missing frames as NA distances, so we know it at least went through those years and wasn't a mistake!

python3  

import numpy as np
import pandas as pd
import pygplates
pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)


# INPUT #

# change your information from here to years!
Number_of_Lands = 3

rotation_model = pygplates.RotationModel('/Users/wbla447/Desktop/Files/Gplates/gplates_frogshapes/Matthews_etal_GPC_2016_410-0Ma_GK07.rot')
# Load the features.
# If this shape file has more than 1 Plate ID, then you theoretically should get multiple results? testing first with India/Africa!
# Going to have to figure out naming section!

features1 = pygplates.FeatureCollection('/Users/wbla447/Desktop/Files/Gplates/gplates_frogshapes/lm_India.gpml')
features2 = pygplates.FeatureCollection('/Users/wbla447/Desktop/Files/Gplates/gplates_frogshapes/lm_Africa.gpml')
features3 = pygplates.FeatureCollection('/Users/wbla447/Desktop/Files/Gplates/gplates_frogshapes/lm_Madagascar.gpml')

Features = [features1, features2, features3] #put all the features in array to access later
LandNames = ['India', 'Africa', 'Madagascar'] #put all the names in an array to access later
# Make sure the features and their names are lines up! Feel free to add more feature files, just make sure you change the Number_of_Lands variable above to reflect this!

Years = [] # leave blank
for y in range(101): # range needs to be one above your max year
    if y%20 == 0: # every million years for 100 million years
        Years.append(y)



print(Years)


###### Dont touch below! ####

# CALCULATE #

LandsMatrixSize = Number_of_Lands*Number_of_Lands

Lands = np.arange(0,LandsMatrixSize)
Lands_df = pd.DataFrame(Lands.reshape(-1, Number_of_Lands)) #creates a dataframe so we can later access these numbers (is there an easier way to do this than we did? probably...)
print(Lands_df) # make sure it's correct.

def getIndexes(dfObj, value): 
      
    # Empty list 
    listOfPos = [] 
      
    # isin() method will return a dataframe with  
    # boolean values, True at the positions     
    # where element exists 
    result = dfObj.isin([value]) 
      
    # any() method will return  
    # a boolean series 
    seriesObj = result.any() 
  
    # Get list of column names where  
    # element exists 
    columnNames = list(seriesObj[seriesObj == True].index) 
     
    # Iterate over the list of columns and 
    # extract the row index where element exists 
    for col in columnNames: 
        rows = list(result[col][result[col] == True].index) 
  
        for row in rows: 
            listOfPos.append(row)
            listOfPos.append(col)
              
    # This list contains a list tuples with  
    # the index of element in the dataframe 
    return listOfPos

Land1 = []
Land2 = []
Land1Lat = []
Land1Lon = []
Land2Lat = []
Land2Lon = []
Reconstruction_Time_Ma = []
Closest_Distance = []

for l in range(LandsMatrixSize): 
    # by iterating over this range, we will pull the same numbers that are in our matrix
    # getIndexes will then use this iterated number, finding it's col/row within the matrix, which will point at which features file and land name will be pulled by pygplates
     
    rowcol = getIndexes(Lands_df, l) # will result in a size2 array with (row, col)
    row = rowcol[0] # separate out the rows and columns
    col = rowcol[1]
    
    if row < col: #note if using row != col then you end up with x2 of everything as it will do both sets [a,b] and [b,a]
        for i in Years: # matrix will be sorted by set compared and THEN time (just looks a lot cleaner)
            reconstruction_time = i
            distance_reconstructed = 50000 #reset distance reconstructed each round (NOTE: earth is 40000 km, when I was setting this to 0, is would messup when the actual reconstructed distance was 0 if there were more than 1 polygones reconstructed!) 
            reconstructed_feature_geometries1 = []
            reconstructed_feature_geometries2 = []
            pygplates.reconstruct(Features[row], rotation_model, reconstructed_feature_geometries1, reconstruction_time) 
            pygplates.reconstruct(Features[col], rotation_model, reconstructed_feature_geometries2, reconstruction_time)
             
            for reconstructed_feature_geometry1 in reconstructed_feature_geometries1:
                for reconstructed_feature_geometry2 in reconstructed_feature_geometries2:
                   
                    if reconstruction_time == 0: # if at the current time use 'get_present_day_geometry'
                        distance_reconstructed_next, Land1_Location_Next, Land2_Location_Next = pygplates.GeometryOnSphere.distance(
                            reconstructed_feature_geometry1.get_present_day_geometry(), reconstructed_feature_geometry2.get_present_day_geometry(), return_closest_positions = True) # Will return a tuple giving us distance and locations)
                         
                        if distance_reconstructed == 50000:  distance_reconstructed = distance_reconstructed_next # if first round, set distance_rec.
                        elif distance_reconstructed_next < distance_reconstructed: # otherwise if the new reconstruction is smaller, take that one
                            distance_reconstructed = distance_reconstructed_next
                            Land1_Location_Save = Land1_Location_Next # If saving that distance, also save the points used
                            Land2_Location_Save = Land2_Location_Next
                       
                    else: # all besides current time
                        distance_reconstructed_next, Land1_Location_Next, Land2_Location_Next = pygplates.GeometryOnSphere.distance(
                            reconstructed_feature_geometry1.get_reconstructed_geometry(), reconstructed_feature_geometry2.get_reconstructed_geometry(), return_closest_positions = True)    
                            
                        if distance_reconstructed == 50000: distance_reconstructed =  distance_reconstructed_next
                        elif distance_reconstructed_next < distance_reconstructed: #If distance is shorter than whats saved so far, use the new distance
                            distance_reconstructed = distance_reconstructed_next
                            Land1_Location_Save = Land1_Location_Next # If using the new distance, also save that location
                            Land2_Location_Save = Land2_Location_Next
             
            distance_reconstructed_in_kms = distance_reconstructed * pygplates.Earth.mean_radius_in_kms # turn it from pygplates automatic measurment to Kms
              
            # Then add it all to our dataframe!
            Land1.append(LandNames[row]) # Lined up with reconstruction 1 (This is why the LandNames and Features arrays need to be in the same order!)
            Land2.append(LandNames[col]) # Lined up with Reconstruction 2
            Reconstruction_Time_Ma.append(i)
             
            if len(reconstructed_feature_geometries1) == 0:
                Closest_Distance.append('NA')
                Land1Lat.append('NA')
                Land1Lon.append('NA')
                Land2Lat.append(Land2_Location_Save.to_lat_lon_list()[0][0])
                Land2Lon.append(Land2_Location_Save.to_lat_lon_list()[0][1])
            elif len(reconstructed_feature_geometries2) == 0:  
                Closest_Distance.append('NA')
                Land1Lat.append(Land1_Location_Save.to_lat_lon_list()[0][0])
                Land1Lon.append(Land1_Location_Save.to_lat_lon_list()[0][1])
                Land2Lat.append('NA')
                Land2Lon.append('NA')
            else:   
                Closest_Distance.append(distance_reconstructed_in_kms)
                Land1Lat.append(Land1_Location_Save.to_lat_lon_list()[0][0])
                Land1Lon.append(Land1_Location_Save.to_lat_lon_list()[0][1])
                Land2Lat.append(Land2_Location_Save.to_lat_lon_list()[0][0])
                Land2Lon.append(Land2_Location_Save.to_lat_lon_list()[0][1])

len(Land1)
len(Land2)
len(Land1Lat)
len(Land1Lon)
len(Land2Lat)
len(Land2Lon)
len(Reconstruction_Time_Ma)
len(Closest_Distance)

Data = {'Reconstruction_Time_Ma': Reconstruction_Time_Ma, 'Land1': Land1, 'Land2': Land2, 'Closest_Distance_km': Closest_Distance, 'Point1_Lat': Land1Lat, 'Point1_Lon': Land1Lon, 'Point2_Lat': Land2Lat, 'Point2_Lon': Land2Lon}
Reconstruction_df = pd.DataFrame(Data, columns = ['Reconstruction_Time_Ma','Land1','Land2','Closest_Distance_km','Point1_Lat', 'Point1_Lon', 'Point2_Lat', 'Point2_Lon'])

print(Reconstruction_df)

header = ["Reconstruction_Time_Ma","Land1","Land2","Closest_Distance_km","Point1_Lat", "Point1_Lon", "Point2_Lat", "Point2_Lon"] # Heads match parser



# OUTPUT # 

# CSV version
Reconstruction_df.to_csv('Desktop/School/Rhacophoridae/Gplates/output.csv') # Change to what you want your file name to be

# tab delimited txt version
Reconstruction_df.to_csv('Desktop/School/Rhacophoridae/Gplates/output.txt', sep='\t', index=False) # Change to what you want your file name to be

"""
The above cuts them apart at the halfway like below

    Reconstruction_Time_Ma   Land1       Land2 Closest_Distance (km)  \
0                        0   India      Africa               1995.38   
1                       20   India      Africa               2131.67   
2                       40   India      Africa               2620.13   
3                       60   India      Africa               1475.88   
4                       80   India      Africa               1258.68   
5                      100   India      Africa               1022.91   
6                      120   India      Africa               1077.51   
7                      140   India      Africa               894.089   
8                      160   India      Africa               33.1628   
9                      180   India      Africa                     0   
10                     200   India      Africa                     0   
11                       0   India  Madagascar               3517.07   
12                      20   India  Madagascar               3150.88   
13                      40   India  Madagascar               2611.85   
14                      60   India  Madagascar               1543.92   
15                      80   India  Madagascar               24.3555   
16                     100   India  Madagascar                     0   
17                     120   India  Madagascar                     0   
18                     140   India  Madagascar                     0   
19                     160   India  Madagascar                     0   
20                     180   India  Madagascar                    NA   
21                     200   India  Madagascar                    NA   
22                       0  Africa  Madagascar               296.155   
23                      20  Africa  Madagascar               463.959   
24                      40  Africa  Madagascar               567.415   
25                      60  Africa  Madagascar               860.101   
26                      80  Africa  Madagascar               860.101   
27                     100  Africa  Madagascar               678.315   
28                     120  Africa  Madagascar               806.702   
29                     140  Africa  Madagascar               894.089   
30                     160  Africa  Madagascar               29.3166   
31                     180  Africa  Madagascar                    NA   
32                     200  Africa  Madagascar                    NA   

    Point1 (Lat)  Point1 (Lon) Point2 (Lat) Point2 (Lon)  
0      23.580400     66.354400        12.66      51.3376  
1      13.611693     65.573805      6.61886      47.4276  
2       9.446947     62.175879      3.02778      39.3545  
3      -9.503571     50.028995     -4.52081      37.6297  
4     -29.954622     38.162792     -28.6459      25.2626  
5     -40.694750     34.311966     -36.1945      24.0548  
6     -42.378000     29.989759     -37.1118      19.3927  
7     -35.506867     33.755651     -28.1553      29.9161  
8     -44.150288     19.414394     -44.1802      19.0007  
9     -32.859359     14.344113     -32.8594      14.3441  
10    -33.845713     14.504332     -33.8457      14.5043  
11    -33.845713     14.504332     -33.8457      14.5043  
12    -33.845713     14.504332     -33.8457      14.5043  
13    -33.845713     14.504332     -33.8457      14.5043  
14    -14.081942     50.887945     -23.1818      39.8034  
15    -29.954622     38.162792     -29.9034      37.9171  
16    -29.954622     38.162792     -29.9034      37.9171  
17    -29.954622     38.162792     -29.9034      37.9171  
18    -29.954622     38.162792     -29.9034      37.9171  
19    -29.954622     38.162792     -29.9034      37.9171  
20    -29.954622     38.162792           NA           NA  
21    -29.954622     38.162792           NA           NA  
22    -15.857849     40.579859     -16.9094      43.1305  
23    -19.498249     38.359863     -20.9882      42.5141  
24    -22.333590     33.594976     -24.2539      38.7429  
25    -22.333590     33.594976     -24.2539      38.7429  
26    -22.333590     33.594976     -24.2539      38.7429  
27    -33.939812     24.830372     -35.0467      32.1108  
28    -37.583760     18.940404     -42.3784      26.0542  
29    -28.155349     29.916144     -35.5069      33.7557  
30    -43.967127     18.971672     -43.9409      19.3361  
31    -43.967127     18.971672           NA           NA  
32    -43.967127     18.971672           NA           NA  

"""


# build the dataframe without using pandas

f = open("Desktop/School/Rhacophoridae/Gplates/handwriteoutput.txt", "a") # Change to what you want your file name to be
print("\t".join(header) + "\n" + \
    "\n".join(str(Reconstruction_Time_Ma[i]) + "\t" + \
        Land1[i] + "\t" + \
            Land2[i] + "\t" + \
                str(Closest_Distance[i]) + "\t" + \
                    str(Land1Lat[i]) + "\t" + \
                        str(Land1Lon[i]) + "\t" + \
                            str(Land2Lat[i]) + "\t" + \
                                str(Land2Lon[i]) + "\t" \
                                    for i in range(len(Land1))), file=f)
f.close()

