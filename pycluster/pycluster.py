#!/usr/bin/python
"""
This file is part of Pycluster API.

Pycluster is a free software API for automating keyword searches in 
Google's Places API and finding and displaying geographic cluster
information using the K-means clusterization technique. With Pycluster
API any user can search for specific keywords on Google Maps for a given
region and have the results saved into a .csv file. This file
can later be used for geographic clusterization using K-means.

Pycluster is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

Pycluster is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with Pycluster.  If not, see <https://www.gnu.org/licenses/>.
"""

import pandas as pd, numpy as np
import folium as fl
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from sklearn.cluster import KMeans
import requests
import json
import time
import os

PYCLUSTER_COLORS = ['maroon','blue','green','orange','red','purple', 'violet', 'brown', 'teal', 'peru','navy','springgreen','slategray','chocolate','darkolivegreen']

def searchPlacesAPI(coordinates_list, keywords_list,
                    radius_m, output_file_name, 
                    API_KEY=None):
    '''
    This function searches Places API on a given region
    defined by a central pair of coordinates, a radius and
    a set of keywords. The function saves a .csv file containing
    the search result points.
    @param:
        coordinates_list: list of [longitude, latitude] pairs
                          there are used geographic center for
                          each keyword provided;
        keywords_list:    list of keywords to be used in the 
                          search;
        radius_m:         search radius in meters;
        output_file_name: output .csv file name;
        API_KEY:          Google API key;
    '''

    if API_KEY is None:
        print("Please, provide a valid Google API Key for Places API search!")
        return -1

    final_data = []
    search_cnt=0
    for coordinate in coordinates_list:
        search_cnt = search_cnt+1
        for keyword in keywords_list:
            url = 'https://maps.googleapis.com/maps/api/place/nearbysearch/json?location=\
                '+coordinate+'&radius='+str(radius_m)+'&keyword='+str(keyword)+'&key='+str(API_KEY)
            while True:
                respon = requests.get(url)
                jj = json.loads(respon.text)
                results = jj['results']
                res_cnt=0
                for result in results:
                    name = result['name']
                    place_id = result ['place_id']
                    lat = result['geometry']['location']['lat']
                    lng = result['geometry']['location']['lng']
                    rating = result['rating']
                    types = result['types'][0]
                    vicinity = result['vicinity']
                    data = [name, place_id, lat, lng, rating, types, vicinity]
                    final_data.append(data)
                    res_cnt=res_cnt+1
                time.sleep(1)
                if 'next_page_token' not in jj:
                    break
                else:
                    next_page_token = jj['next_page_token']
                url = 'https://maps.googleapis.com/maps/api/place/nearbysearch/json?key='+str(api_key)+'&pagetoken='+str(next_page_token)

    labels = ['Name','ID', 'Latitude', 'Longitude', 'Rating', 'Types', 'Vicinity']
    df = pd.DataFrame.from_records(final_data, columns=labels)
    df = df.drop_duplicates()
    df.to_csv(output_file_name)

def findClusters(coords, n_clusters=3):
    ''' findClusters
        Groups a set of coordinates into 'n_clusters' clusters,
        also calculating the median and index of each group.
        Inputs:
            coords:     array of input coordinates
            n_clusters: desired number of clusters (default = 3)
        Output:
            [clusters, clindexes]
            clusters:   set of clusters found
            clindexes:  index of clusters found
    '''
    # adjust the coordinates in 'n_clusters' groups (clusters)
    km_coords = KMeans(n_clusters=n_clusters).fit(coords)
    clusters = km_coords.cluster_centers_ # center of each cluster
    clindexes = km_coords.predict(coords) # index of each cluster
    return clusters, clindexes

def plotClusterCharts(coords, coords_mean,
                      ellipse_mean, clusters, clindexes, 
                      n_clusters=3, regionName=None, saveDir=None):
    ''' plotClusterCharts()
        Plot centroides and mean in 2D graphics.
        Each cluster is plotted with a different color
        Inputs:
            coords:         coordinate array (lat, longi)
            coords_mean:    coordinates mean
            ellipse_mean:   ellipse mean
            clusters:       clusters
            clindexes:      clusters index (vector)
        Outputs:
            2D colored plot of clusters
    '''
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(20,10))
    ax1.scatter(coords[:, 0], coords[:, 1])
    ax1.scatter(coords_mean[0], coords_mean[1], s=200, color = 'purple')
    ax1.add_patch(ellipse_mean)
    ax1.set(title='MEAN AND STANDARD DEVIATION', xlabel='LONGITUDE',ylabel='LATITUDE')
    ax1.text(coords_mean[0], coords_mean[1], 'MEAN')

    # selects each cluster's color
    for z in range(0,n_clusters):
        dataZ = np.where(clindexes==z)
        newDataX = []
        newDataY = []
        for k in dataZ[0]:
            newDataX.append(coords[k][0])
            newDataY.append(coords[k][1])
        color = PYCLUSTER_COLORS[z] # gets the color for that point
        ax2.scatter(newDataX, newDataY, color=color)

    ax2.scatter(clusters[:,0], clusters[:,1], s=500, color = 'black')
    ax2.set(title='CLUSTERS - nclusters = %d'%(n_clusters), xlabel='LONGITUDE', ylabel='LATITUDE')

    ax1.ticklabel_format(useOffset=False)
    ax2.ticklabel_format(useOffset=False)
    plt.ticklabel_format(style='plain')
    if regionName is not None:
        if saveDir is not None:
            plt.savefig(saveDir+"/"+regionName+".png")
        else:
            plt.savefig(regionName+".png")

def printClusterCoords(clusters, rName):
    ''' printClusterCoords
        Prints the central coordinates of the clusters
        in an X / Y table
    '''
    print("-----------------------------------------------")
    print("| Clusters para a região: '%s' "%(rName))
    print("-----------------------------------------------")
    print("| Idx  |       Long       |         Lat       |")
    print("-----------------------------------------------")
    cnt=1
    for c in clusters:
        print("| %3d  |  %2.10f  |   %2.10f  |"%(cnt, c[0], c[1]))
        cnt=cnt+1
    print("-----------------------------------------------")

def drawOnMap(mapCenterLat, mapCenterLong, coords, clusters, clindexes,                          n_clusters=3):
    ''' drawOnMap
        Draws all coordinates and their respective clusters
         on a map centered on (mapCenterLat, mapCenterLong).
         Each cluster (and all points that belong to it)
         receives a specific color chosen from its index
         Note: the number of clusters cannot be greater than length
         the color pattern list (PYCLUSTER_COLORS)
         Inputs:
             mapCenterLat, mapCenterLong: center of the map;
             coords: list of points (lat, long)
             clusters: list of clusters
             clindexes: clusters index
             n_clusters: number of clusters (default = 3)
         Outputs:
             mymap: drawn map
    '''
    mymap = fl.Map(location=[mapCenterLat, mapCenterLong], zoom_start=11)
    cnt = 0
    for longi, lat in coords:
        # select the color of the point from the cluster index
        color = PYCLUSTER_COLORS[clindexes[cnt]%n_clusters]
        # creates the bookmark
        fl.CircleMarker(location=[lat, longi],
                            radius=10,
                            fill_color=color,
                            color=color,
                            opacity=0.5
                        ).add_to(mymap)
        cnt = cnt+1

    # plot the centroids
    color = 'black'
    for centroid in clusters:
        longi = centroid[0]
        lat = centroid[1]
        fl.CircleMarker(location=[lat, longi],
                            radius=16,
                            fill_color=color,
                            color=color,
                            opacity=0.9
                        ).add_to(mymap)
   
    return mymap

def kmeansSearchAndDisplay(regionName, regionCSVFile, nclusters, 
                           centerLongitude, centerLatitude, outdir="tmp/"):

    # tests whether the file exists
    if os.path.isfile(regionCSVFile) is False:
        print("> CSV da região '%s' não encontrado!"%(regionName))

    # reads the csv file with the input coordinates
    coords = pd.read_csv(regionCSVFile)
    # separates the latitude and longitude vectors
    lat = coords['Latitude']
    longi = coords['Longitude']
    # stacks the coordinate array in a numpy array
    coords = np.column_stack((longi, lat))

    # calculates the mean of the coordinate set
    coords_mean = np.mean(coords, 0)
    # calculates the standard deviation of the coordinates
    coords_std = np.std(coords, 0) 
    # calculates the mean of the circumscribed ellipse
    ellipse_mean = patches.Ellipse([coords_mean[0], coords_mean[1]], coords_std[0]*2, coords_std[1]*2, color = 'purple', alpha = 0.1)

    # finds the nclusters for the coordinate set
    clusters, clindexes = findClusters(coords, nclusters)

    # prints the coordinates of the clusters in a table
    printClusterCoords(clusters, regionName)
    

    # plot in 2D graphics
    plotClusterCharts(coords, coords_mean, ellipse_mean, clusters, 
    clindexes, nclusters, regionName, outdir)

    # creates a map and draws the clusters and coordinates on it
    mymap = drawOnMap(centerLongitude, centerLatitude, coords, clusters, clindexes, nclusters)
    
    # saves the map to an HTML file
    if outdir is not None:
        mapFile = outdir+"/map_%s.html"%(regionName)
    else:
        mapFile = "map_%s.html"%(regionName)
    mymap.save(mapFile)


if __name__ == '__main__':
    start_time = time.time()
    print("> Search Keywords in SSA...")
    # Parameters
    coordinates = ['-12.980955, -38.482528', '-13.005153, -38.523342', '-13.007734, -38.492368', 
                '-12.932579, -38.479067', '-12.903903, -38.420924', '-12.920871, -38.353468',
                '-12.886063, -38.319027', '-12.862494, -38.291141',
                '-12.845420, -38.446882', '-12.827465, -38.472917',
                '-12.787177, -38.400631', '-12.983141, -38.466163'] # all in ssa

    # keywords = ['gas+station']
    keywords = ['gas+station', 'posto+gasolina']
    radius = '2000'
    api_key = '' # insert your Places API here
    searchPlacesAPI(coordinates, keywords, radius, 'tmp/outSSA.csv', api_key)
    end_time = time.time()
    print("> Search time: %.3f s"%(1000*(end_time - start_time)))

    print("> Now display on maps: ")
    # coordenadas do centro de Salvador
    MAP_CENTER = [-12.970748228016372, -38.4761939536224]
    nclusters=5
    kmeansSearchAndDisplay("SalvadorCity", 'tmp/outSSA.csv', nclusters, 
                            MAP_CENTER[0], MAP_CENTER[1], outdir="tmp/")
    print("> Finished. Please check tmp/ dir to see the results")
