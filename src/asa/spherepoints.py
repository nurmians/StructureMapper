#!/usr/bin/python
# coding=UTF-8
# -*- coding: UTF-8 -*-

# This file is part of the StructureMapper algorithm.
# Please cite the authors if you find this software useful
#
# https://academic.oup.com/bioinformatics/advance-article/doi/10.1093/bioinformatics/bty086/4857361
 
# MIT License
#
# Copyright 2018 Anssi Nurminen and Vesa P. Hytönen
# Faculty of Medicine and Life Sciences and BioMediTech, University of Tampere, Arvo Ylpön katu 34, 33520 Tampere, Finland 
# Fimlab Laboratories, Biokatu 4, 33520 Tampere, Finland
# Tampere, Finland
# 16-02-2018

# Permission is hereby granted, free of charge, to any person obtaining
# a copy of this software and associated documentation files (the "Software"),
# to deal in the Software without restriction, including without limitation the
# rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
# sell copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all 
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
# THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.


import math
    
def NumberOfIcosahedricPoints( multiplier = 3):
    return ( 10 * pow(4, multiplier-1) + 2 )
    
def IcosahedronPointsOnSphere( multiplier = 3, verbose=False):

    if multiplier < 1 or multiplier > 20:
        multiplier = 3;
        print "Multiplier for generating points should be in range 1-20, using default multiplier %i" % multiplier


    #set corners of the icosahedron
    tao = 1.61803399; #golden ratio

    #Points of an icosahedron (12)
    #        0          1           2          3          4            5           6           7          8           9            10          11
    icos = [ [0,1,tao], [0,-1,tao], [tao,0,1], [1,tao,0], [-1,tao,0], [-tao,0,1], [1,-tao,0], [-1,-tao,0], [tao,0,-1], [-tao,0,-1], [0,1,-tao], [0,-1,-tao] ];
                    
    #Normalize vectors to a length of 1
    for i in range( len(icos)):
        icos[ i] = NormalizeVector( icos[ i])

    #all the lines in the icosahedron (30)
    #the points of the lines are kept as indexes into icos to eliminate the need to compare floating point coordinate equalities
    #          0       1       2       3       4
    lines =  [( 0, 1),( 2, 0),( 3, 0),( 0, 4),( 5, 0),
    #          5       6       7       8       9
              ( 6,11),(11, 7),(11, 8),( 9,11),(10,11),
    #          10      11      12      13      14
              ( 1, 2),( 2, 3),( 3, 4),( 4, 5),( 5, 1),
    #          15      16      17      18      19
              ( 6, 8),( 8,10),(10, 9),( 9, 7),( 7, 6),            
    #          20      21      22      23      24            
              ( 7, 1),( 1, 6),( 6, 2),( 2, 8),( 8, 3),
    #          25      26      27      28      29            
              ( 3,10),(10, 4),( 4, 9),( 9, 5),( 5, 7)]
            
    #lines that together form a triangular side of the icosahedron  (20)
    triangles = [ [0,10,1],  [1,11,2],  [2,12,3],  [3,13,4],  [4,14,0],
                  [5,15,7],  [7,16,9],  [9,17,8],  [8,18,6],  [6,19,5],
                  [21,22,10],[23,24,11],[25,26,12],[27,28,13],[29,20,14],
                  [22,23,15],[24,25,16],[26,27,17],[28,29,18],[20,21,19]]
            
    
    if verbose:
        print "Multiplier %i: Generating %i icosaqhedric sphere points..." % ( multiplier, NumberOfIcosahedricPoints( multiplier))
    
    while multiplier > 1:
        
        newlines = []
        newpoints = []
          
        num_of_old_points = len( icos)
        
        for l in range( len(lines)):
            #Add midway point of each line and set its distance from the origo to 1.0 (normalize)
            icos.append( NormalizeVector( MidwayPoint( icos[ lines[ l][ 0]], icos[ lines[ l][ 1]])))
            #line is split into two lines
            newlines.append( (lines[ l][ 0], num_of_old_points + l))
            newlines.append( (lines[ l][ 1], num_of_old_points + l))           
        
        #are we done yet?
        multiplier -= 1
        if multiplier <= 1: break            
            
        #Every triangle is split into 4 new triangles         
        #
        #    /\        /\
        #   /__\      /1 \
        #  /\  /\    / 4  \
        # /__\/__\  /_2__3_\
        #
        
        newtriangles = []      
        triangle_index = 0;          
        
        for t in range( len(triangles)):
                                                                   
            
            for i in range( 3):        
                line_A_index = triangles[ t][ i]
                line_B_index = triangles[ t][ (i+1) % 3]
          
                shared_A_index = 0 if lines[ line_A_index][ 0] == lines[ line_B_index][ 0] or lines[ line_A_index][ 0] == lines[ line_B_index][ 1] else 1
                shared_B_index = 0 if lines[ line_B_index][ 0] == lines[ line_A_index][ 0] or lines[ line_B_index][ 0] == lines[ line_A_index][ 1] else 1
                
                newlines.append( (line_A_index + num_of_old_points, line_B_index + num_of_old_points))
                #Insert corner triangle (x3)
                newtriangles.append( [ line_A_index*2+shared_A_index, line_B_index*2+shared_B_index, len( newlines)-1 ])
                
            #Insert middle triangle
            size = len( newlines)
            newtriangles.append( range( size-3, size))     
        

        triangles = newtriangles
        lines = newlines        
        
    if verbose:
        print "Generated %i icosaqhedric sphere points\n" % len( icos)
        
    #return generated points
    return icos


      
def MidwayPoint( vector1, vector2):
    return [(vector2[ 0] + vector1[ 0])/2, (vector2[ 1] + vector1[ 1])/2, (vector2[ 2] + vector1[ 2])/2]         
  
def VectorLength( vector):
    return math.sqrt( vector[ 0]*vector[ 0] + vector[ 1]*vector[ 1] + vector[ 2]*vector[ 2])  
  
def NormalizeVector( vector):
    vl = VectorLength( vector)
    return [vector[ 0]/vl, vector[ 1]/vl, vector[ 2]/vl]
  

    
#Octave is a open source matlab equivivalent
#This methods prints out generates points in a Octave plottable format
def PrintOctavePlot( pts):

    
    matrix = "M = ["
    for p in range( len(pts)):
        matrix += ("%f, %f, %f;\n" % (pts[ p][ 0], pts[ p][ 1], pts[ p][2]))
    matrix += "];"    
    
    print matrix
    
    print "clf;"
    print "scatter3(M(:,1),M(:,2),M(:,3));"
    #print "plotmatrix( M, \"x\", \"color\", \"blue\");"
    
    #print "clf; hold on;"
    
    #for p in range( len(pts)):
    #    if p < 12:
    #        print "plot3( %f, %f, %f , \"o\", \"color\", \"blue\");" % (pts[ p][ 0], pts[ p][ 1], pts[ p][2])
    #    print "plot3( %f, %f, %f , \"x\", \"color\", \"blue\");" % (pts[ p][ 0], pts[ p][ 1], pts[ p][2])
        #print "text(  %f, %f, %f , \"%i\", \"color\", \"black\");" % (pts[ p][ 0], pts[ p][ 1], pts[ p][2], p)  


#create sphere points
#Vogel method
def GoldenSpiralPointsOnSphere( numberOfPoints=920 ): 

    N = float( numberOfPoints) 
    pts = []   
    inc = math.pi * (3 - math.sqrt( 5)) 
    off = 2 / N 
    
    for k in range( 0, numberOfPoints):
         y = k * off - 1 + (off / 2)
         r = math.sqrt( 1 - y*y)
         phi = k * inc 
         pts.append( [math.cos(phi)*r, y, math.sin(phi)*r])   
         
    return pts        


def main():

    multiplier = 4;

    pts = IcosahedronPointsOnSphere( multiplier)
    #pts = GoldenSpiralPointsOnSphere( NumberOfIcosahedricPoints( multiplier))

    PrintOctavePlot( pts)
     
    #A=[1,tao,0; -1,tao,0; 1,-tao,0; -1,-tao,0; 0,1,tao; 0,-1,tao; 0,1,-tao; 0,-1,-tao; tao,0,1; -tao,0,1; tao,0,-1; -tao,0,-1];


if __name__ == "__main__":
  main()