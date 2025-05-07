""" volInt.c 

Volume Integration, translated into Python from: 

	  volInt.c                                            
	                                                      
	  This code computes volume integrals needed for      
	  determining mass properties of polyhedral bodies.   
	                                                      
	  For more information, see the accompanying README   
	  file, and the paper                                 
	                                                      
	  Brian Mirtich, "Fast and Accurate Computation of    
	  Polyhedral Mass Properties," journal of graphics    
	  tools, volume 1, number 1, 1996.                    
	                                                      
	  This source code is public domain, and may be used  
	  in any way, shape or form, free of charge.          
	                                                      
	  Copyright 1995 by Brian Mirtich                     
	                                                      
	  mirtich@cs.berkeley.edu                             
	  http://www.cs.berkeley.edu/~mirtich                 
"""

import numpy as np


def read_polyhedron(File):
    #######################
    # Load ,txt data file 
    data = np.loadtxt(File, delimiter=' ', dtype=str) 
    ###############################################################
    # Set Vertex/Faces denoted as v or f in .obj format to array 
    vertex_faces = data[:,0]                                      
    # Get Length of the Vertex/Faces array for range counting     
    V_F_Range = vertex_faces.size                                 
    # Define variable for number of vertices & faces              
    numb_vert = 0                                                 
    numb_face = 0                                                 
    # Scan Data for v & f and count the numbers of each.          
    #  Used for sorting x, y, & z as vertices                     
    for i in range(0,V_F_Range):                                  
        if vertex_faces[i] == 'v':                                
            numb_vert += 1                                        
        else:                                                     
            numb_face += 1                                        
    #########################
    # Assigning Vertex Data 
    ###################################
    # Assign 2nd row of .txt as x input    
    x_input = data[range(0,numb_vert),1]   
    # Assign 3rd row of .txt as y input    
    y_input = data[range(0,numb_vert),2]   
    # Assign 4th row of .txt as z input    
    z_input = data[range(0,numb_vert),3]   
    # Convert Vertices data to float type  
    x = x_input.astype(float)             
    y = y_input.astype(float)            
    z = z_input.astype(float)           
    #########################
    # Assigning Face Data 
    #############################################
    # Range count for face data                 
    row_tot = numb_face + numb_vert             
    # Assign 2nd row of .txt as x input         
    fx_input = data[range(numb_vert,row_tot),1] 
    # Assign 3rd row of .txt as y input         
    fy_input = data[range(numb_vert,row_tot),2] 
    # Assign 4th row of .txt as z input         
    fz_input = data[range(numb_vert,row_tot),3] 
    # Convert Vertices data to float type       
    fx = fx_input.astype(str)                   
    fy = fy_input.astype(str)                  
    fz = fz_input.astype(str)                 
    # Define: number of vertices on ith face 
    #  - used for .C program                 
    ith_face = []                            
    for j in range(numb_vert,row_tot):       
        ith_face.append(3)                   
    ith_array = np.array(ith_face)           
    #########################################
    # Number of Vertex set to array         
    numb_vert_array = []                    
    numb_vert_array.append(numb_vert)       
    # Number of Faces set to array          
    numb_face_array = []                    
    numb_face_array.append(numb_face)       
    # Stacking Columns of Vertex Data      
    Vert_Data_Out_0 = np.column_stack((x, y))               
    Vert_Data_Out   = np.column_stack((Vert_Data_Out_0, z)) 
    # Stacking Columns of Face Data                         
    Face_Data_Out_0 = np.column_stack((ith_array, fx))      
    Face_Data_Out_Y = np.column_stack((Face_Data_Out_0,fy)) 
    Face_Data_Out   = np.column_stack((Face_Data_Out_Y,fz)) 
    ########################################################
    return Vert_Data_Out, Face_Data_Out, numb_vert, numb_face


def compFaceNorm(V,F):
    nF = len(F)
    norm = np.zeros((3,nF))
    
    for it in range(nF): 
        dx1 = V[F[it,1],0] - V[F[it,0],0]
        dy1 = V[F[it,1],1] - V[F[it,0],1]
        dz1 = V[F[it,1],2] - V[F[it,0],2]
        dx2 = V[F[it,2],0] - V[F[it,1],0]
        dy2 = V[F[it,2],1] - V[F[it,1],1]
        dz2 = V[F[it,2],2] - V[F[it,1],2]
        #####
        # Face normal vector
        nx = dy1*dz2 - dy2*dz1
        ny = dz1*dx2 - dz2*dx1
        nz = dx1*dy2 - dx2*dy1 
        ####
        # magnitude
        n_mag = np.sqrt(nx**2 + ny**2 + nz**2)
        #######
        # norms
        # 
        # x 
        norm[it,0] = nx/n_mag
        # y
        norm[it,1] = ny/n_mag
        # z 
        norm[it,2] = nz/n_mag 
        #####################
        # 
        w = -norm[it,0]*V[F[it,0],0] \
            -norm[it,0]*V[F[it,0],1] \
            -norm[it,0]*V[F[it,0],2] 
        return norm, w




def comp_project_int(V,F):
    nV = len(V)
	###
	# Constants
    P1 = Pa = Pb = Paa = Pab = Pbb = Paaa = Paab = Pabb = Pbbb = 0.0

    for it2 in range(nV):
        a0 = 0.0# Filler
        b0 = 0.0# This looks like verts index, 
        a1 = 0.0# but varying in volInt.c from 
        b1 = 0.0# X,Y,Z array notation 
        ###
        da   = a1 - a0
        db   = b1 - b0
        #
        a0_2 = a0*a0
        a0_3 = a0_2*a0
        a0_4 = a0_3*a0
        #
        b0_2 = b0*b0
        b0_3 = b0_2*b0
        b0_4 = b0_3*b0
        # 
        a1_2 = a1*a1
        a1_3 = a1_2*a1
        #
        b1_2 = b1*b1
        b1_3 = b1_2*b1
        ##############
        # 
        C1   = a1 + a0
        Ca   = a1*C1 + a0_2 
        Caa  = a1*Ca + a0_3
        Caaa = a1*Caa + a0_4
        # 
        Cb   = b1*(b1+b0) + b0_2
        Cbb  = b1*Cb + b0_3 
        Cbbb = b1*Cbb + b0_4
        # 
        Cab  = 3*a1_2 + 2*a1*a0 + a0_2
        Kab  = a1_2 + 2*a1*a0 + 3*a0_2
        #
        Caab = a0*Cab + 4*a1_3
        Kaab = a1*Kab + 4*a0_3
        #
        Cabb = 4*b1_3 + 3*b1_2*b0 + 2*b1*b0_2 + b0_3
        Kabb = b1_3 + 2*b1_2*b0 + 3*b1*b0_2 + 4*b0_3
        ###################
        # Summ the coeff 
        P1   += db*C1
        Pa   += db*Ca
        Paa  += db*Caa
        Paaa += db*Caaa
        Pb   += da*Cb
        Pbb  += da*Cbb
        Pbbb += da*Cbbb
        Pab  += db*(b1*Cab + b0*Kab)
        Paab += db*(b1*Caab + b0*Kaab)
        Pabb += da*(a1*Cabb + a0*Kabb)
    ###########
    P1   = P1/2.0
    Pa   = Pa/6.0
    Paa  = Paa/12.0
    Paaa = Paaa/20.0
    Pb   = Pb/-6.0
    Pbb  = Pbb/-12.0
    Pbbb = Pbbb/-20.0
    #################
    return P1, Pa, Pb, Paa, Pab, Pbb, Paaa, Paab, Pabb, Pbbb




if __name__ == "__main__":
    ####
    # Replace with user file selection
    File = 'Apophis.obj'
    #####
    # Read in file 
    V, F, nV, nF = read_polyhedron(File)
    #####
    # Compute Face Normals
    norm , w = compFaceNorm(V,F)
    ####
    # 
    