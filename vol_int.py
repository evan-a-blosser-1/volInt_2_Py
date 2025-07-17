import numpy as np
from dataclasses import dataclass
import sys

# Constants
X, Y, Z = 0, 1, 2
MAX_VERTS = 8192
MAX_FACES = 16380
MAX_POLYGON_SZ = 10

@dataclass
class Face:
    numVerts: int = 0
    norm: np.ndarray = np.zeros(3)
    w: float = 0.0
    verts: list = None
    
    def __post_init__(self):
        if self.verts is None:
            self.verts = []

@dataclass
class Polyhedron:
    numVerts: int = 0
    numFaces: int = 0
    verts: np.ndarray = None
    faces: list = None
    
    def __post_init__(self):
        if self.verts is None:
            self.verts = np.zeros((MAX_VERTS, 3))
        if self.faces is None:
            self.faces = [Face() for _ in range(MAX_FACES)]


def read_polyhedron(File):
    """Read polyhedron data and return structured object"""
    # Create polyhedron object
    p = Polyhedron()
    
    # Load data file
    data = np.loadtxt(File, delimiter=' ', dtype=str)
    vertex_faces = data[:,0]
    V_F_Range = vertex_faces.size
    
    # Count vertices and faces
    p.numVerts = sum(1 for x in vertex_faces if x == 'v')
    p.numFaces = V_F_Range - p.numVerts
    
    # Read vertices
    vert_data = data[:p.numVerts, 1:4].astype(float)
    p.verts = vert_data
    
    # Read faces and compute normals
    face_data = data[p.numVerts:, 1:4].astype(int)
    for i in range(p.numFaces):
        face = p.faces[i]
        face.numVerts = 3  # Assuming triangular faces
        face.verts = face_data[i].tolist()
        
        # Get vertices for normal computation
        v0 = p.verts[face.verts[0]-1]
        v1 = p.verts[face.verts[1]-1]
        v2 = p.verts[face.verts[2]-1]
        
        # Compute edges
        d1 = v1 - v0
        d2 = v2 - v1
        
        # Compute normal
        nx = d1[1]*d2[2] - d2[1]*d1[2]
        ny = d1[2]*d2[0] - d2[2]*d1[0]
        nz = d1[0]*d2[1] - d2[0]*d1[1]
        
        # Normalize
        length = np.sqrt(nx*nx + ny*ny + nz*nz)
        face.norm = np.array([nx, ny, nz]) / length
        
        # Compute w offset
        face.w = -np.dot(face.norm, v0)
    print(f"Read {p.numVerts} vertices and {p.numFaces} faces from {File}")
    return p
                

def comp_projection_integrals(face, verts):
    """Compute various integrations over projection of face"""
    global P1, Pa, Pb, Paa, Pab, Pbb, Paaa, Paab, Pabb, Pbbb
    
    P1 = Pa = Pb = Paa = Pab = Pbb = Paaa = Paab = Pabb = Pbbb = 0.0
    
    for i in range(face.numVerts):
        j = (i + 1) % face.numVerts
        
        a0 = verts[face.verts[i]-1][A]
        b0 = verts[face.verts[i]-1][B]
        a1 = verts[face.verts[j]-1][A]
        b1 = verts[face.verts[j]-1][B]
        
        da = a1 - a0
        db = b1 - b0
        
        # Powers
        a0_2 = a0 * a0
        a0_3 = a0_2 * a0
        a0_4 = a0_3 * a0
        b0_2 = b0 * b0 
        b0_3 = b0_2 * b0
        b0_4 = b0_3 * b0
        a1_2 = a1 * a1
        a1_3 = a1_2 * a1
        b1_2 = b1 * b1
        b1_3 = b1_2 * b1
        
        # Coefficients
        C1 = a1 + a0
        Ca = a1*C1 + a0_2
        Caa = a1*Ca + a0_3
        Caaa = a1*Caa + a0_4
        Cb = b1*(b1 + b0) + b0_2
        Cbb = b1*Cb + b0_3
        Cbbb = b1*Cbb + b0_4
        Cab = 3*a1_2 + 2*a1*a0 + a0_2
        Kab = a1_2 + 2*a1*a0 + 3*a0_2
        Caab = a0*Cab + 4*a1_3
        Kaab = a1*Kab + 4*a0_3
        Cabb = 4*b1_3 + 3*b1_2*b0 + 2*b1*b0_2 + b0_3
        Kabb = b1_3 + 2*b1_2*b0 + 3*b1*b0_2 + 4*b0_3
        
        # Update integrals
        P1 += db*C1
        Pa += db*Ca 
        Paa += db*Caa
        Paaa += db*Caaa
        Pb += da*Cb
        Pbb += da*Cbb
        Pbbb += da*Cbbb
        Pab += db*(b1*Cab + b0*Kab)
        Paab += db*(b1*Caab + b0*Kaab)
        Pabb += da*(a1*Cabb + a0*Kabb)
    
    # Final scaling
    P1 /= 2.0
    Pa /= 6.0 
    Paa /= 12.0
    Paaa /= 20.0
    Pb /= -6.0
    Pbb /= -12.0
    Pbbb /= -20.0
    Pab /= 24.0
    Paab /= 60.0
    Pabb /= -60.0
    
    
    
def comp_face_integrals(face, norm, w):
    """Compute face integrals"""
    global P1, Pa, Pb, Paa, Pab, Pbb, Paaa, Paab, Pabb, Pbbb
    global Fa, Fb, Fc, Faa, Fbb, Fcc, Faaa, Fbbb, Fccc, Faab, Fbbc, Fcca
    
    comp_projection_integrals(face, face.verts)
    
    k1 = 1 / norm[C] 
    k2 = k1 * k1
    k3 = k2 * k1 
    k4 = k3 * k1

    Fa = k1 * Pa
    Fb = k1 * Pb
    Fc = -k2 * (norm[A]*Pa + norm[B]*Pb + w*P1)

    Faa = k1 * Paa
    Fbb = k1 * Pbb
    Fcc = k3 * (np.power(norm[A], 2)*Paa + 2*norm[A]*norm[B]*Pab + 
                np.power(norm[B], 2)*Pbb + w*(2*(norm[A]*Pa + norm[B]*Pb) + w*P1))

    Faaa = k1 * Paaa
    Fbbb = k1 * Pbbb
    Fccc = -k4 * (np.power(norm[A], 3)*Paaa + 
                  3*np.power(norm[A], 2)*norm[B]*Paab +
                  3*norm[A]*np.power(norm[B], 2)*Pabb + 
                  np.power(norm[B], 3)*Pbbb +
                  3*w*(np.power(norm[A], 2)*Paa + 
                       2*norm[A]*norm[B]*Pab + 
                       np.power(norm[B], 2)*Pbb) +
                  w*w*(3*(norm[A]*Pa + norm[B]*Pb) + w*P1))

    Faab = k1 * Paab
    Fbbc = -k2 * (norm[A]*Pabb + norm[B]*Pbbb + w*Pbb)
    Fcca = k3 * (np.power(norm[A], 2)*Paaa + 
                 2*norm[A]*norm[B]*Paab + 
                 np.power(norm[B], 2)*Pabb +
                 w*(2*(norm[A]*Paa + norm[B]*Pab) + w*Pa))

def comp_volume_integrals(polyhedron):
    """Compute volume integrals for polyhedron"""
    global A, B, C, T0, T1, T2, TP
    
    # Initialize arrays
    T0 = 0
    T1 = np.zeros(3)
    T2 = np.zeros(3)
    TP = np.zeros(3)

    for face in polyhedron.faces:
        # Get absolute values of normal components
        nx = abs(face.norm[X])
        ny = abs(face.norm[Y])
        nz = abs(face.norm[Z])
        
        # Choose axis with largest normal component
        if nx > ny and nx > nz:
            C = X
        else:
            C = Y if ny > nz else Z
        A = (C + 1) % 3
        B = (A + 1) % 3

        comp_face_integrals(face, face.norm, face.w)

        T0 += face.norm[X] * (Fa if A == X else (Fb if B == X else Fc))
        
        T1[A] += face.norm[A] * Faa
        T1[B] += face.norm[B] * Fbb
        T1[C] += face.norm[C] * Fcc
        T2[A] += face.norm[A] * Faaa
        T2[B] += face.norm[B] * Fbbb
        T2[C] += face.norm[C] * Fccc
        TP[A] += face.norm[A] * Faab
        TP[B] += face.norm[B] * Fbbc
        TP[C] += face.norm[C] * Fcca

    # Scale the results
    T1 = T1 / 2.0
    T2 = T2 / 3.0
    TP = TP / 2.0
    
    return T0, T1, T2, TP

def main():
    """Main function"""
    file = input("Enter the polyhedron file name: ") 
    # Read polyhedron from file
    polyhedron = read_polyhedron(file)
    
    # Compute volume integrals
    T0, T1, T2, TP = comp_volume_integrals(polyhedron)

    # Print results
    print(f"\nT1 = {T0:+20.6f}\n")
    print(f"Tx = {T1[X]:+20.6f}")
    print(f"Ty = {T1[Y]:+20.6f}")
    print(f"Tz = {T1[Z]:+20.6f}\n")
    
    print(f"Txx = {T2[X]:+20.6f}")
    print(f"Tyy = {T2[Y]:+20.6f}")
    print(f"Tzz = {T2[Z]:+20.6f}\n")
    
    print(f"Txy = {TP[X]:+20.6f}")
    print(f"Tyz = {TP[Y]:+20.6f}")
    print(f"Tzx = {TP[Z]:+20.6f}\n")

    # Compute mass properties
    density = 1.0  # assume unit density
    mass = density * T0

    # Compute center of mass
    r = T1 / T0

    # Compute inertia tensor
    J = np.zeros((3, 3))
    J[X][X] = density * (T2[Y] + T2[Z])
    J[Y][Y] = density * (T2[Z] + T2[X])
    J[Z][Z] = density * (T2[X] + T2[Y])
    J[X][Y] = J[Y][X] = -density * TP[X]
    J[Y][Z] = J[Z][Y] = -density * TP[Y]
    J[Z][X] = J[X][Z] = -density * TP[Z]

    # Translate inertia tensor to center of mass
    J[X][X] -= mass * (r[Y]*r[Y] + r[Z]*r[Z])
    J[Y][Y] -= mass * (r[Z]*r[Z] + r[X]*r[X])
    J[Z][Z] -= mass * (r[X]*r[X] + r[Y]*r[Y])
    # Correct way to update symmetric inertia tensor elements
    J[X][Y] += mass * r[X] * r[Y]
    J[Y][X] = J[X][Y]  
    J[Y][Z] += mass * r[Y] * r[Z]
    J[Z][Y] = J[Y][Z] 
    J[Z][X] += mass * r[Z] * r[X]
    J[X][Z] = J[Z][X]  

    print(f"center of mass: ({r[X]:+12.6f},{r[Y]:+12.6f},{r[Z]:+12.6f})\n")
    print("inertia tensor with origin at c.o.m.:")
    print(f"{J[X][X]:+15.6f} {J[X][Y]:+15.6f} {J[X][Z]:+15.6f}")
    print(f"{J[Y][X]:+15.6f} {J[Y][Y]:+15.6f} {J[Y][Z]:+15.6f}")
    print(f"{J[Z][X]:+15.6f} {J[Z][Y]:+15.6f} {J[Z][Z]:+15.6f}\n")

if __name__ == "__main__":
    main()