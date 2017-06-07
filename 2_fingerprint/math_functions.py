# BINANA is released under the GNU General Public License (see http://www.gnu.org/licenses/gpl.html).
# If you have any questions, comments, or suggestions, please don't hesitate to contact me, Jacob
# Durrant, at jdurrant [at] ucsd [dot] edu. If you use BINANA in your work, please cite 
# Durrant, J. D. and J. A. McCammon (2011). "BINANA: A novel algorithm for ligand-binding
# characterization." J Mol Graph Model 29(6): 888-893.

# The below is extracted from BINANA with stylistic changes made by Joe Paggi.
# I make no claims about the correctness of the below code.

import math
from point_atom import Point

def sigmoid(x, center, compression):
    e = math.exp(compression * (x - center))
    return e / (1.0 + e) 

def vector_subtraction(vector1, vector2):
    return Point(vector1.x - vector2.x, vector1.y - vector2.y, vector1.z - vector2.z)
    
def CrossProduct(Pt1, Pt2): # never tested
    Response = Point(0,0,0)
    
    Response.x = Pt1.y * Pt2.z - Pt1.z * Pt2.y
    Response.y = Pt1.z * Pt2.x - Pt1.x * Pt2.z
    Response.z = Pt1.x * Pt2.y - Pt1.y * Pt2.x
    
    return Response;
    
def vector_scalar_multiply(vector, scalar):
    return Point(vector.x * scalar, vector.y * scalar, vector.z * scalar)
    
def dot_product(point1, point2):
    return point1.x * point2.x + point1.y * point2.y + point1.z * point2.z
    
def dihedral(point1, point2, point3, point4): # never tested
    b1 = vector_subtraction(point2, point1)
    b2 = vector_subtraction(point3, point2)
    b3 = vector_subtraction(point4, point3)
    
    b2Xb3 = CrossProduct(b2,b3)
    b1Xb2 = CrossProduct(b1,b2)
    
    b1XMagb2 = vector_scalar_multiply(b1,b2.magnitude())
    radians = math.atan2(dot_product(b1XMagb2,b2Xb3), dot_product(b1Xb2,b2Xb3))
    return radians
    
def angle_between_three_points(point1, point2, point3): # As in three connected atoms
    vector1 = vector_subtraction(point1, point2)
    vector2 = vector_subtraction(point3, point2)
    return angle_between_points(vector1, vector2)

def return_normalized_vector(vector):
    dist = distance(Point(0,0,0), vector)
    return Point(vector.x/dist, vector.y/dist, vector.z/dist)

    
def angle_between_points(point1, point2):
    new_point1 = return_normalized_vector(point1)
    new_point2 = return_normalized_vector(point2)
    dot_prod = dot_product(new_point1, new_point2)
    if dot_prod > 1.0: dot_prod = 1.0 # to prevent errors that can rarely occur
    if dot_prod < -1.0: dot_prod = -1.0
    return math.acos(dot_prod)
    
def distance(point1, point2):
    deltax = point1.x - point2.x
    deltay = point1.y - point2.y
    deltaz = point1.z - point2.z
    
    return math.sqrt(math.pow(deltax,2) + math.pow(deltay,2) + math.pow(deltaz,2))
        
def project_point_onto_plane(apoint, plane_coefficients): # essentially finds the point on the plane that is closest to the specified point
    # the plane_coefficients are [a,b,c,d], where the plane is ax + by + cz = d
    
    # First, define a plane using cooeficients a, b, c, d such that ax + by + cz = d
    a = plane_coefficients[0]
    b = plane_coefficients[1]
    c = plane_coefficients[2]
    d = plane_coefficients[3]
        
    # Now, define a point in space (s,u,v)
    s = apoint.x
    u = apoint.y
    v = apoint.z
    
    t = (d - a*s - b*u - c*v) / (a*a + b*b + c*c)
        
    # here's the point closest on the plane
    x = s + a*t
    y = u + b*t
    z = v + c*t
        
    return Point(x,y,z)
