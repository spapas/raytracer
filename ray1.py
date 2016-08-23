import numpy as np
from PIL import Image

class Ray(object, ):
    A = None
    B = None

    def __init__(self, a=None, b= None):
        self.A = a
        self.B = b

    @property
    def origin(self):
        return self.A

    @property
    def direction(self):
        return self.B

    def point_at_parameter(self, t):
        return self.A + t*self.B

class HitRecord(object, ):

    def __init__(self, t, p, normal):
        self.t = t
        self.p = p
        self.normal = normal

class Hitable(object, ):
    def hit(self, r, t_min, t_max, ):
        pass

class Sphere(Hitable, ):

    def __init__(self, c, r):
        self.center = c
        self.radius = r

    def hit(self, r, t_min, t_max, ):
        oc = r.origin - self.center
        a = r.direction.dot(r.direction)
        b = oc.dot(r.direction)
        c = oc.dot(oc) - self.radius * self.radius
        discriminant = b*b - a*c

        if discriminant > 0:
            temp = ( -b - np.sqrt(b*b-a*c))/a
            if temp < t_max and temp > t_min:
                p = r.point_at_parameter(temp)
                return True, {
                    't': temp,
                    'p': p,
                    'normal': (p-self.center)/self.radius
                }
            temp = ( -b + np.sqrt(b*b-a*c))/a
            if temp < t_max and temp > t_min:
                p = r.point_at_parameter(temp)
                return True, {
                    't': temp,
                    'p': p,
                    'normal': (p-self.center) / self.radius
                }
        return False, {}

class HitableList(Hitable, ):
    hitables = []
    
    def __init__(self, hitables):
        self.hitables = hitables

    def hit(self, r, t_min, t_max, ):
        hit_anything = False
        closest_so_far = t_max
        rec = {}
        for h in self.hitables:
            res, tmp_rec = h.hit(r, t_min, closest_so_far)
            if res:
                hit_anything = True
                closest_so_far = tmp_rec['t']
                rec = tmp_rec
                
        return hit_anything, rec


def unit_vector(v):
    mag = np.sqrt(v.dot(v))
    return np.array([v[0]/mag, v[1]/mag, v[2]/mag])

def color(r, world):
    
    res, rec = world.hit(r, 0.0, 999999)
    if(res):
        return 0.5 * np.array([
            rec['normal'][0] + 1,
            rec['normal'][1] + 1,
            rec['normal'][2] + 1,
        ])
    else:
        unit_direction = unit_vector(r.direction)
        t = 0.5*(unit_direction[1] + 1.0)
        return (1.0-t)*np.array([1.0, 1.0, 1.0])+ t * np.array([0.5, 0.7, 1.0])

def hit_sphere11(center, radius, r):
    oc = r.origin - center
    a = r.direction.dot(r.direction)
    b = 2.0 * oc.dot(r.direction)
    c = oc.dot(oc) - radius * radius
    discriminant = b * b - 4 * a * c
    if discriminant < 0 :
        return -1.0
    else:
        return (-b - np.sqrt(discriminant)) / (2.0 * a)



def create_image(stream, nx=200, ny=100):
    stream.write("P6\n")
    stream.write("{0} {1}\n".format(nx,ny))
    stream.write("255\n")

    lower_left_corner = np.array([-2.0, -1.0, -1.0])
    horizontal = np.array([4.0, 0.0, 0.0])
    vertical = np.array([0.0, 2.0, 0.0])
    origin = np.array([0.0, 0.0, 0.0])

    world = HitableList([
        Sphere(np.array([0, 0, -1]), 0.5),
        Sphere(np.array([1, 1,-1.5]), 0.5),
        Sphere(np.array([ 0.5, -1.5,-3]), 0.5),
        #Sphere(np.array([ -0.5, -0.5, -1]), 1.0),
        #Sphere(np.array([0,-100,5]), 100.5),
        
    ])
    
    for y in range(ny):
        for x in range(nx):

            u = x*1.0 / nx
            v = y*1.0 / ny
            r = Ray(origin, lower_left_corner+ u*horizontal + v*vertical)
            col = color(r, world)

            ir = int(255.99*col[0])
            ig = int(255.99*col[1])
            ib = int(255.99*col[2])

            stream.write(bytearray([ir, ig, ib]))
        #stream.write('\n')
            

if __name__ == '__main__':
    f = open('f1.ppm', "wb")
    create_image(f, 600, 300)
    f.close()
    
    im = Image.open("f1.ppm")
    print im.format, im.size, im.mode
    im.save("f1.png")