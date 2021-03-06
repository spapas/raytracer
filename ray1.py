import sys
import numpy as np
import random
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
    def __init__(self, t, p, normal, material):
        self.t = t
        self.p = p
        self.normal = normal
        self.material = material

class Camera(object, ):
    lower_left_corner = np.array([-2.0, -1.0, -1.0])
    horizontal = np.array([4.0, 0.0, 0.0])
    vertical = np.array([0.0, 2.0, 0.0])
    origin = np.array([0.0, 0.0, 0.0])

    def get_ray(self, u, v):
        return Ray(self.origin, self.lower_left_corner+ u*self.horizontal + v*self.vertical)

class Hitable(object, ):
    def hit(self, r, t_min, t_max, ):
        pass

class Material(object, ):
    def scatter(self, r_in, hit_record, ):
        pass # return boolean, attenuation, scattered

class Lambertian(Material, ):
    def __init__(self, a):
        self.albedo = a

    def scatter(self, r_in, hit_record, ):
        target = hit_record.p + hit_record.normal + random_in_unit_sphere()
        scattered = Ray(hit_record.p, target-hit_record.p)
        return True, self.albedo, scattered

class Metal(Material, ):
    def __init__(self, a, f):
        self.albedo = a
        if f<1:
            self.fuzz = f
        else:
            self.fuzz = 1

    def scatter(self, r_in, hit_record, ):
        reflected = reflect(unit_vector(r_in.direction), hit_record.normal)
        scattered = Ray(hit_record.p, reflected + self.fuzz*random_in_unit_sphere() )
        return scattered.direction.dot(hit_record.normal)>0, self.albedo, scattered

class Sphere(Hitable, ):

    def __init__(self, c, r, m):
        self.center = c
        self.radius = r
        self.material = m

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
                return True, HitRecord(temp, p, (p-self.center)/self.radius, self.material)

            temp = ( -b + np.sqrt(b*b-a*c))/a
            if temp < t_max and temp > t_min:
                p = r.point_at_parameter(temp)
                return True, HitRecord(temp, p, (p-self.center)/self.radius, self.material)
        return False, None

class HitableList(Hitable, ):
    hitables = []

    def __init__(self, hitables):
        self.hitables = hitables

    def hit(self, r, t_min, t_max, ):
        hit_anything = False
        closest_so_far = t_max
        rec = None
        for h in self.hitables:
            res, tmp_rec = h.hit(r, t_min, closest_so_far)
            if res:
                hit_anything = True
                closest_so_far = tmp_rec.t
                rec = tmp_rec

        return hit_anything, rec

def reflect(v, n):
    return v - 2*v.dot(n)*n

def random_in_unit_sphere():
    while True:
        p = 2.0* np.array(
            [random.random(), random.random(), random.random(), ]
        ) - np.ones(3)
        if p.dot(p) >=1.0:
            break

    return p

def unit_vector(v):
    mag = np.sqrt(v.dot(v))
    return np.array([v[0]/mag, v[1]/mag, v[2]/mag])

def color(r, world, depth):

    res, rec = world.hit(r, 0.001, 999999999)
    if(res):
        res, attenuation, scattered = rec.material.scatter(r, rec, )
        #print scattered, type(scattered)
        if depth<50 and res:
            return attenuation*color(scattered, world, depth+1)

        else:
            return np.array([0.0, 0.0, 0.0])

    else:
        unit_direction = unit_vector(r.direction)
        t = 0.5*(unit_direction[1] + 1.0)
        return (1.0-t)*np.array([1,1,1]) + t * np.array([0.5, 0.7, 1.0])


def create_image(stream, nx=200, ny=100):
    stream.write("P6\n")
    stream.write("{0} {1}\n".format(nx,ny))
    stream.write("255\n")

    world = HitableList([
        Sphere(np.array([0, 0, -1]), 0.5, Lambertian(np.array([0.8, 0.3, 0.3]))),
        Sphere(np.array([0, -100.5,-1]), 100, Lambertian(np.array([0.8, 0.8, 0.0]))),
        Sphere(np.array([1, 0,-1]), 0.5, Metal(np.array([0.8, 0.6, 0.2]), 1.0)),
        Sphere(np.array([-1, 0,-1]), 0.5, Metal(np.array([0.8, 0.8, 0.8]), 0.3)),
        #Sphere(np.array([ 0.5, -1.5,-3]), 0.5),
        #Sphere(np.array([ -0.5, -0.5, -1]), 1.0),
        #Sphere(np.array([0,-100,5]), 100.5),

    ])

    cam = Camera()
    ns = 25
    color_function = color
    np_array = np.array
    random_function = random.random
    get_ray_function = cam.get_ray
    np_sqrt = np.sqrt
    bytes = []

    for ty in range(ny, 0, -1):
        y = ty - 1 #
        for x in range(nx):
            col = np_array([0., 0., 0.])
            for s in range(ns):
                u = (1.0*x+random_function())/nx
                v = (1.0*y+random_function())/ny
                r = get_ray_function(u, v)
                col += color_function(r, world, 0)
            col /= (1.0)*ns
            col = np_array([np_sqrt(col[0]), np_sqrt(col[1]), np_sqrt(col[2]),])

            ir = int(255.99*col[0])
            ig = int(255.99*col[1])
            ib = int(255.99*col[2])

            bytes.append(bytearray([ir, ig, ib]))

    for b in bytes:
        stream.write(b)
        #stream.write('\n')


if __name__ == '__main__':
    sys.setrecursionlimit(9999999)

    f = open('f1.ppm', "wb")

    import time
    start = time.clock()

    create_image(f, 200, 100)
    print time.clock() - start

    f.close()

    im = Image.open("f1.ppm")
    print im.format, im.size, im.mode
    im.save("f1.png")
    im.save("f1.jpg")
    im.close()
