#Asteroid Orbital Model

from numpy import*
import matplotlib.pyplot as plt
from numpy import linalg
import de421
from jplephem import Ephemeris
 
from visual import *

AU = 149597870.7 # km/AU
epsilon = radians(23.43687)   # obliquity of the Ecliptic
k = 0.01720209895 # Gaussian gravitational constant
R_e = 6378 / AU

JD=[2458313.604632,2458316.68749,2458329.735759]

RA_measured=array([[20, 56, 08.762],[20, 59, 20.410],[21,8,43.156]])

Dec_measured=array([[16,29,32.52],[18, 14, 10.95],[22, 28, 36.27]])

t1=JD[0]
t2=JD[1]
t3=JD[2]
list_test=[JD[0],JD[2]]

tau2=0
tau1=k*(t1-t2)
tau3=k*(t3-t2)

R=[]
for i in JD:
    eph = Ephemeris(de421)
    barycenter = eph.position('earthmoon', i)
    moonvector = eph.position('moon', i)
    earth = barycenter - moonvector * eph.earth_share
    R0 = vector(earth)/AU # This is the Sun to Earth center vector in Equatorial system
 
    R_geocentric = -1.*R0
    R.append(array(R_geocentric))
print 'R_geocentric', R
print

def RAconverter(hrs,minutes,seconds):
    return (hrs/24.+minutes/(24.*60)+seconds/(24.*3600))*360.
def Dec_converter(degrees,minutes,seconds):
    return(degrees+minutes/60.+seconds/3600.)

RA_decimal=[]
for i in range(len(RA_measured)):
    x=RA_measured[i][0]
    y=RA_measured[i][1]
    z=RA_measured[i][2]
    RA_decimal.append(RAconverter(x,y,z))

Dec_decimal=[]
for i in range(len(Dec_measured)):
    x=Dec_measured[i][0]
    y=Dec_measured[i][1]
    z=Dec_measured[i][2]
    Dec_decimal.append(Dec_converter(x,y,z))


print "RA", RA_decimal
print
print "Dec", Dec_decimal
print

def get_rohat(RA,Dec):
    x=(math.cos(math.radians(Dec))*math.cos(math.radians(RA)))
    y=(math.cos(math.radians(Dec))*math.sin(math.radians(RA)))
    z=(math.sin(math.radians(Dec)))
    ro_unit=vector(x,y,z)
    return ro_unit

ro_hat1=get_rohat(RA_decimal[JD.index(t1)],Dec_decimal[JD.index(t1)])
print "rohat1",ro_hat1
print
ro_hat2=get_rohat(RA_decimal[JD.index(t2)],Dec_decimal[JD.index(t2)])
print "rohat2",ro_hat2
print
ro_hat3=get_rohat(RA_decimal[JD.index(t3)],Dec_decimal[JD.index(t3)])
print "rohat3",ro_hat3
print


mu=1.0

m1=.005
m2=0

def acceleration(position):
            return (-1.0*mu*position)/(mag(position)**3)

combined_chi=5
best_chi=combined_chi
r2=vector(0.90188324,0.36779556,0.16562503)
v2=vector(.6,.7,.8)



while combined_chi>1E-6:

    RA_Line=[]
    Dec_Line=[]
    sum_res_RA_squared=0
    sum_res_Dec_squared=0

    for t_test in list_test:
        

        i=0

        while i<=10:
            #f and g series
            f1=1-((1/(2.*((mag(r2))**3)))*(tau1**2))+(((dot(r2,v2))/(2.*(mag(r2))**5))*tau1**3)
            
            g1=tau1-((1./(6*(mag(r2))**3))*tau1**3)
            
            f3=1-((1/(2.*((mag(r2))**3)))*(tau3**2))+(((dot(r2,v2))/(2.*(mag(r2))**5))*tau3**3)
            
            g3=tau3-((1./(6*(mag(r2))**3))*tau3**3)
            
            a1=float(g3/(f1*g3-f3*g1))
            a3=float(-g1/(f1*g3-f3*g1))

            #Get the "Ro"s

            Ro1=(a1*dot(cross(R[JD.index(t1)],ro_hat2),ro_hat3)-dot(cross(R[JD.index(t2)],ro_hat2),ro_hat3)+a3*dot(cross(R[JD.index(t3)],ro_hat2),ro_hat3))/(a1*dot(cross(ro_hat1,ro_hat2),ro_hat3))

            Ro2=(a1*dot(cross(ro_hat1,R[JD.index(t1)]),ro_hat3)-dot(cross(ro_hat1,R[JD.index(t2)]),ro_hat3)+a3*dot(cross(ro_hat1,R[JD.index(t3)]),ro_hat3))/(-1.*dot(cross(ro_hat1,ro_hat2),ro_hat3))

            Ro3=(a1*dot(cross(ro_hat2,R[JD.index(t1)]),ro_hat1)-dot(cross(ro_hat2,R[JD.index(t2)]),ro_hat1)+a3*dot(cross(ro_hat2,R[JD.index(t3)]),ro_hat1))/(a3*dot(cross(ro_hat2,ro_hat3),ro_hat1))

            #Get the "r"s
            r1matrix=array(Ro1*ro_hat1)-R[JD.index(t1)]
            r1=vector(r1matrix[0],r1matrix[1],r1matrix[2])
            
            r2matrix=array(Ro2*ro_hat2)-R[JD.index(t2)]
            r2=vector(r2matrix[0],r2matrix[1],r2matrix[2])

            r3matrix=array(Ro3*ro_hat3)-R[JD.index(t3)]
            r3=vector(r3matrix[0],r3matrix[1],r3matrix[2])

            i+=1
            
        best_r2=r2

        v2=(f3/(g1*f3-g3*f1))*r1-(f1/(g1*f3-g3*f1))*r3
        
        best_v2=v2

        #Verlet Integrator
    
        tau=(t_test-t2)*k 

        f=1-((1/(2.*((mag(r2))**3)))*(tau**2))+(((dot(r2,v2))/(2.*(mag(r2))**5))*tau**3)
  
        g=tau-((1./(6*(mag(r2))**3))*tau**3)

        new_r=f*r2+g*v2

        Ro=R[JD.index(t_test)]+new_r

        #Solving for the Ro Vector

        Ro_Vector=vector(Ro[0],Ro[1],Ro[2])
        MagRo=mag(Ro_Vector)
        
        Rohat=[]
        for n in Ro:
            Rohat.append(float(n/MagRo))
            #Getting the Ro unit vector
        
        Dec=math.degrees(math.asin((Rohat[2])))
        #Finding declination

        RA=360+math.degrees(math.atan(((Rohat[1])/Rohat[0])))
        #Finding RA

        RA_Line.append(RA)
        Dec_Line.append(Dec)


        Residual_RA=(RA_decimal[JD.index(t_test)]-RA)
        Residual_Dec=(Dec_decimal[JD.index(t_test)]-Dec)
    
    sum_res_RA_squared+=(Residual_RA)**2
            
    sum_res_Dec_squared+=(Residual_Dec)**2
    
    chi_squared_RA=float(sum_res_RA_squared)
    chi_squared_Dec=float(sum_res_Dec_squared)
    
    combined_chi=float(sqrt(chi_squared_RA**2+chi_squared_Dec**2))

    #Optimization Code
    r2_prev=new_r
    v2_prev=v2
    sigma=.1
    
    if combined_chi<best_chi:
        best_chi=combined_chi
        best_r2=r2_prev
        best_v2=v2_prev
        sigma=sigma/10.
        print "RESIDUALS of RA, Dec:", Residual_RA, Residual_Dec
        print "BEST CHI_SQUARES:", best_chi
        print "Optimized position:", best_r2
        print "Optimized velocity:", best_v2
        
    else:
        best_chi=best_chi
        best_r2=best_r2
        best_v2=best_v2

    r2=best_r2+vector(random.normal(0.0,sigma,3))
    v2=best_v2+vector(random.normal(0.0,sigma,3))

scene = display(height=600, width=600, range=1.75, up=(0,0,1), forward=(-1,-1,-1))

x_axis=arrow(pos=[0,0,0], axis=(1,0,0),shaftwidth=0.00001, color=color.red)

y_axis=arrow(pos=[0,0,0], axis=(0,1,0),shaftwidth=0.00001, color=color.green)

z_axis=arrow(pos=[0,0,0], axis=(0,0,1), shaftwidth=0.0001, color=color.blue)

sun=sphere(position=[0,0,0],radius=.1, color=color.yellow)
earth=sphere(position=R[1], radius=.2, color=color.blue)
asteroid=sphere(position=best_r2,radius=.05, color=color.blue)
orbit=curve(color=color.blue)
ecliptic=box(pos=(0,0,0), size=(20,20,0.0001), opacity=0.1, color=color.blue)

r_add1=0
tau=0.0
dtau=0.001

m1=.005
m2=0
        
r_n_sub1=best_r2

r=r_n_sub1+v2*dtau+.5*acceleration(r_n_sub1)*(dtau**2)

while True:
    rate(.1)
    accel_next=acceleration(r)
    r_add1=2.*r-r_n_sub1+accel_next*(dtau**2)
    asteroid.pos=r_add1
    orbit.append(pos=r_add1)
    temp=r
    r=r_add1
    r_n_sub1=temp
    tau+=dtau