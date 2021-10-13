# -*- coding: utf-8 -*-
"""

# -*- coding: utf-8 -*-
""""""""""""
To Estimate Rocket Launches, Static and Changing Orbits on different planets in the solar system

Includes:

Classic Orbital Elements Calculator
Hoffman Transfer Calculator (including Interplanetary Transfers)
Rocket Equation

Created to solve tutorial and examination problems in a satellite systems course.

Created on Fri Oct  9 01:00:23 2020

@author: Steven Watts
"""
import numpy as np
import math
import sig as s

# useful constants
aGEO = 42164.157145331665 # earth geostationary orbit in kilometers from earth center
g0 = 9.81 # acceleration due to gravity on earth
R0 = 8314.41 # I don't remember what this is
sidereal = 86164.1 # seconds (sidereal day on earth)

Sig = 4 # number of significant figures for final outputs

# Useful functions - how good are these comments
def ssig(x):

    if isinstance(x,list):
        output = []
        for i in range(len(x)):
            output.append(s.sig(x[i],Sig))
        return output

    else:
        return s.sig(x,Sig)
def Rad(angle):
    return angle * np.pi / 180
def Theta(angle):
    return angle * 180 / np.pi
def Omega(a):
    return np.sqrt(muPlanet['sun'] / (a ** 3))

# The basis of everything
planet = input("Planet ? = ")


planetArray = ['earth', 'moon', 'mars', 'jupiter', 'mercury', 'venus', 'saturn', 'uranus', 'neptune', 'pluto', 'ceres', 'eris', 'sun' ]

AU = 149597871 # astronomical unit in kilometers
year = 365*24*3600 # year in seconds

muArray = [398600, 4903, 42828, 126686000, 22032, 324859, 37931187, 5793939, 6836529, 871, 63, 1108, 132712000000]
radiusArray = [6378, 1737, 3396, 71492, 2439, 6052, 58232, 25362, 24622, 1188, 473, 1163, 696000]
aArray = [1*AU, 384748, 1.5241*AU, 5.2040*AU, 0.3871*AU, 0.7226*AU, 9.5388*AU, 19.1914*AU, 30.0611*AU, 39.5294*AU, 2.76596*AU, 67.6681*AU]
periodArray = [1*year, 2360621, 1.8809*year, 11.862*year, 0.2408*year, 0.6152*year, 29.458*year, 84.01*year, 164.79*year, 248.54*year, 4.599*year, 557*year]

mu,pRadius,aPlanet,TPlanet = np.zeros_like(planetArray),np.zeros_like(planetArray),np.zeros_like(planetArray),np.zeros_like(planetArray)
muPlanet, pRadius,aPlanet,TPlanet = dict(),dict(),dict(),dict()

for i in range(len(planetArray)):
    muPlanet[planetArray[i]] = muArray[i] #km^3/s^2
    pRadius[planetArray[i]] = radiusArray[i] #km

for i in range(len(aArray)):
    aPlanet[planetArray[i]] = aArray[i] #km
    TPlanet[planetArray[i]] = periodArray[i] #s

def mass(planet):
    planet = planet.lower()
    if planet == "earth":
        return 5.974e24
    if planet == "moon":
        return 73.48e21
    if planet == "mars":
        return 641.9e21
    if planet == "jupiter":
        return 1.899e27
    if planet == "mercury":
        return 3.285e23
    if planet == "venus":
        return 4.867e24
    if planet == "saturn":
        return 5.683e26
    if planet == "uranus":
        return 8.681e25
    if planet == "neptune":
        return 1.024e26
    if planet == "pluto":
        return 1.30900e22
    if planet == "ceres":
        return 9.39e20
    if planet == "eris":
        return 1.647e22
    if planet == "sun":
        return 1.989e30
    else:
        print("Unlisted Planet (For Mass).")

if planet != 'sun':
    a0 = aPlanet[planet]
    T0 = TPlanet[planet]
    omega0 = Omega(a0)

mu = muPlanet[planet]
r0 = pRadius[planet]
m0 = mass(planet)

def velocity(r,epsilon,mu=mu):
    return np.sqrt(2) * np.sqrt(mu/r + epsilon)

""""""""""
COE Calculator for given semi-major axis or other quantity

epsilon = -mu/2a
epsilon = v^2/2 - mu/r
v = sqrt(2*(mu/r + epsilon))

Ellipse:
a = (r_p + r_a)/2 = h^2/mu *(1-e^2)^(-1)
e = (r_a - r_p)/(r_a + r_p)
r = a*(1-e^2)/(1+ecos(theta))
costheta = (a*(1-e^2)/r -1)/e


Circular equations:
v_circ = sqrt(mu/r)
T = 2*pi*r/sqrt(mu/r) = 2*pi*sqrt(r**3/mu)

Parabolic:
epsilon=0
v^2/2 - mu/r =0
-> v = sqrt(2mu/r)
v_esc = sqrt(2)*v_circ

Hyperbolic:
v^2/2 - mu/r = mu/2a
as r --> inf
v_inf = sqrt(mu/a)
- > v^2/2 - mu/r = v_inf^2/2
- > v^2 = v_esc^2 + v_inf^2
"""""""""""""""
def tconv(T):

    #works for negative values
    T = abs(T)

    second = ssig(T % 60)

    # checks to see if T duration is more than 1 year
    if T >= (3600 * 24 * 365):
        year = int(T // (3600 * 24 * 365))
        day = int((T / (3600 * 24 * 365) % 1) * 365 // 1)
        hour = int((T / (3600 * 24) % 1) * 24 // 1)
        minute = int(T // 60 % 60)
        return str(year) + ' year(s) ' + str(day) + ' days, ' + str(hour) + ' hours, ' + str(minute) + \
            ' minutes & ' + str(second) + ' seconds'

    # checks to see if T duration is more than a day
    if T >= (3600 * 24) and T <(3600 * 24 * 365):
        day = int(T // (3600*24))
        hour = int((T / (3600*24) % 1 ) * 24 // 1)
        minute = int(T // 60 % 60)
        return str(day) + ' day(s), ' + str(hour) + ' hours, ' + str(minute) + \
                                  ' minutes & ' + str(second) + ' seconds'

    if T >= 3600 and T < (3600 * 24): # for T greater than 1 hour
        hour = int(T // 3600 % 60)
        minute = int(T // 60 % 60)
        return str(hour) + ' hour(s), ' + str(minute) + \
                                               ' minutes & ' + str(second) + ' seconds'

    if T >= 60 and T < 3600: # greater than 1 minute
        minute = int(T // 60 % 60)
        return str(minute) + ' minute(s) & ' + str(second) + ' seconds'

    else:
        return str(second) + ' seconds'

def perturbation(a,i,e,deltaOmega=None,period=None):

    if a < r0:
        a+=r0

    i = i*np.pi/180

    omegaDot = s.sig((-2.06474*10**14)*a**(-7/2)*math.cos(i)*(1-e**2)**(-2),5)

    messageDefault = "Ascending Node Decay = "
    
    if omegaDot > 0:
        direction = "eastward"
    else:
        direction = "westward"
    
    #default == day
    if period =="second":
        print(messageDefault,s.sig(omegaDot/(24*3600),5),"deg/s",direction)
        return s.sig(omegaDot/(24*3600),5)
    
    if period =="minute":
        omegaDot = s.sig(omegaDot/(24*60),5)
        print(messageDefault,omegaDot,"deg/min",direction)
    
    if period =="hour":
        omegaDot = s.sig(omegaDot/24,5)
        print(messageDefault,omegaDot,"deg/hr",direction)
    if period == "week":
        omegaDot = s.sig(omegaDot*7,5)
        print(messageDefault,omegaDot,"deg/week",direction)
    if period == "year":
        omegaDot = s.sig(omegaDot*365,5)
        print(messageDefault,omegaDot,"deg/yr",direction)
    else:
        print(messageDefault,omegaDot,"deg/day",direction)


    if isinstance(deltaOmega,float) or isinstance(deltaOmega,int):
        t = s.sig(abs(deltaOmega/omegaDot),5)
        print('Time to move ',deltaOmega,' degrees: ',t, ' days')

    print('-------------------------------------------\n')
    return omegaDot

def inclination(a,e,RAAN):
    if a < r0:
        a+=r0
    i1 = 1/(-2.06474*10**14)
    i2 = a**(7/2)*(1 - e ** 2)**(2)
    i = math.acos(i1*i2*RAAN)
    i = ssig(Theta(i))
    print('Angle of Inclination for a RAAN of ',RAAN,': ',i,' degrees')
    return i

def hyperbola(r_p,v_p,trueAnomaly=None):
    if r_p < r0:
        r_p +=r0
    #epsilon = -mu / (2 * a)
    epsilon = v_p**2/2 - mu/r_p

    a = abs(-mu/(2*epsilon))

    v_inf = np.sqrt(mu/a)

    e = 1 - r_p/(-a)


    a,epsilon,v_inf,e = sigg(a),sigg(epsilon),sigg(v_inf),sigg(e)

    print('\n---------------------------------')
    print('Semimajor Axis: ',a, 'km')
    print('Hyerbolic Excess Velocity: ',v_inf,' km/s')
    print('Eccentricty: ',e)

    if isinstance(trueAnomaly,float) or isinstance(trueAnomaly,int):
        trueAnomaly = trueAnomaly*np.pi/180 #radians conversion

        r = a*(e**2-1)/(1+e*np.cos(trueAnomaly))
        v_r = np.sqrt(2*(mu/r + epsilon))
        h = np.sqrt(r*mu*(1+e*np.cos(trueAnomaly)))
        v_aizmuth = h/r

        r,v_r,h,v_aizmuth = sig(r),sig(v_r),sig(h),sig(v_aizmuth)

        print('Radius at True Anomaly of: ',trueAnomaly,' Radians: ',r, 'km')
        print('Radial Velocity at Subject True Anomaly: ',v_r)
        print('Angular Momentum at Subject True Anomaly: ',h)
        print('Aizmuth Velocity at Subject True Anomaly: ',v_aizmuth)
        print('\n')

        return a,v_inf,e,trueAnomaly,v_r,h,v_aizmuth
    else:
        return a, v_inf, e


def orbit(r_a,r_p,vP=None,vA=None,T=None,anomalyDistance=None,trueAnomaly=None,parabolic='N'):

    #auto detects if altitude or not
    if isinstance(r_a,float) or isinstance(r_a,int):
        if r_a < r0:
            r_a = r_a+r0
    if isinstance(r_p,float) or isinstance(r_p,int):
        if r_p < r0:
            r_p = r_p+r0

    # Executes a chain of alternate equations only if r_a and r_p not defined.
    a=0
    e=0
    if isinstance(r_a,str) or isinstance(r_p,str):

        if isinstance(trueAnomaly,list) and isinstance(anomalyDistance,list):

            for r in range(len(anomalyDistance)):
                if anomalyDistance[r] < r0:
                    anomalyDistance[r] +=r0

            r1,r2 = anomalyDistance[0],anomalyDistance[1]
            theta1,theta2 = np.pi*trueAnomaly[0]/180,trueAnomaly[1]*np.pi/180


            e = (r1/r2 - 1)*(np.cos(theta2) - r1/r2*np.cos(theta1))**(-1)
            a = r1*(1+e*np.cos(theta1))/(1-e**2)
            r_p = a*(1-e)
            r_a = 2*a - r_p

        if isinstance(T,int) or isinstance(T,float):
            if T!=0:
                a = ((T/(2*np.pi))**2 *mu)**(1/3)
                if isinstance(r_a, float) or isinstance(r_a, int):
                    r_p = 2*a - r_a
                if isinstance(r_p, float) or isinstance(r_p, int):
                    r_a = 2*a - r_p

                if not isinstance(vP,int) and not isinstance(vA,int):
                    if isinstance(r_a,str) and isinstance(r_p,str):
                        print('a: ',a, 'km')
                        return a

        if isinstance(vA,int) or isinstance(vA,float):
            if vA !=0:
                if a==0:
                    a = (mu * r_a) / (-vA ** 2 * r_a + 2 * mu)
                else:
                    r_a = (2 * mu * a) / (vA ** 2 * a + mu)
                r_p = 2*a - r_a

        if isinstance(vP,int) or isinstance(vP,float):
            if vP !=0:
                if a==0:
                    a = (mu * r_p) / (-vP ** 2 * r_p + 2 * mu)
                else:
                    r_p = (2*mu*a)/(vP**2*a + mu)
                r_a = 2*a - r_p

    if a ==0:
        a = (r_p + r_a) / 2
    if e==0:
        e = (r_a - r_p)/(r_a+r_p)
    if not isinstance(T,float):
        T = 2*np.pi*np.sqrt(a**3/mu)

    epsilon = -mu / (2 * a)

    v = velocity(a,epsilon)
    vA = velocity(r_a,epsilon)
    vP = velocity(r_p,epsilon)

    h = np.sqrt(-mu**2/(2*epsilon)*(1-e**2))


    def sig(x):
        return s.sig(x,5)

    print('--------------------------------------------------')

    if isinstance(anomalyDistance,int) == True or isinstance(anomalyDistance,float) == True:
        if anomalyDistance !=0:
            #costheta = (a * (1 - e ^ 2) / r - 1) / e
            if anomalyDistance < r0:
                r = anomalyDistance+r0
            trueAnomaly = math.acos((a*(1-e**2))/(e*r) - 1/e)
            v_r = mu / h * np.sin(trueAnomaly)
            #conversion
            trueAnomaly = 180/np.pi * trueAnomaly
            v_a = h / r
            v_r, v_a, trueAnomaly = sig(v_r), sig(v_a), sig(trueAnomaly)
            print('True Anomaly at ',anomalyDistance,'km: ',trueAnomaly, ' degrees')
            print('Azimuth Velocity: ', v_a, 'km/s')
            print('Angular Velocity: ', v_r, 'km/s')

    parabolic = parabolic.lower()
    if parabolic == 'y':
        v_parabolic = np.sqrt(2)*v
        v_excess = v*(np.sqrt(2)-1)

        v_parabolic,v_excess = sig(v_parabolic),sig(v_excess)
        print('Velocity for Parabolic Trajectory: ',v_parabolic, 'km/s')
        print('Extra Velocity needed for Parabolic Trajectory: ',v_excess,'km/s')

    a,e,v,vA,vP,epsilon,h,r_a,r_p = sig(a),sig(e),sig(v),sig(vA),sig(vP),sig(epsilon),sig(h),sig(r_a),sig(r_p)

    if r_a != r_p:
        print('Semimajor Axis: ',a, 'km')
        print('Average Velocity: ',v, 'km/s')
        print('Apo Velocity: ', vA, 'km/s' )
        print('Peri Velocity: ',vP, 'km/s')
        print('Radius at Apogee: ', r_a)
        print('Radius at Perigee: ',r_p)
    else:
        print('Velocity: ',v, 'Kilometers / second')




    print('Specific Angular Momentum: ',h)

    print('Specific energy: ',epsilon)
    print('Eccentricty: ',e)
    print('T: ('+str(T)+') ',tconv(T))
    print("")
    print("Data accurate to "+str(Sig)+" significant figures.")
    print('--------------------------------------------------')
    return a,e,v,vA,vP,T

def coe(r,v,planet):
    
    planet = planet.lower()
    
    # finds mu from dictionary
    mu = muArray[planet]
    
    R = np.sum(r**2)**(1/2)
    V = np.sum(v**2)**(1/2)
    
    energy = V**2/2 - mu/R
    
    a = s.sig(-mu/(2*energy),5)
    
    e = 1/mu*((V**2 - mu/R)*r - np.dot(r,v)*v)
    
    
    E = s.sig(np.sum(e**2)**(1/2),5)
    
    h = np.cross(r,v)
    H = np.sum(h**2)**(1/2)
    
    i = s.sig(math.acos(h[2]/H)*180/np.pi,5)
    
    
    n = np.cross(np.array([0,0,1]),h)
    N = np.sum(n**2)**(1/2)
    
  
    RAAN = s.sig(math.acos(n[0]/N)*180/np.pi,5)
    if n[1] < 0:
        RAAN = s.sig(360 -RAAN,4)
        
    perigee = s.sig(math.acos(np.dot(n,e)/(N*E))*180/np.pi,5)
    if e[2] < 0:
        perigee = s.sig(360 -perigee,5)
    
    trueAnomaly = s.sig(math.acos(np.dot(e,r)/(E*R))*180/np.pi,5)
    boundState = "OUT BOUND"
    if np.dot(r,v) < 0:
        trueAnomaly = s.sig(360 -trueAnomaly,5)
        boundState = "IN BOUND"
    # Rotational Period
    T = 2*np.pi*np.sqrt(a**3/mu)

    # Orbit Analysis
    GEO,LEO,SSO = 0,0,0
    
    orbitType = ""

    ######## EARTH ORBITS #############
    if planet == "earth":
        if a <= 1000+6378:
            LEO = 1
            if i > 90 and i <= 105:
                if E < 0.4:
                    SSO = 1
        
        if T >= 86164.1/3600 -.2 and T <= 86164.1/3600 + .2:
            if i == 0:
                GEO = 1
            else:
                GEO = 2
                
                
        if i > 90:
            orbitType = "RETROGRADE"
            
        if i == 0:
            orbitType = "EQUATORIAL"
        
        if i == 90:
            orbitType = "POLAR"
        
        if i >= 63.43 - 1 and i <= 63.43 + 1:
            if T >= 43082-500 and T <= 43082+500:
                orbitType = "MOLNIYA"
            

        
        if LEO == 1:
            orbitType = "LEO" + " " + orbitType
            
        if SSO == 1:
            orbitType = "SUN-SYNCHRONOUS (SSO)"
            
        if GEO == 1:
            orbitType = "GEOSTATIONARY"
        if GEO == 2:
            orbitType = "GEOSYNCHRONOUS"

    
    ##################################
    
    if planet =="sun":
        orbitType = "SOLAR"
        
    if planet == "moon":
        orbitType = "LUNAR"
        
    if planet == "mars":
        orbitType = "MARTIAN"
        
    if planet == "jupiter":
        orbitType = "JOVIAN"
    

    print("")
    print(" --- 6 COEs --- ")
    print("semi-major axis = ",a, " kilometers")
    print("eccentricity = ",E)
    print("inclination = ",i," degrees")
    print("RAAN = ",RAAN," degrees")
    print("Argument of Perigee = ",perigee," degrees")
    print("True Anomaly = ",trueAnomaly, " degrees")
    print("Orbit Type: ",orbitType)
    print("Satellite is: ",boundState)
    print("")
    print(" --- More Information --- ")
    print("")
    print("Angular Momentum Vector = ",h)
    print("Angular Momentum = ",H)
    print("Ascending Node Vector = ",n)
    print("Ascending Node = ",N)
    print("eccentricity Vector = ",e)
    print("Radius = ",R)
    print("Velocity = ",V)
    if T < 3600:
        tconv(T,'minute')
    else:
        tconv(T,'hour')
    print("")
    print("Data accurate to 5 significant figures.")
    
    
    return a,E,i,RAAN,perigee,trueAnomaly,h,H,n,N,e,R,V,T

def Transfer(h1a, h1p, h2a, h2p, i,f,theta=None,launch='apogee',impulse=None,massInt=None,exhaust=None,massBurn=None):

    r1a = h1a + r0
    r1p = h1p + r0
    r2a = h2a + r0
    r2p = h2p + r0

    a1 = (r1a + r1p) / 2
    eps1 = -mu / (2 * a1)
    v1 = velocity(r1p, eps1)
    omega1 = np.sqrt(mu / a1 ** 3)

    # radians conversion
    i = np.pi / 180 * i
    f = np.pi / 180 * f

    delta = i - f

    print('\n'+'------------------------------------------------------')
    # checks to see if final altitude is the same. If yes, performs simple plane change.
    if h1a == h2a and h1p == h2p:
        deltaV2 = 2 * v1 * np.sin(delta / 2)

        # for calc outside if statement
        vf = v1
        deltaV1 = v1
        vi = v1

        deltaV = ssig(deltaV2)



        ### CO-ORBITAL RENDEZVOUS
        if isinstance(theta, float) or isinstance(theta, int):
            print('----RENDEZVOUS----')
            theta = theta * np.pi / 180 #rad conversion
            thetaTravel = 2 * np.pi - theta

            TOF = thetaTravel / omega1

            aPhase = ( mu * ( thetaTravel / (2 * np.pi * omega1) ) ** 2 ) ** (1/3)
            aTransfer = (aPhase + a1) / 2
            epsTransfer = - mu / (2 * aTransfer)
            deltaV1 = velocity(a1, epsTransfer)
            deltaV2 = velocity(aPhase,epsTransfer)
            deltaV = abs(deltaV1 - deltaV2)

            # formatting
            aPhase,aTransfer,epsTransfer = ssig(aPhase),ssig(aTransfer),ssig(epsTransfer)

            print('TOF = (' + str(TOF) + ') ' + tconv(TOF))
            print('Semi-Major Axis of Phasing Orbit: ',aPhase, 'km')
            print('Semi-Major Axis of Transfer Orbit: ',aTransfer, 'km')
            print('Specific Energy of Transfer Orbit: ',epsTransfer, 'km^2/s^2')

        # formatting
        deltaV,v1 = ssig(deltaV),ssig(v1)
        print('DeltaV: ', deltaV, ' km/sec')
        print('Velocity of Main Orbit: ',v1, 'km/sec')



    else:
        a2 = (r2a + r2p) / 2
        aTransfer = (r1p + r2a) / 2
        eps2 = -mu / (2*a2)
        if r1a > r2a:
            epsTransfer = -mu / (r1a + r2p)
        else:    
            epsTransfer = - mu / (2*aTransfer)

        v1Transfer = velocity(r1p, epsTransfer)
        v2Transfer = velocity(r2a,epsTransfer)

        omega2 = np.sqrt(mu / a2 ** 3)
        v2 = velocity(r2a,eps2)

        if i == 0 and f == 0:
            print('Hohmann Transfer\n')
            deltaV1 = abs(v1Transfer- v1)
            deltaV2 = abs(v2Transfer - v2)

            deltaV = deltaV1 + deltaV2

            TOF = np.pi * np.sqrt(aTransfer ** 3 / mu)



            alpha = omega2 * TOF  # lead angle

            thetaF = np.pi - alpha  # final angle


            ### CO-PLANAR RENDEZVOUS
            if isinstance(theta, float) or isinstance(theta, int):
                thetaRad = theta * np.pi / 180  # rad conversion

                alpha = omega2 * TOF  # lead angle
                alphaDeg = ssig(alpha * 180 / np.pi)  # degree conversion
                
                thetaF = np.pi - alpha  # final angle
                thetaF_deg = alpha * 180 / np.pi

                waitTime = (thetaF - thetaRad) / (omega2 - omega1)

                alphaDeg, thetaF_deg, waitTime = ssig(alphaDeg), ssig(thetaF_deg), ssig(waitTime)

                print('----RENDEZVOUS----')
                print('Initial Angle: ', theta, 'deg', ' / ', thetaRad, 'rad')
                print('Final Angle: ', thetaF_deg, 'deg'' / ', thetaF, 'rad')
                print('Lead Angle: ', alphaDeg, 'deg'' / ', alpha, 'rad')
                print('Wait Time: (' + str(waitTime) + ') ' + tconv(waitTime))


        else:

            if launch == 'perigee':
                print('Plane Change at Perigee \n')
                vi = v1
                vf = v1Transfer
                deltaV1 = np.sqrt(vi ** 2 + vf ** 2 - 2 * vi * vf * np.cos(delta))
                deltaV2 = abs(v2 - v2Transfer)

            if launch == 'apogee':
                print('Plane Change at Apogee \n')
                
                # Special Plane Change Type - always go from highest altitude (see Q1,2018. Midterm,2020.)
                if r1a > r2a:
                    vi = velocity(r1a,epsTransfer)
                    vf = velocity(r1a,eps1)
                    v2Transfer = velocity(r2p,epsTransfer)
                    deltaV1 = np.sqrt(vi ** 2 + vf ** 2 - 2 * vi * vf * np.cos(delta))
                    deltaV2 = abs(v2Transfer - v2)
                    
                else:
                    vi = v2Transfer
                    vf = v2
                    deltaV1 = abs(v1Transfer - v1)
                    deltaV2 = np.sqrt(vi ** 2 + vf ** 2 - 2 * vi * vf * np.cos(delta))
                    
                


        eps2,v2,v2Transfer,omega2 = ssig(eps2),ssig(v2),ssig(v2Transfer),ssig(omega2)

        deltaV = deltaV1 + deltaV2

        # formatting
        deltaV1,deltaV2,deltaV,v1Transfer,epsTransfer = \
            ssig(deltaV1),ssig(deltaV2),ssig(deltaV),ssig(v1Transfer), \
            ssig(epsTransfer)

        print('Delta V1: ', deltaV1, ' km/sec')
        print('Delta V2: ', deltaV2, ' km/sec')
        print('Total DeltaV: ', deltaV, ' km/sec')



        TOF = np.pi * np.sqrt(aTransfer ** 3 / mu)

    if i !=0 or f !=0:
        deltaV2X = vf * np.cos(delta) - vi
        deltaV2Y = vf * np.sin(delta)

        # formatting
        deltaV2X, deltaV2Y, vi, vf = ssig(deltaV2X), ssig(deltaV2Y), ssig(vi), ssig(vf)


    # formatting
    v1,eps1,omega1 = ssig(v1),ssig(eps1),ssig(omega1)





    print('\n')

    if i != 0 or f != 0:
        print('Delta V2 X component: ', deltaV2X, ' km/sec')
        print('Delta V2 Y component: ', deltaV2Y, ' km/sec')
    print('Specific Energy of Orbit 1: ',eps1, 'km^2/s^2')
    if h1a != h2a or h1p != h2p: #and
        print('Specific Energy of Orbit 2: ',eps2, 'km^2/s^2')
        print('TOF = (' + str(TOF) + ') ' + tconv(TOF))
        print('Specific Energy of Transfer Orbit: ',epsTransfer, 'km^2/s^2')
        if i != 0 or f != 0:
            print('Initial Velocity of Plane Change: ', vi, 'km/sec')
            print('Final Velocity of Plane Change: ', vf, 'km/sec')
        print('Velocity of Orbit 1: ',v1, 'km/sec')
        print('Velocity of Orbit 2: ',v2, 'km/sec')
        print('Velocity of Transfer Orbit from Perigee: ',v1Transfer, 'km/sec')
        print('Velocity of Transfer Orbit from Apogee: ', v2Transfer, 'km/sec')
        if h1a != h2a or h1p != h2p: #and
            print('Angular Velocity of High Orbit: ', omega2, 'rad/sec')

    print('Angular Velocity of Low Orbit: ', omega1, 'rad/sec')

    print('------------------------------------------------------'+'\n')



def Impulse(gamma,mol,temperature):

    R = R0/mol

    a0 = np.sqrt(gamma * R * temperature)

    C = a0 / (gamma * (2 / (gamma + 1))**((gamma + 1)/(2 * gamma - 2)) )

    I = C / g0 * gamma * ((2/(gamma - 1)) * (2/(gamma + 1))**((gamma+1)/(gamma-1)) )**(1/2)

    C, I = ssig(C),ssig(I)
    print('Characteristic Exhaust Velocity: ',C, 'km/sec')
    print('Specific Impulse: ('+str(I)+')', tconv(I))

    return R,a0,C,I

def Rocket(deltaV,I=None,massInt=None,massBurn=None,exhaustV=None):

    if massInt == None and massBurn == None:
        return print('Need a Mass Variable to Compute.')


    if isinstance(I,float) or isinstance(I,int):
        C = I*g0
    elif isinstance(exhaustV,float) or isinstance(exhaustV,int):
        C = exhaustV
    else:
        return print('Need an Exhaust Velocity or Specific Impulse quantity.')

    ratio = np.exp(deltaV / C)

    if isinstance(massInt,float) or isinstance(massInt,int):
        massF = (massInt * (ratio - 1) ) / ratio
    else:
        massF = massBurn * (ratio - 1)

    massF = ssig(massF)

    print('Propellant Mass: ',massF, 'kg')
    return massF

def Launch(L0, i,Node='AN',a=None,vLoss=0,phi=0,beta=0,vBurn=0): #NOTE: not designed for non circular and only for earth
    
    L0,i, phi = Rad(L0), Rad(i), Rad(phi)
    alpha = i

    #if i == L0:
    #   phi = np.pi / 2

    gamma = math.asin(math.cos(alpha)/math.cos(L0))

    if beta == 0:
        if L0 < 0:
            if Node == 'AN':
                beta = np.pi - gamma
            if Node == 'DN':
                beta = gamma
        if i > 0:
            if Node == 'AN':
                beta = gamma
            if Node == 'DN':
                beta = np.pi - gamma
    else:
        beta = Rad(beta)
    vLaunch = 0.4651 * np.cos(L0)

    if isinstance(a,float) or isinstance(a,int):
        if a < r0:
            a += r0
        vG = np.sqrt((2 * mu * (a - r0)) / (a*r0)) #zenith

        if vBurn == 0:
            vBurn = np.sqrt(mu/a)
        vBurnS = -vBurn * np.cos(phi) * np.cos(beta)
        vBurnE = vBurn * np.cos(phi) * np.sin(beta)
        vBurnZ = vBurn * np.sin(phi)

        deltaVS = vBurnS
        deltaVE = vBurnE - vLaunch
        deltaVZ = vG + vBurnZ

        deltaV = np.sqrt(deltaVS**2 + deltaVE**2 + deltaVZ**2)

        deltaVDesign = deltaV + vLoss

        # formatting
        vBurnS,vBurnE,vBurnZ,deltaVS,deltaVE,deltaVZ,deltaV,deltaVDesign,vG, = \
        ssig(vBurnS),ssig(vBurnE),ssig(vBurnZ),ssig(deltaVS),ssig(deltaVE), \
        ssig(deltaVZ),ssig(deltaV),ssig(deltaVDesign),ssig(vG)

        print('\n---------------------------------------------------------')
        print('Delta V Needed: ',deltaV, 'km/sec')
        if vLoss > 0:
            print('Delta V Design: ',deltaVDesign, 'km/sec')

        print('Velocity of Gravity Loss: ', vG, 'km/sec (Z)')
        print('Burnout Velocities (S,E,Z): ',vBurnS,vBurnE,vBurnZ, 'km/sec')
        print('deltaV Components (S,E,Z): ',deltaVS,deltaVE,deltaVZ, 'km/sec')

    iDeg, alphaDeg, gammaDeg, betaDeg = ssig(Theta(i)), \
    ssig(Theta(alpha)),ssig(Theta(gamma)),ssig(Theta(beta))

    vLaunch = ssig(vLaunch)
    print('Launch Velocity: ', vLaunch, 'km/sec (Z)')
    print('---------')
    print('Angles:')

    print('Alpha: ', alphaDeg, 'deg')
    print('Gamma: ', gammaDeg, 'deg')
    print('Beta: ', betaDeg, 'deg')

def Interplanetary(target,theta=0,parkPlanet=None,parkTarget=None):

    if target == planet:
        return print('Target can not be same as starting planet.')
    
    if target == 'sun':
        return print('Program not Currently able to compute Sun Orbits.')
    
    # starting formalities
    Planet = planet[0].upper() + planet[1:]
    Target = target[0].upper() + target[1:]
    print('\n---------------------------------------------------------')
    print('INTERPLANETARY TRANSFER FROM ' + Planet.upper() + ' to ' + Target.upper() + '\n')

    muSun = muPlanet['sun']

    
    aTarget = aPlanet[target]
    aTransfer = (aTarget + a0) / 2
    epsTransfer = -muSun/ (2 * aTransfer)
    epsPlanet = -muSun / (2 * a0)
    epsTarget = -muSun / (2 * aTarget)
    omegaTarget = Omega(aTarget)
    TTarget = TPlanet[target]

    vPlanet = velocity(a0,epsPlanet,mu=muSun)
    vTransfer = velocity(a0,epsTransfer,mu=muSun)
    vInf = abs(vPlanet - vTransfer)
    vTarget = velocity(aTarget,epsTarget,mu=muSun)
    vTransferTarget = velocity(aTarget,epsTransfer,mu=muSun)
    vTargetInf = abs(vTarget - vTransferTarget)

    epsInf = (vInf ** 2) / 2
    epsInfTarget = (vTargetInf ** 2) / 2

    if isinstance(parkPlanet, int) or isinstance(parkPlanet, float):
        if parkPlanet < r0:
            parkPlanet += r0
        if parkTarget < pRadius[target]:
            parkTarget += pRadius[target]

        vHyp = velocity(parkPlanet,epsInf)
        vPark = np.sqrt(mu / parkPlanet)
        vBoost = abs(vHyp - vPark)
        vHypTarget = velocity(parkTarget,epsInfTarget,mu=muPlanet[target])
        vParkTarget = np.sqrt(muPlanet[target] / parkTarget)
        vRetro = abs(vParkTarget - vHypTarget)

        vMission = vBoost + vRetro

        # formatting
        vBoost,vRetro,vMission = ssig(vBoost),ssig(vRetro),ssig(vMission)

        print('Delta V Mission: ',vMission, 'km/sec')
        print('Delta V Boost: ',vBoost, 'km/sec')
        print('Delta V Retro: ',vRetro, 'km/sec')
        print('Parking Velocities (Planet, Target): ',vPark, vParkTarget, 'km/sec')
        print('Hyperbolic Excess Velocities (Planet, Target): ',vHyp, vHypTarget, 'km/sec')

    print('Tangential Velocities (Planet, Target: ',vPlanet,vTarget, 'km/sec')
    print('Departure Velocities (Planet, Target): ',vInf,vTargetInf, 'km/sec')
    print('Transfer Velocities (Planet, Target): ',vTransfer,vTransferTarget, 'km/sec')
    print('')
    print('-----')

    SOI = a0 * (m0 / mass('sun')) ** (2 / 5)
    if target == 'moon':
        SOITarget = aTarget * (mass(target) / m0 ) ** (2 / 5)
    else:
        SOITarget = aTarget * (mass(target) / mass('sun')) ** (2 / 5)

    synodic = (2 * np.pi) / abs(omega0 - omegaTarget)

    TOF = np.pi * np.sqrt(aTransfer ** 3 / muPlanet['sun'])
    alpha = omegaTarget * TOF
    alphaDeg = Theta(alpha)
    phase = np.pi - alpha
    phaseDeg = Theta(phase)

    thetaRad = Rad(theta)
    waitTime = abs(phase - thetaRad) / abs(omegaTarget - omega0)

    #formatting
    thetaRad,waitTime = ssig(thetaRad),ssig(waitTime)

    # formatting
    SOI,SOITarget,synodic,TOF = ssig(SOI),ssig(SOITarget),ssig(synodic),ssig(TOF)
    alpha,alphaDeg,phase,phaseDeg = ssig(alpha),ssig(alphaDeg),ssig(phase),ssig(phaseDeg)
    aTransfer, epsTransfer = ssig(aTransfer),ssig(epsTransfer)

    print('Sphere of Influence of ' + Target + ': ',SOITarget, 'km')
    print('Sphere of Influence of ' + Planet + ': ', SOI, 'km')
    print('Semi-Major Axis of Transfer Orbit: ',aTransfer, 'km')
    print('Specific Energy of Transfer Orbit: ',epsTransfer, 'km^2/s^2')
    print('Lead Angle: ',alphaDeg, 'deg / ',alpha, 'rad')
    print('Inital Angle: ', theta, 'deg / ', thetaRad, 'rad')
    print('Final Angle: ', phaseDeg, 'deg / ', phase, 'rad')
    print('TOF (' + str(TOF) + ')',tconv(TOF))
    if isinstance(theta, float) or isinstance(theta, int):
        print('Wait Time (' + str(waitTime) + ')', tconv(waitTime))

    print('Synodic Period (' + str(synodic) + ')',tconv(synodic))
    print('----')
    print('Specific Energies (Planet, Target): ',epsPlanet,epsTarget, 'km^2/s^2')
    print('Angular Velocities (Planet, Target): ',omega0,omegaTarget, 'rad/sec')
    print('Semi Major Axes (Planet, Target): ',a0,aTarget, 'km')


    print('---------------------------------------------------------')
    return SOI, SOITarget,vBoost,vMission

runcheck = input("Run COE Program (Y/N) ? ")
if runcheck == "Y" or runcheck == "y":
    strR = input("radius ? = ")
    strV = input("velocity ? = ")

    r = np.array([float(strR.split()[0]), float(strR.split()[1]), float(strR.split()[2])])
    v = np.array([float(strV.split()[0]), float(strV.split()[1]), float(strV.split()[2])])
    coe(r, v, planet)

misc = False
if misc == True:
    #works:
    # orbit(1600,600)
    # orbit(300,250)
    # orbit(42000,250)
    # orbit('?',150,T=5400)
    #orbit(70000,7000,anomalyDistance=1000)
    #orbit('?',640,vP=8)
    #orbit('?','?',vP=8,T=7200)
    #orbit('?','?',vA=6.1825,T=7200)
    #orbit('?','?',trueAnomaly=[130,50],anomalyDistance=[1700,500])
    #hyperbola(300,15,trueAnomaly=100)
    #orbit(1000,1000,parabolic='y')
    #orbit(25000,25000,parabolic='y')
    # orbit('?','?',T=43082)
    pass



    ######################
    #tutorial 4 questions

    #Transfer(300,300,3000,3000,0,0) #- Good
    #Transfer(800,480,16000,16000,0,0) # - delta V1 slightly off ~ 0.1 km/s
    #Transfer(500,500,1000,1000,0,0) #- GOOD
    #Transfer(1.496e8-r0,1.496e8-r0,2.278e8-r0,2.278e8-r0,0,0) #planet == sun - GOOD

    #Transfer(250,250,250,250,28,57) #- Good
    #Transfer(600,600,600,600,28,20) #- Good
    #Transfer(130,130,130,130,57,90) # GOOD
    #Transfer(130,130,130,130,90,125) # GOOD. Negative. TODO: fix negatives, or confirm if should be negative -- need STK
    #Transfer(150,150,20000,20000,45, 28) # GOOD.
    #Transfer(500,500,20469,20341,28.5,55) # GOOD.
    # - MISSING: LEO TO GEO plane change then hohmann
    # - MISSING: LEO TO GEO hohmann then plane change
    #Transfer(300,300,aGEO-r0,aGEO-r0,28.5,0,launch='perigee') # GOOD.
    #Transfer(300,300,aGEO-r0,aGEO-r0,28.5,0) # GOOD.

    #Transfer(192, 192, 35782, 35782, 0, 0, theta=180) #- GOOD
    #Transfer(120, 120, 240, 240, 0, 0, theta=135) #- GOOD (Wait time off by ~ 20 minutes)
    #Transfer(200,200,410,410,0,0,theta=30) #- GOOD

    #Transfer(240,240,240,240,0,0,theta=-35) #- GOOD
    #Transfer(350,350,350,350,0,0,theta=5.1096) # - GOOD
    #Transfer(350,350,350,350,0,0,theta=-5.1096) # - GOOD

    # ( That is all of tutorial 4 ) #


#tutorial 5a questions (perturbations) (and tutorial 7)
#perturbation(700,30,0)
#perturbation(7278,50,0.001,deltaOmega=60)
#inclination(400,0,0.986)





#tutorial 7 questions (propulsion)
    #Rocket(10,exhaustV=1000,massInt=1000) # (Q7) - GOOD
    # Q8 - not sure how to do
    #Impulse(1.225,21.79,3415) # (Q9) GOOD
    #Rocket(1831, I=290,massInt=None, massBurn=907) # (Q10) GOOD
    #Rocket(3880,massBurn=1300,I=340) # (Q10) GOOD

    ################
    #tutorial 8 questions (launch)

    #Launch(28.5,41) # - GOOD
    #Launch(28.5,28.5,a=400,vLoss=1.5) # - GREAT
    #Launch(9.05,9.05,a=500,vLoss=1) # - GREAT
    #Launch(34.6,55,a=500,Node='DN',vLoss=0.8)
    #Transfer(500,500,500,500,55,90,)
    #orbit('?',500,T=sidereal/2)  ###a = 26561.754569239343, r_a = 46246
    #Transfer(500,500,46246-r0,500,55,63.4) # INCORRECT (maybe) - PERFECT
    #Transfer(210,210,aGEO-r0,aGEO-r0,5.2,0) # For vBurn -> = 10.23
    #Launch(5.2,5.2,beta=87,phi=0,a=210,vLoss=0.8,Node='AN',vBurn=10.23) #- Great
    #Launch(-39.2,i=90,a=450,vLoss=1) # First try. I love this code.
    #Rocket(3923, massBurn=850, I=290) - #off by 5 kg
    #Rocket(3880,I=340,massBurn=1300) # - off by 4.5 kg - learnt this was because Naomi used 9.8, not 9.81
    #Transfer(500,500,aGEO-r0,aGEO-r0,0,0) #Vburn = 3.816
    #Rocket(3816,massBurn=1200,I=300) #GOOD. bit hard to understand at first. Q only wanted manoveure

    # End Tutorial 8

## tutorial 9 ##


## Tutorial 9 ##
#Interplanetary('mars') # good
#Interplanetary('jupiter') # good
#Interplanetary('venus') # good
#Interplanetary('moon') # good enough
#Interplanetary('saturn',theta=180) # GREAT
#Interplanetary('venus',parkPlanet = 6600, parkTarget = 6400) # First fucking try.
#Interplanetary('mars',parkPlanet = 6697, parkTarget = 3580, theta = 50) # wait time off by ~  20 hours
#Interplanetary('saturn',parkPlanet = 180, parkTarget = 1000) # off by 0.01
#Interplanetary('jupiter',parkPlanet = 200, parkTarget = 1000) # good
#Interplanetary('mercury', parkPlanet = 180, parkTarget = 200) # slightly off. maybe because need to consider actual ellipse

## Past Exams ##

# S1 2017 White:

#Q1: (all right except maybe last deltaV from transfer)
#orbit('?',4200,T=14*(24*3600)) #jupiter
#Interplanetary('jupiter')
#orbit('?','?',T=53.5*24*3600) #jupiter - a = 4092933
#Transfer(3273300-r0,4200,4092933-r0,4092933-r0,0,0)

#Q2:
#inclination(925,0,-360/90) #. a)L: i = 49.85
#inclination(925,0,-360/135) # i = 64.54

#Transfer(925,925,925,925,inclination(925,0,-360/90),inclination(925,0,-360/135)) # DeltaV = 1.889.
# ^ only need this one line of code. Amazing.

#Q3:



print('\n')
print('Angular Velocity of '+planet.upper()+': ',ssig(omega0), 'rad/sec')
print('Semi Major Axis of '+planet.upper()+': ',ssig(a0), 'km')
print('MU of '+planet.upper()+': ',ssig(mu), 'm^3 / s^2')
print('MU of SUN: ', ssig(muPlanet['sun']), 'm^3 / s^2')

#assignment
#orbit(70000+r0,3000+r0) #moon - gateway






# RETIRED FUNCTIONS:
def transfer(h1a,h1p,h2a,h2p,theta=None,launch='perigee'):

    r1a = h1a + r0
    r1p = h1p + r0
    r2a = h2a + r0
    r2p = h2p + r0

    a1 = (r1a+r1p)/2
    a2 = (r2a+r2p)/2

    epsilon1 = -mu / (2 * a1)
    epsilon2 = -mu / (2 * a2)

    if launch == 'apogee':
        a_trans = (r1a + r2p)/2
        epsilon_trans = -mu / (2 * a_trans)
        v_trans1 = velocity(r1a,epsilon_trans)
        v1 = velocity(r1a, epsilon1)
    else:
        a_trans=(r1p+r2p)/2
        epsilon_trans = -mu / (2 * a_trans)
        v_trans1 = velocity(r1p, epsilon_trans)
        v1 = velocity(r1p, epsilon1)




    v2 = velocity(a2, epsilon2)

    v_trans2 = velocity(r2p, epsilon_trans)

    deltaV1 = abs(v_trans1 - v1)
    deltaV2 = abs(v_trans2 - v2)

    deltaV = deltaV1 + deltaV2

    TOF = np.pi*np.sqrt(a_trans**3/mu)

    omega1 = np.sqrt(mu / a1**3)
    omega2 = np.sqrt(mu / a2**3)

    alpha = omega2*TOF #lead angle

    thetaF = np.pi - alpha #final angle

    # Rendezvous component
    if isinstance(theta,float) or isinstance(theta,int):
        thetaRad = theta * np.pi / 180 #rad conversion

        alpha = omega2 * TOF  # lead angle
        alphaDeg = ssig(alpha * 180 / np.pi) # degree conversion

        thetaF = np.pi - alpha  # final angle
        thetaF_deg = alpha * 180 / np.pi

        waitTime = (thetaF - thetaRad) / (omega2 - omega1)

        alphaDeg,thetaF_deg,waitTime = ssig(alphaDeg), ssig(thetaF_deg), ssig(waitTime)

        print('\n')
        print('----RENDEZVOUS----')
        print('Initial Angle: ',theta, 'deg', ' / ',thetaRad, 'rad')
        print('Final Angle: ',thetaF_deg, 'deg'' / ',thetaF, 'rad')
        print('Lead Angle: ',alphaDeg, 'deg'' / ',alpha, 'rad')
        print('Wait Time: ('+str(waitTime)+') '+tconv(waitTime))

    print('\n')

    omega1,omega2,deltaV,deltaV1,deltaV2,v1,v2 = ssig([omega1,omega2,deltaV,deltaV1,deltaV2,v1,v2])

    print('DeltaV: ', deltaV, ' km/s')
    print('TOF = ('+str(TOF)+') '+tconv(TOF))
    print('DeltaV1: ', deltaV1, ' km/s')
    print('DeltaV2: ', deltaV2, ' km/s')
    print('v1: ',v1, 'km/s')
    print('v2: ',v2, 'km/s')
    print('v1 Transfer: ', v_trans1)
    print('v2 Transfer: ', v_trans2)
    print('Angular Velocity of Low Orbit: ', omega1, 'rad/sec')
    print('Angular Velocity of High Orbit: ', omega2, 'rad/sec')

    return v_trans1,v_trans2