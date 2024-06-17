import numpy as np
import scipy as scp
import matplotlib.pyplot as plt



class Screen_Delta:

    def __init__(self):

        self.init_constants()


    def init_constants(self):

        self.x = "hello world"

        # self.beta0 = 
        self.alphaHe = 1 / 3
        self.alphaE = - 1 / 3

        self.Rearth = 3.2322 * 10**13 # [eV^-1]
        # self.rhoE = self.beta0 / (5 * (2 self.alphaHe))
        self.rhoE = 0.00971433
        self.rhoEarth = 0.0485717

        self.A1d = 1


    def phShIntnewA(self,rMax,r,A1d=1):

        return np.real((2 * r**2 - rMax**2) - np.emath.sqrt((2 * r**2 - rMax**2)**2 + 32 * (r**4 - rMax**2 * r**2))) / (16 * r**2 * A1d)


    def external_phi_to_total(self,func):

        return lambda r: self.rhoE * np.heaviside(self.Rearth * (1-r),0) + func(self.Rearth * r) * np.heaviside((r-1) * self.Rearth,0) # [eV^-1]
    

    def Analytic_delta_screening(self):

        for A1d in [1,2,3,5,6,7]:

            rMax = scp.optimize.fsolve(lambda rMax, r : self.phShIntnewA(rMax,r,A1d) - self.rhoE,self.Rearth,(self.Rearth),maxfev=300)[0]
            # print(rMax)

            f = self.external_phi_to_total(lambda r: self.phShIntnewA(rMax,r,A1d))

            rs = np.linspace(0.995,1.02,200)

            fig, ax = plt.subplots()
            ax.plot(rs,rs**2 * f(rs))
            plt.ylim(0,0.01)
            plt.show()
            

    def smallrETerm(self,RMax,r,phi,phiRMax):

        f = lambda r: np.real(np.sqrt(1 + phi) * np.emath.sqrt(1 + phi - RMax**2 / r**2 * (1 + phiRMax))) * np.heaviside(RMax - r,1)

        return f(r)
 

    def Numeric_delta_screening(self,rs):

        self.rhoConstVList = []

        i = 0

        while True:

            if i%100 == 0:
                print(i)
            
            fphi = lambda phi: self.alphaHe * (1 - 2 * phi) + self.alphaE * ( 1 + phi - self.smallrETerm(1,rs[i],phi,0))

            try:
                
                root = scp.optimize.brentq(fphi,10e-5,self.rhoEarth)

            except ValueError: # brentq requires function to evaluate to different signed values on the boundary
                
                i+=1 

                continue

            self.rhoConstVList.append([rs[i],root])

            if root >= self.rhoEarth - 10e-4:

                break

            i += 1
                
        print("done {}:".format(i),self.rhoConstVList[-1])

        self.rhoConstVList = np.array(self.rhoConstVList)

        fig, ax = plt.subplots()
        ax.plot(np.linspace(10e-5,self.rhoEarth),fphi(np.linspace(10e-5,self.rhoEarth)))
        plt.show()

        fig, ax = plt.subplots()
        ax.plot(self.rhoConstVList[:,0],self.rhoConstVList[:,1])
        plt.show()

        

class Screen_MB:

    def __init__(self):

        self.init_constants()

        self.get_MB_factor()


    def init_constants(self):

        self.vInfList = 0.05 * 1.2**np.arange(25) # [] v_inf / v_0
        self.nvs = np.size(self.vInfList)

        self.phiList = [[1,0]]
        self.RMaxList = [[1,0]]

        self.minPoint = 0.9
        # self.pointNumber = int(1e4)
        self.pointNumber = int(1e2)
        self.scanPoints = [1 - i * (1- self.minPoint) / self.pointNumber for i in range(self.pointNumber)]

        self.phiE =  0.0485717


    def get_MB_factor(self):

        """
        i - [] element of self.vInfList to evaluate MB factor on 
        """

        def MB_arg(j):

            return self.vInfList[j]**2 * np.exp(-3 / 2 * self.vInfList[j]**2)

        # if not 0 <= i < self.nvs:
            
        #     print("Index passed is out of range. Must be i < {}".format(self.nvs))

        #     return

        # MBNorm = sum([self.vInfList[j]**2 * np.exp(-self.vInfList[j]**2) for j in range(self.nvs)])
        MBNorm = sum([MB_arg(j) for j in range(self.nvs)])
       
        self.MB_factor = np.array([MB_arg(j) for j in range(self.nvs)] / MBNorm)


    def complexF(self,i,RMax,r):

        ##### WHERE IS THIS DEFINED?

        return RMax**2 / r**2 * (self.vInfList[i]**2 + self.RMaxList[i][1]**2) - self.vInfList[i]**2


    def rootF(self,i,RMax,r,phi,n):

        return (-4 * phi + 4 * phi**2 - 14 / 3 * phi**3) / n + \
            self.MB_factor[i] / self.vInfList[i] * np.real(np.emath.sqrt(phi - self.complexF(i,RMax,r)))


    def derivF(self,i,RMax,r,phi,n):

        np.real(((-3 + 9 * phi) - \
                 np.sum(self.MB_factor / (2 * np.emath.sqrt(self.vInfList**2+phi) * self.vInfList))) / n \
                    + (self.MB_factor[i] / (2 * np.emath.sqrt(phi - self.complexF(i,RMax,r)))))


    def Get_phi_MB(self):

        earthFlag = False
        zoomFlag = False

        for j in range(self.nvs):
            
            phiListWorking = []
            delNum = 0
            zoomCount = 0

            if earthFlag or zoomFlag:
                break

            print(j)

            for k in range(self.pointNumber):

                if zoomCount > 12:

                    zoomFlag = True

                    print("You've zoomed in more than 12 times, problem?")

                    break

                print(self.RMaxList)
                print(self.scanPoints)
                # print([self.RMaxList[l][0],self.scanPoints[k]])
                print([self.complexF(l,self.RMaxList[l][0],self.scanPoints[k]) for l in range(j+1)])

                phiMin = 1e-15 + np.amax([self.complexF(l,self.RMaxList[l][0],self.scanPoints[k]) for l in range(j+1)])

                for l in range(j):

                    if np.emath.sqrt(self.vInfList[l]**2 + phiMin - \
                                     self.RMaxList[l][0]**2 / self.scanPoints[k]**2 \
                                        * (self.vInfList[l]**2 + self.RMaxList[l][1])) <= 1e-8:
                        
                        phiMin += 1e-12

                # get phiMax
                try:
                
                    derivFsum = lambda phi: np.sum([self.derivF(i,self.RMaxList[i][0],self.scanPoints[k],phi,j) for i in range(j)])

                    phiMax = scp.optimize.brentq(derivFsum,phiMin,0.1)

                except ValueError: # brentq requires function to evaluate to different signed values on the boundary
                    
                    k+=1 

                    continue

                # get res
                try:
                
                    rootFsum = lambda phi: np.sum([self.rootF(i,self.RMaxList[i][0],self.scanPoints[k],phi,j) for i in range(j)])

                    res = scp.optimize.brentq(rootFsum,phiMin - 1e-12,phiMax)

                    Failed = False

                except ValueError: # brentq requires function to evaluate to different signed values on the boundary
                    
                    Failed = True

                if not Failed:

                    phiNew = res

                    if phiNew > 2 / 3 * self.phiE:

                        earthFlag = True

                        self.phiList.extend(phiListWorking)

                        print("Reached Earth phi")

                        break

                    phiListWorking.append([self.scanPoints[k],phiNew])

                    print("philist is")
                    print(self.phiList)
                    print(phiListWorking)

                    self.phiList.extend(phiListWorking)
                    phiListarray = np.array(self.phiList)

                    print(phiListarray)

                    phifunc = scp.interpolate.CubicSpline(phiListarray[:][0],phiListarray[:][1])


if __name__=="__main__":
    
    scrn = Screen_Delta()

    # scrn.Analytic_delta_screening()

    # scrn.Numeric_delta_screening(1 - np.array([1e-5 * i for i in range(2000)])) # reproduces Ina's code
    
    scrn_MB = Screen_MB()

    # print(scrn_MB.vInfList[3])

    # print(scrn_MB.MB_factor)

    # print(scrn_MB.rootF(0,scrn_MB.RMaxList[0][0],scrn_MB.scanPoints[1],0,2))

    scrn_MB.Get_phi_MB()