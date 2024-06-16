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
        self.alphaE = - 1  / 3

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

        

if __name__=="__main__":
    
    scrn = Screen_Delta()

    # scrn.Analytic_delta_screening()

    # scrn.Numeric_delta_screening(1 - np.array([1e-5 * i for i in range(2000)])) # reproduces Ina's code
    
