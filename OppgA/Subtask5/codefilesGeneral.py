import numpy as np
import xlrd
from prettytable import PrettyTable
"""
python program for solving a loadflow problem with the newton raphson method

created by: Erlend Løvstakken, Einar Ingmar Frøhaug, Eivind Roland, Håkon Broch and Eilert Henriksen

"""


def createlines(linedata): #creates a dictonary for all the lines in the system, and saves their address
    linedict={}
    for j in range(1,len(linedata.col_values(0))):
        linedict["line"+str(j)]=line(linedata,j)
    return linedict




def createbus(busdata,linedict): #creates a dictonary for all the buses in the system, and saves their adress
    busdict={}
    for j in range(1,len(busdata.col_values(0))):
        busdict['bus'+str(j)]=bus(busdata,j)
    for q in linedict:
        busdict['bus'+str(linedict[q].fromval)].lines.append(linedict[q])
    return busdict



class bus:
    def __init__(self,exceldata,i):
        self.number=int(exceldata.cell_value(i,0)) #integer for busnr
        self.pspec=exceldata.cell_value(i,1) #float
        self.qspec=exceldata.cell_value(i,2) #float
        self.vspec=exceldata.cell_value(i,3) #float
        self.qmin=exceldata.cell_value(i,4) #float
        self.qmax=exceldata.cell_value(i,5) #float
        self.d=np.nan_to_num(exceldata.cell_value(i,6)) #float
        self.cap=exceldata.cell_value(i,7)*1j #float
        self.pmax=exceldata.cell_value(i,8) # float, not necessary
        self.lines=[]
        self.type=None
        self.exceedlim=False

        


class line:
    def __init__(self,exceldata,i):
        self.fromval=int(exceldata.cell_value(i,0)) #int
        self.to=int(exceldata.cell_value(i,1)) #int
        self.R=exceldata.cell_value(i,2) #float
        self.X=exceldata.cell_value(i,3)*1j #complex
        self.shunt=exceldata.cell_value(i,4)*1j #complex
        self.y=1/(self.R+self.X) #complex


def ybus(linedict,busdict):  #need busdict to get the amount of buses, in a general case the linedict wouldnt always give the right number of buses.
    uptri=np.round((np.zeros((int(len(busdict)),int(len(busdict))),dtype=complex)),2) #creates an empty array
    shunts=np.nan_to_num(np.zeros(len(busdict),dtype=complex)) #nan to numpy is implemented to deal with an issue were some values are defined as NaN.
    for k in linedict:
        uptri[linedict[k].fromval-1,linedict[k].to-1]+=linedict[k].y #creates the upper triangular of the ybus
        shunts[linedict[k].fromval-1]+=linedict[k].shunt/2 #half line model
        shunts[linedict[k].to-1]+=linedict[k].shunt/2  #add all shunts to a list.
    lotri=np.transpose(uptri) #transposes upper triangular, which gives same values, but in the lower triangular
    ybus=-np.add(uptri,lotri) #sum these two and make them negative
    diag=-(ybus.sum(1)) #Sum up the rows, and change sign
    for i in range(len(diag)): #iterate through the diagonal of the ybus
        if busdict['bus'+str(i+1)].cap!=0:
            ybus[i,i]=diag[i]+shunts[i]+1/busdict['bus'+str(i+1)].cap
        else:
            ybus[i,i]=diag[i]+shunts[i] #implement the diagonal+shunts into the ybus
        
    return ybus

def typecheck(busdict): #use the assumption that slack bus is specified in excel file with vspec=1 and d=0
    for k in busdict:
        if bool(busdict[k].pspec) and bool(busdict[k].vspec):
            busdict[k].type='PV'
        else:
            busdict[k].type='PQ'
        if busdict[k].d==0 and busdict[k].vspec==1 and busdict[k].pspec==0:
            busdict[k].type='slack'
    return busdict



def powercalc(Y, buses):
    n_PV = 0      #Number of PV buses.
    n_PQ = 0
    i = 0

    k = 0
    N = np.size(buses) #Number of buses
    indexvectorQ = np.zeros(0,dtype=int) #Vector containing PQ bus indexes
    indexvectorP = np.zeros(0,dtype=int) #Vector containing PV bus indexes
    Vcomp = np.zeros(N,dtype=np.float64) #Voltage matrix with V[0] = VBus1 and so on.
    D = np.zeros(N,dtype=np.float64)      #Angle matrix with d[0] = dBus1 and so on.

    for b in buses:
        Vcomp[k] = b.vspec
        D[k] = b.d
        if b.type == "PV":
            indexvectorP = np.append(indexvectorP, b.number - 1)
            n_PV += 1

        elif b.type == "PQ":
            indexvectorP = np.append(indexvectorP,b.number - 1)
            indexvectorQ = np.append(indexvectorQ,b.number - 1)#Subtracting to account for bus 1 having index 0
            n_PQ += 1

        k += 1
    n = n_PV+2*n_PQ #J_dimension4
    Yabs = abs(Y)
    Yang = np.angle(Y)

    V = np.abs(Vcomp)
    calcpower = np.zeros(n,dtype=np.float64)

    for i in range(n):
        #first calculating powers
        if i <= n_PV+n_PQ-1:
            index_power = indexvectorP[i]
            for k in range(N): #Calculating real power
                calcpower[i] += V[k]*V[index_power]*Yabs[index_power][k]*np.cos(D[index_power]-D[k]-Yang[index_power][k])
        elif i > n_PV+n_PQ-1:
            index_power = indexvectorQ[i-n_PV-n_PQ]
            for k in range(N): #Calculating Reactive Power
                calcpower[i] += V[k]*V[index_power]*Yabs[index_power][k]*np.sin(D[index_power]-D[k]-Yang[index_power][k])
    return calcpower

def jacobianpower(Y,buses):


    n_PV = 0      #Number of PV buses.
    n_PQ = 0
    i = 0
    j = 0
    k = 0
    N = np.size(buses) #Number of buses
    indexvectorQ = np.zeros(0,dtype=int) #Vector containing PQ bus indexes
    indexvectorP = np.zeros(0,dtype=int) #Vector containing PV bus indexes
    Vcomp = np.zeros(N,dtype=np.float64) #Voltage matrix with V[0] = VBus1 and so on.
    D = np.zeros(N,dtype=np.float64)      #Angle matrix with d[0] = dBus1 and so on.

    for b in buses:
        Vcomp[k] = b.vspec
        D[k] = b.d
        if b.type == "PV":
            indexvectorP = np.append(indexvectorP, b.number - 1)
            n_PV += 1

        elif b.type == "PQ":
            indexvectorP = np.append(indexvectorP,b.number - 1)
            indexvectorQ = np.append(indexvectorQ,b.number - 1)#Subtracting to account for bus 1 having index 0
            n_PQ += 1

        k += 1

    n = n_PV+2*n_PQ #J_dimension
    Yabs = abs(Y)
    Yang = np.angle(Y)

    V = np.abs(Vcomp)
    J = np.zeros((n,n),dtype=np.float64)



    for i in range(n):
        #now calculating jacobian
        for j in range(n):

            if i <= n_PV+n_PQ-1:


                if j <=n_PV+n_PQ-1: #J1

                    index_i = indexvectorP[i] #From first P and onwards
                    index_j = indexvectorP[j] #From delta2 and onwards
                    if i is not j:
                        J[i][j] =   V[index_i]*V[index_j]*Yabs[index_i][index_j]*np.sin(D[index_i]-D[index_j]-Yang[index_i][index_j])
                    else:
                        J[i][j] = V[index_j]**2*Yabs[index_j][index_j]*np.sin(-Yang[index_j][index_j]) #Pre-removing where the derivative is zero
                        for k in range(N):
                            J[i][j]-= V[index_j]*V[k]*Yabs[index_j][k]*np.sin(D[index_i]-D[k]-Yang[index_i][k])

                elif j > n_PV+n_PQ-1: #J2
                    index_i = indexvectorP[i]      #From P2 and onwards
                    index_j = indexvectorQ[j-n_PV-n_PQ] #From V_nPV+1 and onwards
                    if i is not j:
                        J[i][j] =   V[index_i]*Yabs[index_i][index_j]*np.cos(D[index_i]-D[index_j]-Yang[index_i][index_j])
                    else:
                        J[i][j] = V[index_j]*Yabs[index_j][index_j]*np.cos(-Yang[index_j][index_j]) #Pre adding because the derivative when k = j is different
                        for k in range(N):
                            J[i][j] += V[k]*Yabs[index_i][k]*np.cos(D[index_i]-D[k]-Yang[index_i][k])

            elif i > n_PV+n_PQ-1:


                 if j <=n_PV+n_PQ-1: #J3
                    index_i = indexvectorQ[i-n_PV-n_PQ]#From Q_nPV+1 and onwards
                    index_j = indexvectorP[j]         #From d2 and onwards
                    if i is not j:
                        J[i][j] = - V[index_i]*V[index_j]*Yabs[index_i][index_j]*np.cos(D[index_i]-D[index_j]-Yang[index_i][index_j])
                    else:
                      J[i][j] = -V[index_j]**2*Yabs[index_j][index_j]*np.cos(-Yang[index_j][index_j]) #Pre-removing where the derivative is zero
                      for k in range(N):
                          J[i][j]+=V[index_i]*V[k]*Yabs[index_i][k]*np.cos(D[index_i]-D[k]-Yang[index_i][k])

                 elif j > n_PV+n_PQ-1: #J4
                    index_i = indexvectorQ[i-n_PV-n_PQ] #From Q_nPV+1 and onwards
                    index_j = indexvectorQ[j-n_PV-n_PQ] #From V_nPV+1 and onwards
                    if i is not j:
                        J[i][j] =   V[index_i]*Yabs[index_i][index_j]*np.sin(D[index_i]-D[index_j]-Yang[index_i][index_j])
                    else:
                        J[i][j] = V[index_j]*Yabs[index_j][index_j]*np.sin(-Yang[index_j][index_j]) #Pre adding because the derivative of k = j is different
                        for k in range(N):
                            J[i][j]+=V[k]*Yabs[index_i][k]*np.sin(D[index_i]-D[k]-Yang[index_i][k])

    return J



def reactivepowercalc(Y, buses):
    

    k = 0
    N = np.size(buses) #Number of buses 
    indexvectorP = np.zeros(0,dtype=int) #Vector containing PV bus indexes
    Vcomp = np.zeros(N,dtype=np.float64) #Voltage matrix with V[0] = VBus1 and so on.
    D = np.zeros(N,dtype=np.float64)      #Angle matrix with d[0] = dBus1 and so on.

    for b in buses:
        Vcomp[k] = b.vspec
        D[k] = b.d
        k+=1
        if b.type == "PV":
            indexvectorP = np.append(indexvectorP, b.number - 1)
                
    Yabs = abs(Y)
    Yang = np.angle(Y)
    
    V = np.abs(Vcomp)    

    
    for i in indexvectorP: 
        #first calculating powers
        buses[i].qspec =0.0
        for k in range(N): #Calculating Reactive Power
            buses[i].qspec += V[k]*V[i]*Yabs[i][k]*np.sin(D[i]-D[k]-Yang[i][k])
    return buses 


def ActualPowers(buses):

    powers = np.zeros(0,dtype=np.float64)
    indexesP = np.zeros(0,dtype=np.int)
    indexesQ = np.zeros(0,dtype=np.int)
    qp = np.zeros(0,dtype=np.float64)
    
    for b in buses:#Vector of powers corresponding to mismatch vector
        if b.type == "PV":
            powers = np.append(powers,b.pspec)
            indexesP = np.append(indexesP,b.number-1)
        elif b.type == "PQ":
            powers = np.append(powers,b.pspec)
            indexesP = np.append(indexesP,b.number-1)
            qp = np.append(qp,b.qspec)
            indexesQ = np.append(indexesQ,b.number-1)
    powers = np.append(powers,qp) #appending reactive powers        
    return powers, indexesP,indexesQ

def slackpower(ybus,buses):
    V = np.zeros(len(buses),dtype=np.float64) #Voltage matrix with V[0] = VBus1 and so on.
    D = np.zeros(len(buses),dtype=np.float64)      #Angle matrix with d[0] = dBus1 and so on.
    Yabs = abs(ybus)
    Yang = np.angle(ybus)
    
    k = 0
    for b in buses:
        V[k] = abs(b.vspec)
        D[k] = b.d
        k = k+1
        
    for b in buses:
        if b.type == 'slack':
            b.pspec = 0
            b.qspec = 0
            index_slack = b.number -1 
            for i in range(len(buses)):
                b.pspec += V[i]*V[index_slack]*Yabs[index_slack][i]*np.cos(D[index_slack]-D[i]-Yang[index_slack][i])
                b.qspec += V[i]*V[index_slack]*Yabs[index_slack][i]*np.sin(D[index_slack]-D[i]-Yang[index_slack][i])
    return buses                