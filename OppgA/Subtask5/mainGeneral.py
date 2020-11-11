from codefilesGeneral import *

#initializes the files from excel
filename=(r"C:\example.xlsx") #Remember to change this to your own path for example.xlsx, it is important to include the first r outside ""
file = xlrd.open_workbook(filename)
busdata=file.sheet_by_index(0)
linedata=file.sheet_by_index(1)
testlines=createlines(linedata)
testbuses=createbus(busdata,testlines)
yBus=ybus(testlines,testbuses)
typecheck(testbuses)
conv = False


buses = np.array(list(testbuses.values())) #change from dict to array...
Vprev=np.zeros(len(buses),dtype=np.float64)
Vorg=Vprev
count=0 #count the iterations needed to complete assignment
for b in buses:
    if b.type=='PQ':
        b.vspec=1.0
    if b.type=='PV':
        Vorg[b.number]=b.vspec      
        
while not conv:
    count+=1
    powers,indexesP,indexesQ=ActualPowers(buses)
    pwrcalc=powercalc(yBus,buses)
    k = 0
    J= jacobianpower(yBus, buses)
    u=powers-pwrcalc #mismatch vector
    deltax=np.linalg.solve(J,u)

    for i in indexesP:
        buses[i].d += deltax[k]
        k += 1

    for i in indexesQ:
        buses[i].vspec += deltax[k]
        k += 1

    #finds remaining powers for limit check
    buses=reactivepowercalc(yBus,buses)
    
    #THIS CHECKS FOR LIMITS ON Qs
    for index,b in enumerate(buses): #check for bus changes
        if b.qspec >= b.qmax and b.type=='PV':
            b.qspec=b.qmax
            b.type='PQ'
            b.exceedlim=True
            Vprev[index]=b.vspec
        elif b.qspec<b.qmin and b.type=='PV':
            b.qspec=b.qmin
            b.type='PQ'
            b.exceedlim=True
            Vprev[index]=b.vspec

#check if changed PQ bus can return to PV bus
    for ind in range(len(Vprev)): 
        if buses[ind].exceedlim and  (buses[ind].qspec==buses[ind].qmax and buses[ind].vspec>Vorg[ind]) or (buses[ind].qspec==buses[ind].qmin and buses[ind].vspec<Vorg[ind]): #legge til under lovlig verdi greie
            buses[ind].type='PV'
            buses[ind].exceedlim=False
            Vprev[ind]=0


    if np.sum(abs(deltax)) < 10**-5:

        conv = True
buses=slackpower(yBus,buses) #calculate the power needed from slack
        
 #This is only visualizations, not important for code.
 
print("found a solution in: ",count,"iterations.")

for b in buses:
    if b.exceedlim:
        print("This system is not realizable")
        break

vectorV=[]
for j in range(len(yBus)):
    vectorV.append(str(j+1))    
ybusvisual=PrettyTable()
ybusvisual.add_column('YBUS',vectorV)
for i in range(len(yBus)):
    ybusvisual.add_column(str(i+1),np.round(yBus[i,:],4))
print(ybusvisual)
Finalspecs=PrettyTable(["Bus","voltage[pu]","angle[rad]","active power[pu]","reactive power[pu]"])
for b in buses:
    Finalspecs.add_row([np.round(b.number,5),np.round(b.vspec,5),np.round(b.d,5),np.round(b.pspec,5),np.round(b.qspec,5)])
    if b.exceedlim:
        print('Bus',b.number,'exceeds legal value for voltage')
print(Finalspecs) 



