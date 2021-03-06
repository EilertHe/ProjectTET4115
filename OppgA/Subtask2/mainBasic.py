from codefilesBasic import *
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
count=0
for b in buses:
    if b.type=='PQ':
        b.vspec=1.0    
        
while not conv:
    count+=1
    powers,indexesP,indexesQ=ActualPowers(buses)
    pwrcalc=powercalc(yBus,buses)
    k = 0
    J= jacobianpower(yBus, buses)
    u=powers-pwrcalc
    deltax=np.linalg.solve(J,u)

    for i in indexesP:
        buses[i].d += deltax[k]
        k += 1

    for i in indexesQ:
        buses[i].vspec += deltax[k]
        k += 1

    #Finds remaining reactive powers
    buses=reactivepowercalc(yBus,buses)
    
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



