import numpy as np

"""
Assumptions:
Bus type doesnt change
Only this 3 bus system
with random variables for Rs,Xs,P2,P3,Q3 and V2

"""
def P(v,y,d,i):
    P=0
    for j in range(3):
        P+=abs(y[i,j])*abs(v[j])*abs(v[i])*np.cos(d[i]- d[j]-np.angle(y[i,j]))
    return P

def Q(v,y,d,i):
    Q=0
    for j in range(3):
        Q+=abs(y[i,j])*abs(v[j])*abs(v[i])*np.sin(d[i]- d[j]-np.angle(y[i,j]))
    return Q

def my_threebus_NR_10101(r12, x12, r13, x13, p2, v2, p3, q3):
    myangles = np.array([0, 0, 0.0]) #changed to numpy array, needs to be float
    myVs = np.array([1, v2, 1]) #flat start, changed v2 to given v2
    uspec=np.array([[p2],[p3],[q3]])
    y12=1/(r12+x12*1j)
    y13=1/(r13+ x13*1j)
    y=np.array([[y12+y13,-y12,-y13],[-y12,y12,0],[-y13,0,y13]]) #since this is a simple, set problem, there is no point in making a more advanced solution to Ybus
    counter=0 #just to check for how many iterations are needed to finish and hinder overflow
    testlimit=20 #max amount of legal iterations before overflow, decided by me
    u=np.array([[P(myVs,y,myangles,1)],[P(myVs,y, myangles,2)],[Q(myVs,y,myangles,2)]])
    x=np.array([[myangles[1]],[myangles[2]],[myVs[2]]])
    eps=0.000001 #legal fault margin

    while abs(uspec[0]-u[0])>eps and abs(uspec[1]-u[1])>eps and abs(uspec[2]-u[2])>eps and counter<testlimit : #last condition is to hinder overflow.
        J11=-abs(y[1,0])*abs(myVs[1])* abs(myVs[0])*np.sin(myangles[1]-myangles[0]-np.angle(y[1,0]))-abs(y[1,2])*abs(myVs[1])* abs(myVs[2])*np.sin(myangles[1]-myangles[2]-np.angle(y[1,2]))
        J12=-abs(y[1,2])*abs(myVs[1])* abs(myVs[2])*np.sin(myangles[1]-myangles[2]-np.angle(y[1,2]))
        J13=abs(y[1,2])*abs(myVs[1])*np.cos(myangles[1]-myangles[2]-np.angle(y[1,2]))
        J21=-abs(y[2,1])*abs(myVs[2])* abs(myVs[1])*np.sin(myangles[2]-myangles[1]-np.angle(y[2,1])) #is 0
        J22=-abs(y[2,0])*abs(myVs[2])* abs(myVs[0])*np.sin(myangles[2]-myangles[0]-np.angle(y[2,0]))-abs(y[2,1])*abs(myVs[2])* abs(myVs[1])*np.sin(myangles[2]-myangles[1]-np.angle(y[2,1]))
        J23=abs(y[2,0])*abs(myVs[0])*np.cos(myangles[2]-myangles[0]-np.angle(y[2,0]))+abs(y[2,1])*abs(myVs[1])*np.cos(myangles[2]-myangles[1]-np.angle(y[2,1]))+2*abs(y[2,2])*abs(myVs[2])*np.cos(-np.angle(y[2,2]))
        J31=abs(y[2,1])*abs(myVs[2])* abs(myVs[1])*np.cos(myangles[2]-myangles[1]-np.angle(y[2,1])) #is 0
        J32=abs(y[2,0])*abs(myVs[2])* abs(myVs[0])*np.cos(myangles[2]-myangles[0]-np.angle(y[2,0]))+abs((y[2,1]))*abs(myVs[2])* abs(myVs[1])*np.cos(myangles[2]-myangles[1]-np.angle(y[2,1]))
        J33=abs(y[2,0])*abs(myVs[0])*np.sin(myangles[2]-myangles[0]-np.angle(y[2,0]))+abs(y[2,1])*abs(myVs[1])*np.sin(myangles[2]-myangles[1]-np.angle(y[2,1]))+2*abs(y[2,2])*abs(myVs[2])*np.sin(-np.angle(y[2,2]))
        #Not the prettiest solution, but probably most time efficent in this task... since there is quite the clear pattern
        J=np.array([[J11,J12,J13],[J21,J22,J23],[J31,J32, J33]])
        #dx=np.linalg.inv(J)@(uspec-u)
        dx=np.linalg.solve(J,uspec-u)
        x=x+dx
        myangles[1]= x[0]
        myangles[2]= x[1]
        myVs[2]=x[2]
        u=np.array([[P(myVs,y,myangles,1)],[P(myVs,y, myangles,2)],[Q(myVs,y,myangles,2)]])

        counter+=1
    print("Finished in " , counter , " iterations")
    #no need to do any pretty visuals for this problem since it is just to check if it is correct.
    return myangles, myVs




if __name__ == '__main__':
    import random
    from newton_raphson import check_3bus

    np.set_printoptions(precision=3)
    np.set_printoptions(suppress=True)

    print('Testing 3-bus NR loadflow solver for candidate: 10101')

    for i in range(0, 5):
        # All values in pu
        R12 = random.randint(3, 6) / 100  # R between 0.03 and 0.06
        X12 = random.randint(20, 30) / 100  # X between 0.2 and 0.3

        R13 = random.randint(2, 4) / 100  # R between 0.02 and 0.04
        X13 = random.randint(10, 20) / 100  # X between 0.1 and 0.2

        P2 = random.randint(1, 10) / 10  # P between 0.1 and 1
        V2 = 1 + random.randint(-5, 5) / 100  # V between 0.95 and 1.05

        P3 = -random.randint(1, 10) / 10  # P between -0.1 and -1
        Q3 = -random.randint(1, 10) / 10  # Q between -0.1 and -1

        myangles, myVs = my_threebus_NR_10101(R12, X12, R13, X13, P2, V2, P3, Q3)

        print('Test {}: '.format(i + 1), end='')
        check_3bus(R12, X12, R13, X13, P2, V2, P3, Q3,myangles, myVs,set_print= False)
