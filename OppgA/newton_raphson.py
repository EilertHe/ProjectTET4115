#
# Power system analysis project: Fully general Newton-Raphson, used as benchmark for 3bus system.
#
# Developed by: Sjur Føyen, Dept. of Electric Power Engineering, NTNU
#
#       September 2020
#

import numpy as np
from typing import List


class Bus:
    def __init__(self,Busdata):
        self.num = int(Busdata[0])
        self.type = Busdata[1]
        self.P_spec = Busdata[2]
        self.Q_spec = Busdata[3]
        self.V_spec = Busdata[4]
        self.Qlim = [Busdata[5], Busdata[6]]
        self.Plim = [Busdata[7], Busdata[8]]
        if not np.isnan(Busdata[9]):
            self.comp = Busdata[9]
        else:
            self.comp = 0

        self.connected_lines = {}  # type: Dict[Bus,Line]
        self.powerflow = {}  # type: Dict[Bus,float]
        # Load flow internal and result variables
        self.n = 0
        self.angle = 0
        if not np.isnan(self.V_spec):
            self.V = self.V_spec
        else:
            self.V = 1
        self.P = 0
        self.Q = 0
        self.j_ind1 = None
        self.j_ind2 = None


    def print(self):

        stri = "{:^6}|{:^6}|{:^8}|{:^8}|{:^8}|{:^8}|".format(self.num,self.type,self.P_spec,self.Q_spec,self.V_spec,
                                                            np.round(self.angle,4))
        stri += "{:^8}|{:^8}|{:^8}".format(np.round(self.V,4),np.round(self.P,4),np.round(self.Q,4))
        print(stri)

    def printpf(self):
        for bus, pf in self.powerflow.items():
            print('From {} to {}:  {}'.format(self.num, bus.num, np.round(pf, 4)))

class Line:
    def __init__(self,Linedata):
        self.nodes = Linedata[0:2]
        self.R = Linedata[2]
        self.X = Linedata[3]
        self.y = 1/(self.R+1j*self.X)


    def other(self,bus: Bus):
        if bus.num != self.nodes[1]:
            return self.nodes[1]
        else:
            return self.nodes[0]


def admittance_matrix(bd):
    """
    Function is for the moment deprecated, as
    the NR solver adopts the line-oriented approach.
    :param bd:
    :return:
    """
    Nbus = len(bd)
    Ybus = 1j*np.zeros((Nbus,Nbus))
    for i in bd:
        selfY = i.comp
        for bus,line in i.connected_lines.items():
            Ybus[i.num,bus.num] = -line.y
            selfY += line.y
        Ybus[i.num,i.num] = selfY # set self-admittance
    #print(Ybus)
    return Ybus


def printsys(bd:List[Bus]):
    str = "{:^6}|{:^6}|{:^8}|{:^8}|{:^8}|{:^8}|{:^8}|{:^8}|{:^8}".format("#","Type","Pspec","Qspec",
                                                                         "Vspec","theta","V","P","Q")
    print(str)

    for i in bd:
        i.print()


def pv_limits(bd:List[Bus]):
    """
    #TODO implement as class method
    Checks PV bus limits on reactive power
    :param bd: list of Bus objects
    :return: none
    """
    for i in bd:
        if i.type == 'PV':
            if i.Q > i.Qlim[1]:
                print('---------------------------------------------------------------------------')
                print('!!!!!!!!!!!!!!!!! Upper Q limit violated at bus {} !!!!!!!!!!!!!!!!!'.format(i.num))
                i.type = 'PQ'
                i.Q_spec = i.Qlim[1]

            if i.Q < i.Qlim[0]:
                print('---------------------------------------------------------------------------')
                print('!!!!!!!!!!!!!!!!! Lower Q limit violated at bus {}!!!!!!!!!!!!!!!!!'.format(i.num))
                i.type = 'PQ'
                i.Q_spec = i.Qlim[0]


def T(i:Bus,j:Bus, y):
    """
    Simplification of load flow equations
    :param i: frombus
    :param j: tobus
    :param y: line admittance
    :return: float
    """
    return np.real(y)*np.cos(i.angle-j.angle)+np.imag(y)*np.sin(i.angle-j.angle)


def U(i:Bus,j:Bus,y):
    """
    Simplification of load flow equations
    :param i: frombus
    :param j: tobus
    :param y: line admittance
    :return: float
    """
    return np.real(y) * np.sin(i.angle - j.angle) - np.imag(y) * np.cos(
        i.angle - j.angle)


def calculate_power(bd:List[Bus]):
    """
    Calculates and sets power for each Bus
    :param bd: List of Bus elements
    :return: none
    """
    for i in bd:
        i.P = 0
        i.Q = 0
        self_y = 1j*i.comp
        for bus, line in i.connected_lines.items():
            self_y += line.y
            i.P -= i.V*bus.V*T(i,bus,line.y)
            i.Q -= i.V*bus.V*U(i,bus,line.y)
        i.P += i.V**2*np.real(self_y)
        i.Q -= i.V**2*np.imag(self_y)

    return None


def calc_jacobian(bd:List[Bus]):
    """
    "Meat" of the load flow solver. Computes the linearised power flow
    equations, popularly known as the Jacobian.
    The complexity of the code stems from the indices issue: bus indices
    are not the same as Jacobian indices.
    :param bd: List of Bus elements
    :return: Jacobian
    """
    j_ind1 = 0
    j_ind2 = 0
    for i in bd:
        if i.type == 'PQ':
            i.j_ind1 = j_ind1
            i.j_ind2 = j_ind2
            j_ind2 += 1
            j_ind1 += 1
        elif i.type == 'PV':
            i.j_ind1 = j_ind1
            j_ind1 += 1
    J1 = np.zeros((j_ind1,j_ind1))
    J2 = np.zeros((j_ind1,j_ind2))
    J3 = np.zeros((j_ind2,j_ind1))
    J4 = np.zeros((j_ind2,j_ind2))

    for i in bd:
        dP1 = 0
        dP2 = 0
        dQ1 = 0
        dQ2 = 0
        self_g = 0
        self_b = i.comp

        #print('{}, type = {}'.format(i.num,i.j_ind1))
        if i.j_ind1 is not None:  # Bus has angle in state variable - means upper part of Jacobian

            for bus, line in i.connected_lines.items():
                self_g += np.real(line.y)
                dP1 += i.V * bus.V * U(i, bus, line.y)
                dP2 -= bus.V * T(i, bus, line.y)
                if bus.j_ind1 is not None:  # Bus has angle in x - this is J1
                    J1[i.j_ind1, bus.j_ind1] = -i.V*bus.V*U(i,bus,line.y)
                if bus.j_ind2 is not None:  # Bus has voltage in x - this is J2
                    J2[i.j_ind1, bus.j_ind2] = -i.V * T(i, bus, line.y)
            dP2 += 2 * i.V * self_g
            J1[i.j_ind1, i.j_ind1] = dP1
            if i.j_ind2 is not None:
                J2[i.j_ind1, i.j_ind2] = dP2

        if i.j_ind2 is not None:  # Bus has voltage in state variable - is lower part
            for bus, line in i.connected_lines.items():
                self_b += np.imag(line.y)
                dQ1 -= i.V * bus.V * T(i, bus, line.y)
                dQ2 -= bus.V * U(i, bus, line.y)
                if bus.j_ind1 is not None:  # Bus has angle in x - this is J3
                    J3[i.j_ind2, bus.j_ind1] = i.V*bus.V*T(i, bus, line.y)
                if bus.j_ind2 is not None:  # Bus has voltage in x - this is J4
                    J4[i.j_ind2, bus.j_ind2] = -i.V * U(i, bus, line.y)
            dQ2 -= 2*i.V*self_b
            J3[i.j_ind2, i.j_ind1] = dQ1
            J4[i.j_ind2, i.j_ind2] = dQ2

    return np.vstack((np.hstack((J1,J2)),np.hstack((J3,J4))))  # Return the full Jacobian matrix


def newton_raphson_iteration(bd:List[Bus],jac):
    """
    A newton-raphson iteration.
    Sets up mismatch vector and solves (not by inversion) the linearised
    equation sets. Finally it updates the state variables.
    :param bd: List of Bus elements
    :param jac: Jacobian matrix
    :return: None
    """
    dP = []
    dQ = []
    a = 0
    for i in bd:
        if i.j_ind1 is not None:
            dP.append(i.P_spec-i.P)
            a += 1
        if i.j_ind2 is not None:
            dQ.append(i.Q_spec-i.Q)

    dS = np.asarray(dP+dQ)
    dx = np.linalg.solve(jac, dS)
    dx = dx
    for i in bd:
        if i.j_ind1 is not None:
            i.angle += dx[i.j_ind1]
        if i.j_ind2 is not None:
            i.V += dx[a+i.j_ind2]


def newton_raphson(bd:List[Bus], set_print = True):
    """
    The iterative scheme of NR load flow.
    #TODO implement convergence check
    :param bd: List of Bus elements
    :return: None
    """
    calculate_power(bd)
    if set_print == True:
        print('---------------------------------------------------------------------------')
        print('Initial conditions')
        printsys(bd)
        print('Initial Jacobian')
        jac = calc_jacobian(bd)
        print(jac)
    for i in range(4):
        pv_limits(bd)
        jac = calc_jacobian(bd)
        #print(jac)
        newton_raphson_iteration(bd, jac)
        calculate_power(bd)

        if set_print == True:
            print('---------------------------------------------------------------------------')
            print('Iteration {}'.format(i + 1))
            printsys(bd)
    if set_print == True:
        powerflows(bd)

    return None


def powerflows(bd:List[Bus]):
    print('---------------------------------------------------------------------------')
    print('Transmission line power flow')
    print('----------------------------')
    print(' Bus indices       P      Q    ')
    Ssum = 0
    for i in bd:
        for bus,line in i.connected_lines.items():
            current = (i.V*np.exp(1j*i.angle)-bus.V*np.exp(1j*bus.angle))*line.y
            i.powerflow[bus] = i.V*np.exp(1j*i.angle)*np.conj(current)
            Ssum += i.powerflow[bus]
        i.printpf()
    print('Total power losses = {}'.format(Ssum))


def check_3bus(R12,X12,R13,X13,P1,V1,P2,Q2,myangles,myVs,set_print = False):
    bus0 = Bus([0,'SB',0,0,1,-100,100,-100,100,0])
    bus1 = Bus([1, 'PV', P1, 0, V1, -100, 100, -100, 100,0])
    bus2 = Bus([2, 'PQ', P2, Q2, 1,-100, 100, -100, 100,0])


    Line1 =  Line([0,1,R12,X12])
    Line2 = Line([0,2,R13,X13])

    bus0.connected_lines[bus1] = Line1
    bus1.connected_lines[bus0] = Line1
    bus0.connected_lines[bus2] = Line2
    bus2.connected_lines[bus0] = Line2

    bus_data = [bus0,bus1,bus2]
    newton_raphson(bus_data, set_print=False)
    angles = [i.angle for i in bus_data]
    Vs = [i.V for i in bus_data]

    if set_print:
        print('')
        print('       Bus 1    Bus 2    Bus 3')
        print("  V:   {:0.2f}     {:0.2f}     {:0.2f}".format(Vs[0], Vs[1], Vs[2]))
        print("myV:   {:0.2f}     {:0.2f}     {:0.2f}".format(myVs[0], myVs[1], myVs[2]))
        print("  d:   {:0.2f}     {:0.2f}     {:0.2f}".format(angles[0], angles[1], angles[2]))
        print("myd:   {:0.2f}     {:0.2f}     {:0.2f}".format(myangles[0], myangles[1], myangles[2]))
    tol = 0.01

    for i in range(0,3):
        if angles[i]-myangles[i]>tol:
            print('Fail')
            return
        if Vs[i]-myVs[i]>tol:
            print('Fail')
            return
    print('Success!')
