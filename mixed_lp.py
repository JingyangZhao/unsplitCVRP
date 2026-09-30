import math
from gurobipy import *
import sys
import time
from fractions import Fraction
################################################################################################################################
### generalflag = 1: general capacity; generalflag = 0: fixed capacity.
### theta = 1/N.  N = 40 gives the 40 x 20 parameter grid.

def mixlp(generalflag, N):
    generalflag = int(generalflag)
    N = int(N)
    K = N // 2
    ############################################################################################################################
    m = Model('calculate ratio')
    y = m.addVar(vtype=GRB.CONTINUOUS, lb=-GRB.INFINITY, name='y')
    # m.setParam('OutputFlag', 0)
    m.setParam('Threads', 4)
    m.setParam('Method', 2)
    m.setParam('BarOrder', 1)
    m.setParam('BarHomogeneous', 1)
    m.setParam('NumericFocus', 2)
    m.setParam('Crossover', 0)
    # m.setParam('BarConvTol', 1e-12)
    # m.setParam('FeasibilityTol', 1e-8)
    # m.setParam('OptimalityTol', 1e-8)
    # m.setParam('TimeLimit', 1800)
    ############################################################################################################################
    ### Parameter choices j in J.  Here |J| = N*(N/2), so N = 40 gives the 40 x 20 grid.
    ### p[j], delta[j], and xi[j] correspond to p_j, delta_j, and xi_j in the paper.
    J = []
    p = {}
    delta = {}
    delta_index = {}
    xi = {}

    j = 0
    for ip in range(1, N + 1):
        for h in range(1, K + 1):
            J.append(j)
            p[j] = ip / N
            delta[j] = h / N
            delta_index[j] = h
            xi[j] = 1.0 if generalflag == 1 else p[j]
            j += 1

    delta_indices = range(1, K + 1)  # distinct values delta = h/N; J itself indexes all parameter pairs
    ############################################################################################################################
    ### Cases s in R: (), (d1), (d1,d2), (d1,d2,d3^q), (d1,d2,d3^q,dr).
    ### A threshold index h < K represents [h/N,(h+1)/N], and h = K represents [1/2,1].
    cases = [(0, 0, ())]

    for d1 in delta_indices:
        cases.append((1, 0, (d1,)))
        for d2 in range(1, d1 + 1):
            if d1 + d2 > N:
                continue
            cases.append((2, 0, (d1, d2)))
            for d3 in range(1, d2 + 1):
                for q in range(1, min(N - 3, (N - d1 - d2) // d3) + 1):
                    cases.append((3, q, (d1, d2, d3)))
                for dr in range(1, d3 + 1):
                    for q in range(1, min(N - 4, (N - d1 - d2 - dr) // d3) + 1):
                        cases.append((4, q, (d1, d2, d3, dr)))

    R = range(len(cases))
    ############################################################################################################################
    ### n_s, m_{s,i}, ell_{s,i}, and u_{s,i} in the paper.
    n = {}
    multiplicity = {}
    interval = {}
    ell = {}
    upper = {}

    for s in R:
        n[s], q, interval[s] = cases[s]

        if n[s] == 0:
            mult = ()
        elif n[s] == 1:
            mult = (1,)
        elif n[s] == 2:
            mult = (1, 1)
        elif n[s] == 3:
            mult = (1, 1, q)
        else:
            mult = (1, 1, q, 1)

        for i in range(1, n[s] + 1):
            multiplicity[s,i] = mult[i-1]
            h = interval[s][i-1]
            ell[s,i] = h / N
            upper[s,i] = (h + 1) / N if h < K else 1.0
    ############################################################################################################################
    ### For each pair of demand intervals, determine whether
    ### (f_delta(d_u)+f_delta(d_v)-1)_+ is always zero, already affine, or can cross zero.
    positive_delta = {}
    ambiguous_delta = {}

    for u in delta_indices:
        for v in range(1, u + 1):
            positive_delta[u,v] = []
            ambiguous_delta[u,v] = []

            du_low = Fraction(u, N)
            du_high = Fraction(u + 1, N) if u < K else Fraction(1)
            dv_low = Fraction(v, N)
            dv_high = Fraction(v + 1, N) if v < K else Fraction(1)

            for h in delta_indices:
                delta_h = Fraction(h, N)

                if u < h:
                    fu_low = du_low / (1 - delta_h)
                    fu_high = du_high / (1 - delta_h)
                elif u < K:
                    fu_low = (2 * du_low - delta_h) / (1 - delta_h)
                    fu_high = (2 * du_high - delta_h) / (1 - delta_h)
                else:
                    fu_low = fu_high = Fraction(1)

                if v < h:
                    fv_low = dv_low / (1 - delta_h)
                    fv_high = dv_high / (1 - delta_h)
                elif v < K:
                    fv_low = (2 * dv_low - delta_h) / (1 - delta_h)
                    fv_high = (2 * dv_high - delta_h) / (1 - delta_h)
                else:
                    fv_low = fv_high = Fraction(1)

                ### Since f_delta is nondecreasing, these are the minimum and maximum
                ### of f_delta(d_u)+f_delta(d_v)-1 over the two demand intervals.
                low = fu_low + fv_low - 1
                high = fu_high + fv_high - 1

                if high <= 0:
                    continue
                elif low >= 0:
                    positive_delta[u,v].append(h)
                else:
                    ambiguous_delta[u,v].append(h)

    positive_B = {h: [] for h in delta_indices}
    positive_C = {h: [] for h in delta_indices}
    ambiguous_B = {h: [] for h in delta_indices}
    ambiguous_C = {h: [] for h in delta_indices}

    for s in R:
        if n[s] >= 2:
            u, v = interval[s][0], interval[s][1]
            for h in positive_delta[u,v]:
                positive_B[h].append(s)
            for h in ambiguous_delta[u,v]:
                ambiguous_B[h].append(s)

        if n[s] >= 3:
            u, v = interval[s][1], interval[s][2]
            for h in positive_delta[u,v]:
                positive_C[h].append(s)
            for h in ambiguous_delta[u,v]:
                ambiguous_C[h].append(s)
    ############################################################################################################################
    ### Variables in LP (C): y, eta_s, x_{s,i}, B_{j,s}, C_{j,s}.
    eta = m.addVars(R, vtype=GRB.CONTINUOUS, lb=0.0, ub=1.0, name='eta')
    x = m.addVars(((s,i) for s in R for i in range(1, n[s] + 1)), vtype=GRB.CONTINUOUS, lb=0.0, ub=1.0, name='x')
    B = m.addVars(((h,s) for h in delta_indices for s in ambiguous_B[h]), vtype=GRB.CONTINUOUS, lb=0.0, ub=1.0, name='B')
    C = m.addVars(((h,s) for h in delta_indices for s in ambiguous_C[h]), vtype=GRB.CONTINUOUS, lb=0.0, ub=1.0, name='C')
    ############################################################################################################################
    ### Auxiliary sum in LP (C): price_delta[h] = sum_s sum_i m_{s,i} w_{j,s,i}
    price_delta = m.addVars(delta_indices, vtype=GRB.CONTINUOUS, lb=0.0, ub=GRB.INFINITY, name='price_delta')

    for h in delta_indices:
        delta_h = h / N
        price_expr = LinExpr()

        for s in R:
            for i in range(1, n[s] + 1):
                t = interval[s][i-1]

                if t < h:
                    w = x[s,i] / (1-delta_h)
                elif t < K:
                    w = (2*x[s,i] - delta_h*eta[s]) / (1-delta_h)
                else:
                    w = eta[s]

                price_expr += multiplicity[s,i]*w

        m.addConstr(price_delta[h] == price_expr)
    ############################################################################################################################
    ### Auxiliary sum in LP (C): demand_sum = sum_s sum_i m_{s,i} x_{s,i}.
    ### It is introduced only to avoid repeating the same large sum in all |J| cost constraints.
    demand_sum = m.addVar(vtype=GRB.CONTINUOUS, lb=0.0, ub=1.0, name='demand_sum')
    m.addConstr(
        demand_sum == quicksum(
            multiplicity[s,i]*x[s,i] for s in R for i in range(1, n[s] + 1)
        )
    )
    ############################################################################################################################
    ### Auxiliary sum in LP (C): B_delta[h] = sum_s B_{j,s} for any j with delta_j = h/N.
    B_delta = m.addVars(delta_indices, vtype=GRB.CONTINUOUS, lb=0.0, ub=1.0, name='B_delta')

    for h in delta_indices:
        delta_h = h / N
        saving_expr = LinExpr()

        for s in ambiguous_B[h]:
            saving_expr += B[h,s]

        for s in positive_B[h]:
            t1, t2 = interval[s][0], interval[s][1]

            if t1 < h:
                w1 = x[s,1] / (1-delta_h)
            elif t1 < K:
                w1 = (2*x[s,1] - delta_h*eta[s]) / (1-delta_h)
            else:
                w1 = eta[s]

            if t2 < h:
                w2 = x[s,2] / (1-delta_h)
            elif t2 < K:
                w2 = (2*x[s,2] - delta_h*eta[s]) / (1-delta_h)
            else:
                w2 = eta[s]

            saving_expr += w1 + w2 - eta[s]

        m.addConstr(B_delta[h] == saving_expr)
    ############################################################################################################################
    ### Auxiliary sum in LP (C): C_delta[h] = sum_s C_{j,s} for any j with delta_j = h/N.
    C_delta = m.addVars(delta_indices, vtype=GRB.CONTINUOUS, lb=0.0, ub=1.0, name='C_delta')

    for h in delta_indices:
        delta_h = h / N
        saving_expr = LinExpr()

        for s in ambiguous_C[h]:
            saving_expr += C[h,s]

        for s in positive_C[h]:
            t2, t3 = interval[s][1], interval[s][2]

            if t2 < h:
                w2 = x[s,2] / (1-delta_h)
            elif t2 < K:
                w2 = (2*x[s,2] - delta_h*eta[s]) / (1-delta_h)
            else:
                w2 = eta[s]

            if t3 < h:
                w3 = x[s,3] / (1-delta_h)
            elif t3 < K:
                w3 = (2*x[s,3] - delta_h*eta[s]) / (1-delta_h)
            else:
                w3 = eta[s]

            saving_expr += w2 + w3 - eta[s]

        m.addConstr(C_delta[h] == saving_expr)
    ############################################################################################################################
    ### LP (C): y <= sum_s G_{j,s}, for every j in J.
    for j in J:
        h = delta_index[j]
        kappa = min(3*p[j]*p[j] - 2*p[j]*p[j]*p[j], p[j]**1.5)

        m.addConstr(
            y <= -math.log(p[j]) + p[j]*price_delta[h] + xi[j]*(1-demand_sum)/(1-delta[j]) - p[j]*p[j]*B_delta[h] - (kappa-p[j]*p[j])*C_delta[h]
        )
    ############################################################################################################################
    ### LP (C): sum_s eta_s = 1.
    m.addConstr(quicksum(eta[s] for s in R) == 1)
    ############################################################################################################################
    ### LP (C): sum_i m_{s,i} x_{s,i} <= eta_s.
    for s in R:
        m.addConstr(
            quicksum(multiplicity[s,i]*x[s,i] for i in range(1, n[s] + 1)) <= eta[s]
        )
    ############################################################################################################################
    ### LP (C): ell_{s,i} eta_s <= x_{s,i} <= u_{s,i} eta_s.
    for s in R:
        for i in range(1, n[s] + 1):
            m.addConstr(x[s,i] >= ell[s,i]*eta[s])
            m.addConstr(x[s,i] <= upper[s,i]*eta[s])
    ############################################################################################################################
    ### LP (C): x_{s,i} >= x_{s,i+1}.
    for s in R:
        for i in range(1, n[s]):
            m.addConstr(x[s,i] >= x[s,i+1])
    ############################################################################################################################
    ### LP (C): B_{j,s} >= w_{j,s,1}+w_{j,s,2}-eta_s.
    for h, s in B.keys():
        delta_h = h / N
        t1, t2 = interval[s][0], interval[s][1]

        if t1 < h:
            w1 = x[s,1] / (1-delta_h)
        elif t1 < K:
            w1 = (2*x[s,1] - delta_h*eta[s]) / (1-delta_h)
        else:
            w1 = eta[s]

        if t2 < h:
            w2 = x[s,2] / (1-delta_h)
        elif t2 < K:
            w2 = (2*x[s,2] - delta_h*eta[s]) / (1-delta_h)
        else:
            w2 = eta[s]

        m.addConstr(B[h,s] >= w1 + w2 - eta[s])
    ############################################################################################################################
    ### LP (C): C_{j,s} >= w_{j,s,2}+w_{j,s,3}-eta_s.
    for h, s in C.keys():
        delta_h = h / N
        t2, t3 = interval[s][1], interval[s][2]

        if t2 < h:
            w2 = x[s,2] / (1-delta_h)
        elif t2 < K:
            w2 = (2*x[s,2] - delta_h*eta[s]) / (1-delta_h)
        else:
            w2 = eta[s]

        if t3 < h:
            w3 = x[s,3] / (1-delta_h)
        elif t3 < K:
            w3 = (2*x[s,3] - delta_h*eta[s]) / (1-delta_h)
        else:
            w3 = eta[s]

        m.addConstr(C[h,s] >= w2 + w3 - eta[s])
    ############################################################################################################################
    ### LP (C): 0 <= B_{j,s}, C_{j,s} <= eta_s.
    for h, s in B.keys():
        m.addConstr(B[h,s] <= eta[s])

    for h, s in C.keys():
        m.addConstr(C[h,s] <= eta[s])
    ############################################################################################################################
    m.setObjective(y, GRB.MAXIMIZE); m.optimize()
    if m.status == GRB.OPTIMAL:
        return m.objVal
    else: return -1


if __name__ == '__main__':
    start = time.time()
    generalflag = sys.argv[1]
    N = sys.argv[2]

    value = mixlp(generalflag, N)
    result = f'generalflag: {generalflag}, N: {N}, value: {value}, runtime: {time.time()-start} seconds'

    print(result)
