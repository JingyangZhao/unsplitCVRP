import math
from gurobipy import *
import sys
import time
from fractions import Fraction
################################################################################################################################
### generalflag = 1: general capacity; generalflag = 0: fixed capacity.
### N = 40 gives the 40 x 20 grid.

def singlelp(p, delta, generalflag, N):
    p = float(p)
    delta = float(delta)
    generalflag = int(generalflag)
    N = int(N)
    K = N // 2
    threshold_index = round(delta * N); delta = threshold_index / N
    ############################################################################################################################
    m = Model('calculate ratio')
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
    ### Cases s in R: (), (d1), (d1,d2), (d1,d2,d3^q), (d1,d2,d3^q,dr).
    ### Since theta = delta, the only intervals are [delta,1/2] and [1/2,1].
    interval_indices = (threshold_index, K) if threshold_index < K else (K,)
    cases = [(0, 0, ())]

    # d1, d2, d3, dr are interval indices; their lower endpoints are d1/N, etc.
    for d1 in interval_indices:
        cases.append((1, 0, (d1,)))
        for d2 in interval_indices:
            if d2 > d1 or d1 + d2 > N:
                continue
            cases.append((2, 0, (d1, d2)))
            for d3 in interval_indices:
                if d3 > d2:
                    continue
                qmax = min((N-1)//threshold_index - 2, (N-d1-d2)//d3)
                for q in range(1, qmax + 1):
                    cases.append((3, q, (d1, d2, d3)))
                for dr in interval_indices:
                    if dr > d3:
                        continue
                    qmax = min((N-1)//threshold_index - 3, (N-d1-d2-dr)//d3)
                    for q in range(1, qmax + 1):
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
            upper[s,i] = 0.5 if h < K else 1.0
    ############################################################################################################################
    ### For each pair of demand intervals, determine whether
    ### (f_delta(d_u)+f_delta(d_v)-1)_+ is always zero, already affine, or can cross zero.
    positive_pairs = set()
    ambiguous_pairs = set()
    delta_exact = Fraction(threshold_index, N)

    for u in interval_indices:
        for v in interval_indices:
            if v > u:
                continue

            du_low = Fraction(u, N)
            du_high = Fraction(1, 2) if u < K else Fraction(1)
            dv_low = Fraction(v, N)
            dv_high = Fraction(1, 2) if v < K else Fraction(1)

            if u < threshold_index:
                fu_low = du_low / (1 - delta_exact)
                fu_high = du_high / (1 - delta_exact)
            elif u < K:
                fu_low = (2 * du_low - delta_exact) / (1 - delta_exact)
                fu_high = (2 * du_high - delta_exact) / (1 - delta_exact)
            else:
                fu_low = fu_high = Fraction(1)

            if v < threshold_index:
                fv_low = dv_low / (1 - delta_exact)
                fv_high = dv_high / (1 - delta_exact)
            elif v < K:
                fv_low = (2 * dv_low - delta_exact) / (1 - delta_exact)
                fv_high = (2 * dv_high - delta_exact) / (1 - delta_exact)
            else:
                fv_low = fv_high = Fraction(1)

            ### Since f_delta is nondecreasing, these are the minimum and maximum
            ### of f_delta(d_u)+f_delta(d_v)-1 over the two demand intervals.
            low = fu_low + fv_low - 1
            high = fu_high + fv_high - 1

            if high <= 0:
                continue
            elif low >= 0:
                positive_pairs.add((u,v))
            else:
                ambiguous_pairs.add((u,v))

    positive_B = []
    positive_C = []
    ambiguous_B = []
    ambiguous_C = []

    for s in R:
        if n[s] >= 2:
            pair = interval[s][0], interval[s][1]
            if pair in positive_pairs:
                positive_B.append(s)
            if pair in ambiguous_pairs:
                ambiguous_B.append(s)

        if n[s] >= 3:
            pair = interval[s][1], interval[s][2]
            if pair in positive_pairs:
                positive_C.append(s)
            if pair in ambiguous_pairs:
                ambiguous_C.append(s)
    ############################################################################################################################
    ### Variables in the single-setting LP (B): eta_s, x_{s,i}, B_s, C_s.
    eta = m.addVars(R, vtype=GRB.CONTINUOUS, lb=0.0, ub=1.0, name='eta')
    x = m.addVars(((s,i) for s in R for i in range(1, n[s] + 1)), vtype=GRB.CONTINUOUS, lb=0.0, ub=1.0, name='x')
    B = m.addVars(ambiguous_B, vtype=GRB.CONTINUOUS, lb=0.0, ub=1.0, name='B')
    C = m.addVars(ambiguous_C, vtype=GRB.CONTINUOUS, lb=0.0, ub=1.0, name='C')
    ############################################################################################################################
    ### Auxiliary sum in LP (B): price_delta = sum_s sum_i m_{s,i} w_{s,i}
    price_delta = m.addVar(vtype=GRB.CONTINUOUS, lb=0.0, ub=GRB.INFINITY, name='price_delta')

    price_expr = LinExpr()

    for s in R:
        for i in range(1, n[s] + 1):
            t = interval[s][i-1]

            if t < threshold_index:
                w = x[s,i] / (1-delta)
            elif t < K:
                w = (2*x[s,i] - delta*eta[s]) / (1-delta)
            else:
                w = eta[s]

            price_expr += multiplicity[s,i]*w

    m.addConstr(price_delta == price_expr)

    ############################################################################################################################
    ### Auxiliary sum in LP (B): demand_sum = sum_s sum_i m_{s,i} x_{s,i}.
    demand_sum = m.addVar(vtype=GRB.CONTINUOUS, lb=0.0, ub=1.0, name='demand_sum')
    m.addConstr(
        demand_sum == quicksum(
            multiplicity[s,i]*x[s,i] for s in R for i in range(1, n[s] + 1)
        )
    )
    ############################################################################################################################
    ### Auxiliary sum in LP (B): B_delta = sum_s B_s.
    B_delta = m.addVar(vtype=GRB.CONTINUOUS, lb=0.0, ub=1.0, name='B_delta')

    saving_expr = LinExpr()

    for s in ambiguous_B:
        saving_expr += B[s]

    for s in positive_B:
        t1, t2 = interval[s][0], interval[s][1]

        if t1 < threshold_index:
            w1 = x[s,1] / (1-delta)
        elif t1 < K:
            w1 = (2*x[s,1] - delta*eta[s]) / (1-delta)
        else:
            w1 = eta[s]

        if t2 < threshold_index:
            w2 = x[s,2] / (1-delta)
        elif t2 < K:
            w2 = (2*x[s,2] - delta*eta[s]) / (1-delta)
        else:
            w2 = eta[s]

        saving_expr += w1 + w2 - eta[s]

    m.addConstr(B_delta == saving_expr)
    ############################################################################################################################
    ### Auxiliary sum in LP (B): C_delta = sum_s C_s.
    C_delta = m.addVar(vtype=GRB.CONTINUOUS, lb=0.0, ub=1.0, name='C_delta')

    saving_expr = LinExpr()

    for s in ambiguous_C:
        saving_expr += C[s]

    for s in positive_C:
        t2, t3 = interval[s][1], interval[s][2]

        if t2 < threshold_index:
            w2 = x[s,2] / (1-delta)
        elif t2 < K:
            w2 = (2*x[s,2] - delta*eta[s]) / (1-delta)
        else:
            w2 = eta[s]

        if t3 < threshold_index:
            w3 = x[s,3] / (1-delta)
        elif t3 < K:
            w3 = (2*x[s,3] - delta*eta[s]) / (1-delta)
        else:
            w3 = eta[s]

        saving_expr += w2 + w3 - eta[s]

    m.addConstr(C_delta == saving_expr)
    ############################################################################################################################
    ### LP (B): sum_s eta_s = 1.
    m.addConstr(quicksum(eta[s] for s in R) == 1)
    ############################################################################################################################
    ### LP (B): sum_i m_{s,i} x_{s,i} <= eta_s.
    for s in R:
        m.addConstr(
            quicksum(multiplicity[s,i]*x[s,i] for i in range(1, n[s] + 1)) <= eta[s]
        )
    ############################################################################################################################
    ### LP (B): ell_{s,i} eta_s <= x_{s,i} <= u_{s,i} eta_s.
    for s in R:
        for i in range(1, n[s] + 1):
            m.addConstr(x[s,i] >= ell[s,i]*eta[s])
            m.addConstr(x[s,i] <= upper[s,i]*eta[s])
    ############################################################################################################################
    ### LP (B): x_{s,i} >= x_{s,i+1}.
    for s in R:
        for i in range(1, n[s]):
            m.addConstr(x[s,i] >= x[s,i+1])
    ############################################################################################################################
    ### LP (B): B_s >= w_{s,1}+w_{s,2}-eta_s.
    for s in B.keys():
        t1, t2 = interval[s][0], interval[s][1]

        if t1 < threshold_index:
            w1 = x[s,1] / (1-delta)
        elif t1 < K:
            w1 = (2*x[s,1] - delta*eta[s]) / (1-delta)
        else:
            w1 = eta[s]

        if t2 < threshold_index:
            w2 = x[s,2] / (1-delta)
        elif t2 < K:
            w2 = (2*x[s,2] - delta*eta[s]) / (1-delta)
        else:
            w2 = eta[s]

        m.addConstr(B[s] >= w1 + w2 - eta[s])
    ############################################################################################################################
    ### LP (B): C_s >= w_{s,2}+w_{s,3}-eta_s.
    for s in C.keys():
        t2, t3 = interval[s][1], interval[s][2]

        if t2 < threshold_index:
            w2 = x[s,2] / (1-delta)
        elif t2 < K:
            w2 = (2*x[s,2] - delta*eta[s]) / (1-delta)
        else:
            w2 = eta[s]

        if t3 < threshold_index:
            w3 = x[s,3] / (1-delta)
        elif t3 < K:
            w3 = (2*x[s,3] - delta*eta[s]) / (1-delta)
        else:
            w3 = eta[s]

        m.addConstr(C[s] >= w2 + w3 - eta[s])
    ############################################################################################################################
    ### LP (B): 0 <= B_s, C_s <= eta_s. 
    for s in B.keys():
        m.addConstr(B[s] <= eta[s])

    for s in C.keys():
        m.addConstr(C[s] <= eta[s])
    ############################################################################################################################
    xi = 1.0 if generalflag == 1 else p
    kappa = min(3*p*p - 2*p*p*p, p**1.5)
    m.setObjective(
        -math.log(p) + p*price_delta + xi*(1-demand_sum)/(1-delta)
        - p*p*B_delta - (kappa-p*p)*C_delta, GRB.MAXIMIZE
    ); m.optimize()
    if m.status == GRB.OPTIMAL:
        return m.objVal
    else: return -1


if __name__ == '__main__':
    start = time.time()
    generalflag = int(sys.argv[1])
    N = int(sys.argv[2])

    best_value = math.inf
    best_p = None
    best_delta = None

    for threshold_index in range(1, N // 2 + 1):
        delta = threshold_index / N
        for ip in range(1, N + 1):
            p = ip / N
            value = singlelp(p, delta, generalflag, N)
            if value < best_value:
                best_value = value
                best_p = p
                best_delta = delta

    value, p, delta = best_value, best_p, best_delta
    result = f'generalflag: {generalflag}, N: {N}, value: {value:.12f}, p: {p:.12g}, delta: {delta:.12g}, runtime: {time.time()-start:.6f} seconds'
    
    print(result)
