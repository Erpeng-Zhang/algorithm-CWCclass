import pulp


# Create a new LP problem
prob = pulp.LpProblem("Job_Shop_Scheduling_Problem", pulp.LpMinimize)

# Sets
J = [...]  # Set of jobs
M = [...]  # Set of machines
V = [...]  # Set of vehicles
f = []
t = []
J_j = [...]  
# n_j + 1
J1_j = [...]  

# Parameters
p_kl = {...}  # Processing time of job k on machine l
tau_kl = {...}  # Transport time from operation (k-1) to k for job l
tau_LU = {...}  # Transport time from LU to job's first machine
v_ij = {...}  # Time vehicle needs to travel between jobs i and j
A = ...  # Total number of vehicles
LN = ...  # A large number

# Decision Variables
w = pulp.LpVariable.dicts("w", ((k, l, i, j) for k in J for l in M for i in J for j in M), 0, 1, pulp.LpBinary)
x = pulp.LpVariable.dicts("x", ((k, l, i, j) for k in J for l in M for i in J for j in M), 0, 1, pulp.LpBinary)
z = pulp.LpVariable.dicts("z", ((i, j) for i in J for j in M), cat='Binary')
u = pulp.LpVariable.dicts("u", ((i, j) for i in J for j in M), cat='Binary')
c = pulp.LpVariable.dicts("c", ((k, l) for k in J for l in M), lowBound=0)
T = pulp.LpVariable("T", lowBound=0)
v = pulp.LpVariable.dicts("v", ((k, l) for k in J for l in M), lowBound=0)

# Objective Function
prob += T

# Constraints
# 1. T is at least the completion time of the last operation
for l in M:
    for k in J:
        prob += T >= c[k, l]

# 2. Each operation is followed by exactly one other operation on the same machine
for l in M:
    for k in J:
        prob += sum(x[k, l, i, j] for i in J for j in M) == 1
        prob += sum(x[i, j, k, l] for i in J for j in M) == 1

# 3. Ensure the same number of first and last tasks for vehicles
for j in J:
    prob += sum(z[i, j] for i in J) == sum(u[i, j] for i in J)

# 4. At most A vehicles available
prob += sum(u[i, j] for i in J for j in M) <= A

# 5. Each task followed and preceded by exactly one other task
for l in M:
    for k in J:
        prob += sum(x[k, l, i, j] for i in J for j in M) == 1
        prob += sum(x[i, j, k, l] for i in J for j in M) == 1

# 6. Completion time constraints
for l in M:
    for k in J:
        prob += c[k, l] >= c[(k-1), l] + p_kl[k, l] + tau_kl[k, (k-1)]
        prob += c[k, l] >= c[(k-1), l] + p_kl[k, l] + tau_LU[k, l]
        prob += c[k, l] >= c[(k-1), l] + p_kl[k, l] + v_ij[i, j]

# 7. Vehicle time constraints
for i in J:
    for j in M:
        prob += v[i, j] >= v[i, j] + tau_kl[(k-1), l] + tau_LU[k, l]

# 8. Define nature of variables
for l in M:
    for k in J:
        prob += x[k, l, i, j] <= 1
        prob += x[i, j, k, l] <= 1

# 9. Non-negative variables
for l in M:
    for k in J:
        prob += c[k, l] >= 0
        prob += v[k, l] >= 0

# Solve the problem
prob.solve()

# Print the results
print(f"Status: {pulp.LpStatus[prob.status]}")
print(f"Objective value: {pulp.value(prob.objective)}")
for v in prob.variables():
    print(f"{v.name} = {v.varValue}")
