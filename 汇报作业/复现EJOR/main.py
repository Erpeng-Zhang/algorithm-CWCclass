import random
import math

# Parameters
w = 1 / (2 * math.log(2))  # Inertia weight
c1 = 0.5 + math.log(2)  # Individual learning factor
c2 = 0.5 + math.log(2)  # Swarm learning factor
V_min = -2  # Minimum velocity
V_max = 2  # Maximum velocity
S_max = 40  # Swarm size

TI = 100  # Maximum number of iterations
TF = 0.01  # Minimum fitness value
R = 50  # Number of repetitions
alpha = 0.9999  # Annealing coefficient

# Problem data
N = 12  # Number of tasks
M = 3  # Number of machines
V = 2  # Number of vehicles
data = []  # Task data - task ID, machine ID, processing time

LB = cal_LB()  # Lower bound of feasible solutions

# Transportation time data
time = []

# Initialize particle swarm
particles = [Model() for _ in range(S_max)]

for i in range(S_max):
    particles[i].x = [random.random() for _ in range(2 * N)]
    particles[i].y = fitness(particles[i].x, data, time)

# Initialize individual best
Pbest = particles.copy()

# Find the index of the initial global best
[_, a] = min([particle.y for particle in particles])

# Initialize global best
Gbest = particles[a]

X_SA = particles[a].x
F_SA = particles[a].y

# Initialize temperature and iteration count
T = TI
g = 0

# Main loop
while T > TF or Gbest.y > LB:
    for r in range(R):
        # TODO: Generate neighboring solutions
        Y = cal_larboslusion(X_SA)
        F_Y = fitness(Y, data, time)

        # Calculate fitness change
        delta_F = F_SA - F_Y

        if delta_F > 0:
            X_SA = Y
            F_SA = F_Y
        else:
            if random.random() < math.exp(delta_F / T):
                X_SA = Y
                F_SA = F_Y

                # Update global best
                if F_SA < Gbest.y:
                    Gbest.x = X_SA
                    Gbest.y = F_SA
            else:
                # Calculate acceptance probability of new solution
                Pr = math.exp(delta_F / T)

                if Pr > random.random():
                    X_SA = Y
                    F_SA = F_Y

    # Update particles
    for s in range(S_max):
        # TODO: Update particle velocity
        particles[s].v = w * particles[s].v + c1 * random.random() * (Pbest[s].x - particles[s].x) + c2 * random.random() * (Gbest.x - particles[s].x)

        # TODO: Update particle position
        particles[s].x = [sum(x) for x in zip(particles[s].x, particles[s].v)]

        # Update particle fitness
        particles[s].y = fitness(particles[s].x, data, time)

        # Update individual best
        if particles[s].y < Pbest[s].y:
            Pbest[s].x = particles[s].x
            Pbest[s].y = particles[s].y

    # Update global best
    [_, a] = min([particle.y for particle in particles])

    if particles[a].y < Gbest.y:
        Gbest.x = particles[a].x
        Gbest.y = particles[a].y

    # Update X_SA and F_SA
    X_SA = Gbest.x
    F_SA = Gbest.y

    # Update temperature and iteration count
    T = alpha * T
    g += 1

    # Print the best solution at each iteration
    print("Iteration:", g, ", Current temperature:", T, ", Best solution:", Gbest.y)
    plt.plot(g, Gbest.y)
    plt.hold(True)

# TODO: Fitness function
def fitness(x, data, time):
    y
