using Convex, Mosek, MosekTools, LinearAlgebra
# Written in Julia Version 1.4.2
# Mosek and Mosektools require a (free) license to be installed.
# Alternatively, one could use the solver SCS

# The code used to obtain the results from the thesis are given below
# If you want to run some code for yourself, you can un-comment a section
# by using CTRL + /

# Variables:
# A = number of answers Alice
# B = number of answers Bob
# S = number of questions Alice
# T = number of questions Bob

# P = given bipartite correlation
# d = entanglement dimension

# E[s][a] is complex dxd matrix for all s in S and a in A
# F[t][b] is complex dxd matrix for all t in T and b in B
# psi = complex unit vector with dimension d^2
# rho = complex (d^2)x(d^2) matrix with trace 1

function findE(F, rho, P::Array,d,A,B,S,T)
    #Solves SDP for optimal E, given state rho, POVM F and desired P
    E = [[ComplexVariable(d,d) for s=1:S] for a=1:A]
    M = Variable(1, Positive())
    objective = M
    s_1 = [[[[Variable(1, Positive()) for t=1:T] for s=1:S] for b=1:B] for a=1:A]
    s_2 = [[[[Variable(1, Positive()) for t=1:T] for s=1:S] for b=1:B] for a=1:A]

    constraints = [E[a][s] in :SDP for s=1:S for a=1:A]
    constraints += [P[a][b][s][t] - tr(kron(E[a][s], F[b][t])*rho) == M - s_1[a][b][s][t] for a=1:A for s=1:S for b=1:B for t=1:T]
    constraints += [P[a][b][s][t] - tr(kron(E[a][s], F[b][t])*rho) == -M + s_2[a][b][s][t] for a=1:A for s=1:S for b=1:B for t=1:T]

    ide = Matrix{Float64}(I, d, d)
    for s in 1:S
        constraints += sum(E[a][s] for a=1:A) == ide
    end

    problem = minimize(objective,constraints)
    solve!(problem, () -> Mosek.Optimizer())

    E = [[E[a][s].value for s=1:S] for a=1:A]

end

function findF(E, rho, P::Array,d,A,B,S,T)
    #Solves SDP for optimal F, given density matrix rho, POVM E and desired P
    F = [[ComplexVariable(d,d) for t=1:T] for b=1:B]
    M = Variable(1, Positive())
    objective = M
    s_1 = [[[[Variable(1, Positive()) for t=1:T] for s=1:S] for b=1:B] for a=1:A]
    s_2 = [[[[Variable(1, Positive()) for t=1:T] for s=1:S] for b=1:B] for a=1:A]

    constraints = [F[b][t] in :SDP for t=1:T for b=1:B]

    if d==1
        # Of course, this doesn't seem neccessary, but there seems to be
        # a bug with the trace when d=1
        constraints += [P[a][b][s][t] - E[a][s]* F[b][t]*rho == M - s_1[a][b][s][t] for a=1:A for s=1:S for b=1:B for t=1:T]
        constraints += [P[a][b][s][t] - E[a][s]* F[b][t]*rho == -M + s_2[a][b][s][t] for a=1:A for s=1:S for b=1:B for t=1:T]

    else
        constraints += [P[a][b][s][t] - tr(kron(E[a][s], F[b][t])*rho) == M - s_1[a][b][s][t] for a=1:A for s=1:S for b=1:B for t=1:T]
        constraints += [P[a][b][s][t] - tr(kron(E[a][s], F[b][t])*rho) == -M + s_2[a][b][s][t] for a=1:A for s=1:S for b=1:B for t=1:T]
    end

    ide = Matrix{Float64}(I, d, d)
    for t in 1:T
        constraints += sum(F[b][t] for b=1:B) == ide
    end

    problem = minimize(objective,constraints)
    solve!(problem, () -> Mosek.Optimizer())

    F = [[F[b][t].value for t=1:T] for b=1:B]
end

function generaterandomrho(d::Int)
    # generates random density matrix
    psi = rand(Float64,(d^2,1)) + im*rand(Float64,(d^2,1))
    psi /= norm(psi, 2)
    rho = psi * psi'
end

function randomE(d,A,S)
    # generates random POVM E, by generating random matrices and solving for
    # the closest SDP matrices

    E = [[rand(Float64,(d,d))+ im*rand(Float64,(d,d)) for s=1:S] for a=1:A]

    E_SDP= [[ComplexVariable(d,d) for s=1:S] for a=1:A]
    M = Variable(1)
    objective = M

    s_1 = [[Variable(1, Positive())  for s=1:S] for a=1:A]
    s_2 = [[Variable(1, Positive()) for s=1:S] for a=1:A]

    constraints = [E_SDP[a][s] in :SDP for a=1:A for s=1:S]
    constraints += M >= 0

    constraints += [norm(E_SDP[a][s] - E[a][s]) == M - s_1[a][s] for a=1:A for s=1:S]

    constraints += [norm(E_SDP[a][s] - E[a][s]) == -M + s_2[a][s] for a=1:A for s=1:S]

    ide = Matrix{Float64}(I, d, d)
    for s in 1:S
        constraints += sum(E_SDP[a][s] for a=1:A) == ide
    end

    problem = minimize(objective,constraints)
    solve!(problem, () -> Mosek.Optimizer())

# Making sure that all matrices are Hermitian
    E = [[(E_SDP[a][s].value+E_SDP[a][s].value')/2 for s=1:S] for a=1:A]

end

function generaterandomoperators(d,A,B,S,T)
    rhostart = generaterandomrho(d)
    ide = Matrix{Float64}(I, d, d)

    Estart = randomE(d,A,S)
    Fstart = randomE(d,B,T)
    rhostart, Estart, Fstart
end

function generateP(rho, E, F,d,A,B,S,T)
    [[[[tr(kron(E[a][s], F[b][t])*rho) for  t=1:T] for s=1:S] for b=1:B] for a=1:A]
end

function testfindE(d,A,B,S,T)
    rho, E, F= generaterandomoperators(d,A,B,S,T)
    P = generateP(rho, E, F,d,A,B,S,T)
    findE(F, rho, P,d,A,B,S,T)
end

function findrho(E,F,P,d,A,B,S,T)
    # solves SDP for optimal density matrix rho, given POVM's E and F,
    # and desired correlation P
    rho = ComplexVariable(d^2,d^2)
    M = Variable(1, Positive())
    objective = M
    s_1 = [[[[Variable(1, Positive()) for t=1:T] for s=1:S] for b=1:B] for a=1:A]
    s_2 = [[[[Variable(1, Positive()) for t=1:T] for s=1:S] for b=1:B] for a=1:A]

    constraints = rho in :SDP
    constraints += [P[a][b][s][t] - tr(kron(E[a][s], F[b][t])*rho) == M - s_1[a][b][s][t] for a=1:A for s=1:S for b=1:B for t=1:T]
    constraints += [P[a][b][s][t] - tr(kron(E[a][s], F[b][t])*rho) == -M + s_2[a][b][s][t] for a=1:A for s=1:S for b=1:B for t=1:T]
    constraints += tr(rho) == 1

    problem = minimize(objective,constraints)
    solve!(problem, () -> Mosek.Optimizer())
    rho.value
end

function randomP_fixed_d(num_iterations,d,A,B,S,T)
    #executes see-saw algorithm of num_iterations iterations
    #returns a list of errors per step
    rho, E, F = generaterandomoperators(d,A,B,S,T)
    P_real = generateP(rho, E, F,d,A,B,S,T) # We want to approximate this P
    errors = zeros(0)
    rho, E, F = generaterandomoperators(d,A,B,S,T) # We initialize with new rho,E,F
    for i in 1:num_iterations
            myerror = maximum(norm(P_real[a][b][s][t] - generateP(rho,E,F,d,A,B,S,T)[a][b][s][t])  for t=1:T for s=1:S for b=1:B for a=1:A)
            push!(errors,myerror)
        E = findE(F,rho,P_real,d,A,B,S,T)
            myerror = maximum(norm(P_real[a][b][s][t] - generateP(rho,E,F,d,A,B,S,T)[a][b][s][t])  for t=1:T for s=1:S for b=1:B for a=1:A)
            push!(errors,myerror)
        F = findF(E,rho,P_real,d,A,B,S,T)
            myerror = maximum(norm(P_real[a][b][s][t] - generateP(rho,E,F,d,A,B,S,T)[a][b][s][t])  for t=1:T for s=1:S for b=1:B for a=1:A)
            push!(errors,myerror)
        rho = findrho(E,F,P_real,d,A,B,S,T)
    end

    errors

end

function randomP_increasing_d(num_iterations,d,d_max,epsilon,A,B,S,T)
    #The complete algorithm, determining the correlation and the entanglement
    # dimension
    # This version generates a P itself
    # In most cases, a solution is found for samples =1. However, there have
    # been instances where samples was higher.
    # Sometimes, here is a bug, in which the POVM's vanish, causing an error
    rho, E, F = generaterandomoperators(d,A,B,S,T)
    P_real = generateP(rho, E, F,d,A,B,S,T) # We want to approximate this P
    d=1
    errors = zeros(0)
    while d<d_max+1
        for sa in 1:samples
            errors = zeros(0)
            rho, E, F = generaterandomoperators(d,A,B,S,T) # We initialize with new rho,E,F
            for i in 1:num_iterations
                    myerror = maximum(norm(P_real[a][b][s][t] - generateP(rho,E,F,d,A,B,S,T)[a][b][s][t])  for t=1:T for s=1:S for b=1:B for a=1:A)
                    push!(errors,myerror)
                    println(errors)
                E = findE(F,rho,P_real,d,A,B,S,T)
                    myerror = maximum(norm(P_real[a][b][s][t] - generateP(rho,E,F,d,A,B,S,T)[a][b][s][t])  for t=1:T for s=1:S for b=1:B for a=1:A)
                    push!(errors,myerror)
                    println(errors)
                F = findF(E,rho,P_real,d,A,B,S,T)
                    myerror = maximum(norm(P_real[a][b][s][t] - generateP(rho,E,F,d,A,B,S,T)[a][b][s][t])  for t=1:T for s=1:S for b=1:B for a=1:A)
                    push!(errors,myerror)
                    println(errors)
                rho = findrho(E,F,P_real,d,A,B,S,T)
                if myerror < epsilon
                    return errors, d,sa,E,F,"Success"
                end
            end
        end
        d+=1
    end
    return errors,"Failed"
end

function givenP_increasing_d(P_real,num_iterations,d_max,epsilon,A,B,S,T)
    #The complete algorithm, determining the correlation and the entanglement
    # dimension
    # This version takes the correlation as input
    # In most cases, a solution is found for samples =1. However, there have
    # been instances where samples was higher.
    # There is a bug, in which the POVM's vanish, causing an error

    d=1
    errors = zeros(0)
    while d<d_max+1
        for sa in 1:samples
            errors = zeros(0)
            rho, E, F = generaterandomoperators(d,A,B,S,T) # We initialize with new rho,E,F
            for i in 1:num_iterations
                    myerror = maximum(norm(P_real[a][b][s][t] - generateP(rho,E,F,d,A,B,S,T)[a][b][s][t])  for t=1:T for s=1:S for b=1:B for a=1:A)
                    push!(errors,myerror)
                E = findE(F,rho,P_real,d,A,B,S,T)
                    myerror = maximum(norm(P_real[a][b][s][t] - generateP(rho,E,F,d,A,B,S,T)[a][b][s][t])  for t=1:T for s=1:S for b=1:B for a=1:A)
                    push!(errors,myerror)
                F = findF(E,rho,P_real,d,A,B,S,T)
                    myerror = maximum(norm(P_real[a][b][s][t] - generateP(rho,E,F,d,A,B,S,T)[a][b][s][t])  for t=1:T for s=1:S for b=1:B for a=1:A)
                    push!(errors,myerror)
                rho = findrho(E,F,P_real,d,A,B,S,T)
                if myerror < epsilon
                    return d,errors,E,F,rho,generateP(rho,E,F,d,A,B,S,T),sa,"Success"
                end
            end
        end
        d+=1
    end
    return errors,"Failed"
end

###############   Deterministic Correlations     ###############

# A=2
# B=2
# S=2
# T=2
# num_iterations = 20
# samples = 5
# d=1
# d_max = 4
# epsilon = 10^-8
#
# P_det = [
# [
# [[0,0],[0,0]],
# [[0,0],[0,0]]
# ],
# [
# [[0,0],[0,0]],
# [[0,0],[0,0]]
# ]
# ]
# P_det[1][1][1][1]=1
# P_det[1][2][1][2]=1
# P_det[2][1][2][1]=1
# P_det[2][2][2][2]=1
#
# P_used = P_det
#
# println("Start")
# println("====================================")
# println(givenP_increasing_d(P_used,num_iterations,d_max,epsilon,A,B,S,T))

#############      Private Randomness Correlations      #######################

# A=4
# B=4
# S=2
# T=2
# num_iterations = 10
# samples = 2
# d_max = 4
# epsilon = 10^-8
#
# function P_A_priv(a,s)
#     if a==1 && s==1
#         return 0.5
#     elseif a==2 && s==1
#         return 0.5
#     elseif a==3 && s==2
#         return 0.5
#     elseif a==4 && s==2
#         return 0.5
#     else
#         return 0
#     end
# end
#
# function P_B_priv(b,t)
#     if b==1 && t==1
#         return 1/4
#     elseif b==2 && t==1
#         return 3/4
#     elseif b==3 && t==2
#         return 0.1
#     elseif b==4 && t==2
#         return 0.9
#     else
#         return 0
#     end
# end
#
#
# P_pcorr = [[[[P_A_priv(a,s)*P_B_priv(b,t) for  t=1:T] for s=1:S] for b=1:B] for a=1:A]
#
# println("Start")
# println("====================================")
# println(givenP_increasing_d(P_pcorr,num_iterations,d_max,epsilon,A,B,S,T))

############Approximated  Shared Randomness Correlations    #################3

# A=2
# B=13
# S=2
# T=2
# num_iterations = 10
# samples = 2
# d_max = 4
# epsilon = 10^-3
#
# #the following empirical result is found using Python, with 10^6 samples.
# python = [[[[0.250679, 0.0], [0.041953, 0.0]], [[0.0, 0.0], [0.083058, 0.0]], [[0.249778, 0.014256], [0.041774, 0.027882]], [[0.0, 0.027541], [0.0, 0.0]], [[0.0, 0.041567], [0.0, 0.027583]], [[0.0, 0.05552], [0.0, 0.0]], [[0.0, 0.069334], [0.0, 0.027813]], [[0.0, 0.083663], [0.0, 0.0]], [[0.0, 0.069553], [0.0, 0.027837]], [[0.0, 0.055741], [0.0, 0.0]], [[0.0, 0.041969], [0.0, 0.02744]], [[0.0, 0.027966], [0.0, 0.0]], [[0.0, 0.014026], [0.0, 0.027884]]], [[[0.0, 0.0], [0.20881, 0.0]], [[0.500309, 0.0], [0.41673, 0.0]], [[0.0, 0.014123], [0.208567, 0.0]], [[0.0, 0.027882], [0.0, 0.055707]], [[0.0, 0.041501], [0.0, 0.055517]], [[0.0, 0.055349], [0.0, 0.11076]], [[0.0, 0.069946], [0.0, 0.111288]], [[0.0, 0.08325], [0.0, 0.166194]], [[0.0, 0.069919], [0.0, 0.111335]], [[0.0, 0.055449], [0.0, 0.110674]], [[0.0, 0.041436], [0.0, 0.055508]], [[0.0, 0.02827], [0.0, 0.055629]], [[0.0, 0.013967], [0.0, 0.0]]]]
#
# println("Start")
# println("====================================")
# println(givenP_increasing_d(python,num_iterations,d_max,epsilon,A,B,S,T))


################# Example Shared Randomness Correlations ####################
#
# A=2
# B=3
# S=2
# T=2
# num_iterations = 5
# samples = 1
# d_max = 7#S*(A-1)+T*(B-1)+S*T*(A-1)*(B-1)
# epsilon = 10^-4
#
# P_shar = [[[[0.0 for  t=1:T] for s=1:S] for b=1:B] for a=1:A]
#
# P_shar[1][1][1][1]=0.25
# P_shar[1][1][1][2]=0.5*5/6*5/6
# P_shar[1][1][2][1]=1/4*1/6
# P_shar[1][1][2][2]=5/36
# P_shar[1][2][1][1]=0.0
# P_shar[1][2][1][2]=0.5*2*1/6*5/6
# P_shar[1][2][2][1]=0.5*1/6
# P_shar[1][2][2][2]=0.0
# P_shar[1][3][1][1]=1/4
# P_shar[1][3][1][2]=0.5*1/36
# P_shar[1][3][2][1]=0.25*1/6
# P_shar[1][3][2][2]=1/36
# P_shar[2][1][1][1]=0.0
# P_shar[2][1][1][2]=0.5*5/6*5/6
# P_shar[2][1][2][1]=1/4*5/6
# P_shar[2][1][2][2]=5/6*4/6
# P_shar[2][2][1][1]=1/2
# P_shar[2][2][1][2]=0.5*2*1/6*5/6
# P_shar[2][2][2][1]=0.5*5/6
# P_shar[2][2][2][2]=2*1/6*5/6
# P_shar[2][3][1][1]=0.0
# P_shar[2][3][1][2]=0.5*1/36
# P_shar[2][3][2][1]=5/6*1/4
# P_shar[2][3][2][2]=0.0
#
# println("Start")
# println("====================================")
# println(givenP_increasing_d(P_shar,num_iterations,d_max,epsilon,A,B,S,T))


################ Given Shared Randomness correlation ##############

# A=2
# B=2
# S=2
# T=2
# num_iterations = 20
# samples = 2
# d_max = S*(A-1)+T*(B-1)+S*T*(A-1)*(B-1)
# epsilon = 10^-8

# P_shar = [[[[0.0 for  t=1:T] for s=1:S] for b=1:B] for a=1:A]
#
# P_shar[1][1][1][1]=0.25
# P_shar[1][1][1][2]=0.5*5/6*5/6
# P_shar[1][1][2][1]=1/4*1/6
# P_shar[1][1][2][2]=5/36
# P_shar[1][2][1][1]=1/4
# P_shar[1][2][1][2]=0.5*2*1/6*5/6 + 0.5*1/36
# P_shar[1][2][2][1]=0.5*1/6 + 0.25*1/6
# P_shar[1][2][2][2]=1/36
# P_shar[2][1][1][1]=0.0
# P_shar[2][1][1][2]=0.5*5/6*5/6
# P_shar[2][1][2][1]=1/4*5/6
# P_shar[2][1][2][2]=5/6*4/6
# P_shar[2][2][1][1]=1/2
# P_shar[2][2][1][2]=0.5*2*1/6*5/6 + 0.5*1/36
# P_shar[2][2][2][1]=0.5*5/6+5/6*1/4
# P_shar[2][2][2][2]=2*1/6*5/6


# println("Start")
# println("====================================")
# println(givenP_increasing_d(P_shar,num_iterations,d_max,epsilon,A,B,S,T))

################# Testing randomly generated classical correlations #########

# function randomshared(A,B,S,T,dets)
#     # Generates random classical correlations (C_loc(Gamma))
#     P_shar2 = [[[[[0.0 for  t=1:T] for s=1:S] for b=1:B] for a=1:A] for i=1:dets]
#     lambda=rand(Float64,(dets,1))
#     sumlambda = 0.0
#     for i=1:dets
#         sumlambda += lambda[i]
#     end
#     for i=1:dets
#         lambda[i]/=sumlambda
#     end
#     for i=1:dets
#         for s=1:S
#             for t= 1:T
#                 a=rand((1,A))
#                 b=rand((1,B))
#                 P_shar2[i][a][b][s][t]=1
#             end
#         end
#     end
#
#     result = [[[[0.0 for  t=1:T] for s=1:S] for b=1:B] for a=1:A]
#     for t=1:T
#         for s=1:S
#             for b=1:B
#                 for a=1:A
#                     res=0.0
#                     for i=1:dets
#                         res+=lambda[i]*P_shar2[i][a][b][s][t]
#                     end
#                     result[a][b][s][t]=res
#                 end
#             end
#         end
#     end
#
#     result = [[[[result[a][b][s][t]  for  t=1:T] for s=1:S] for b=1:B] for a=1:A]
# end
#
# A=2
# B=2
# S=2
# T=2
# num_iterations = 20
# samples = 2
# d_max = S*(A-1)+T*(B-1)+S*T*(A-1)*(B-1)
# epsilon = 10^-4
# dets=5
#
# P_shar3 = randomshared(A,B,S,T,dets)
#
# println("Start")
# println("====================================")
# println(givenP_increasing_d(P_shar3,num_iterations,d_max,epsilon,A,B,S,T))

############  Quantum correlations ###################

A=1
B=2
S=3
T=4
d=3
num_iterations = 20
samples = 5
d_max = d+1
epsilon = 10^-8

println("Start")
println("====================================")
println(randomP_increasing_d(num_iterations,d,d_max,epsilon,A,B,S,T))
#
