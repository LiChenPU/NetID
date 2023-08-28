################

# Run this file for a toy linear programing problem
# solved by cplexAPI

################
library(cplexAPI)

# In the following, an example MIP will beW created and solved: 
# maximize
# z = x1 + 2x2 + 3x3 + x4
# subject to
# −x1 + x2 + x3 + 10x4 ≤ 20
# x1 − 3x2 + x3 ≤ 30
# x2 − 3.5x4 = 0
# With all variables being non-negative, x1 ≤ 40 and x4 ∈ {2, 3, 4} (x4 is integer).

# Open a IBM ILOG CPLEX environment.
env <- openEnvCPLEX()
# Create a problem object.
prob <- initProbCPLEX(env, pname = "example")
# Prepare data structures for the problem object. Number of columns, rows and non-zero elements.
nc <- 4
nr <- 3
nz <- 9
# Objective function.
obj <- c(1.0, 2.0, 3.0, 1.0)
# Right hand side.
rhs <- c(20.0, 30.0, 0.0)
# Sense of the right hand side.
sense <- c("L", "L", "E")
# Vatiable types.
ctype <- c("C", "C", "C", "C")
# Variable lower bounds.
lb <- c(0.0, 0.0, 0.0, 2.0)
# Variable upper bounds.
ub <- c(40.0, CPX_INFBOUND, CPX_INFBOUND, 3.0)
# The constraint matrix is passed in column major order format. Be careful here: all
# indices start with 0! Begin indices of rows.
beg <- c(0, 2, 5, 7)
# Number of non-zero elements per row.
cnt <- c(2, 3, 2, 2)
# Column indices.
ind <- c(0, 1, 0, 1, 2, 0, 1, 0, 2)
# Non-zero elements.
val <- c(-1.0, 1.0, 1.0, -3.0, 1.0, 1.0, 1.0, 10.0, -3.5)
# Load problem data.
copyLpCPLEX(env, prob, nc, nr, CPX_MAX, obj, rhs, sense,
              beg, cnt, ind, val, lb, ub)
# Set Variable types.
copyColTypeCPLEX(env, prob, ctype)
# Solve the problem using MIP.
mipoptCPLEX(env, prob)
# Retrieve solution after optimization.
solutionCPLEX(env, prob)
# All counts start with 0
addRowsCPLEX(env, prob, 0, 1, 1, 0, 0, 1,
             rhs = 30, sense = "E",
             cnames = NULL, rnames = NULL)
delRowsCPLEX(env, prob, 3,3)

delProbCPLEX(env, prob)
closeEnvCPLEX(env)

# checkAddRowsCPLEX(env, lp, ncols, nrows, nnz, matbeg, matind, matval,
#                   rhs = NULL, sense = NULL,
#                   cnames = NULL, rnames = NULL)
# delRowsCPLEX(env, lp, begin, end)

# https://www.ibm.com/support/knowledgecenter/en/SSSA5P_12.5.1/ilog.odms.cplex.help/refcallablelibrary/html/functions/CPXaddrows.html

# writeProbCPLEX(env, prob, "prob.lp")
# lp <- initProbCPLEX(env)
# readCopyProbCPLEX(env, lp, "prob.lp")

# sense[i]	= 'L'	<= constraint
# sense[i]	= 'E'	= constraint
# sense[i]	= 'G'	>= constraint
# sense[i]	= 'R'	ranged constraint

# CPX_CONTINUOUS	'C'	continuous variable j
# CPX_BINARY	'B'	binary variable j
# CPX_INTEGER	'I'	general integer variable j
# CPX_SEMICONT	'S'	semi-continuous variable j
# CPX_SEMIINT	'N'	semi-integer variable j

