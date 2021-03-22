library(markovchain)
library(Matrix)

#States space
phoneStates <- c("Off","On","Lock","Unlock")
phoneStates

#Forming transition matrix
phoneMatrix <- matrix(data = c(0.51/100,33.04/100,59.4/100,7.05/100,45.32/100,2.03/100,0,52.65/100,2.83/100,95.64/100,0,1.53/100,80.50/100,13.58/100,0,5.92/100),
                      byrow = TRUE, nrow = 4,dimnames = list(phoneStates, phoneStates))

mcphone <- new("markovchain", states = phoneStates, byrow = TRUE,
               transitionMatrix = phoneMatrix, name = "Phone")

mcphone

#Markov chain
is(mcphone,"markovchain")

#Transition plot
plot(mcphone,main="Phone Markov Chain")

#Initial conditions: the state it starts in time, t=0
#Start from Off
s0= matrix(c(1,0,0,0),nrow = 1,byrow = T)
s1=s0*mcphone
s1

#Start from On
s0= matrix(c(0,1,0,0),nrow = 1,byrow = T)
s1=s0*mcphone
s1

#Start from Lock
s0= matrix(c(0,0,1,0),nrow = 1,byrow = T)
s1=s0*mcphone
s1

#Start from Unlock
s0= matrix(c(0,0,0,1),nrow = 1,byrow = T)
s1=s0*mcphone
s1

#Steady state
phoneMatrix

fvalue=function(phoneMatrix,n)
{
  pn=phoneMatrix
  
  for(i in 2:n)
  {
    pd=pn
    pn=phoneMatrix%*%pd
    print(i)
    round(pn,4)
    print(pn)
  }
  
}

fvalue(phoneMatrix,20)


#After raising transition probability matrix to a high enough power, which is 20, we can see that the steady-state distribution is 𝛱 = [ 0.31 0.31 0.18 0.2 ].


#Number of visits in 3 steps
v=diag(4)+phoneMatrix+phoneMatrix%*%phoneMatrix+phoneMatrix%*%phoneMatrix%*%phoneMatrix
v

#Suppose the screen is off at initial, then the expected number of total visits to Off, On, Lock and Unlock are 1.6462560, 1.030522, 0.7296850 and 0.5935370 respectively.

#Suppose the screen is On at initial, then the expected number of total visits to Off, On, Lock and Unlock are 1.0582482, 1.657772, 0.5277950 and 0.7561850 respectively.

#Suppose the screen is Lock at initial, then the expected number of total visits to Off, On, Lock and Unlock are 0.8986178, 1.220047, 1.2816756 and 0.5996592 respectively.

#Suppose the screen is Unlock at initial, then the expected number of total visits to Off, On, Lock and Unlock are 1.1639109, 0.930837, 0.5454738 and 1.3597783 respectively.

#Reachability
p<-phoneMatrix
p1=p
p2=p1%*%p1
p3=p2%*%p1
p4=p2%*%p2
p5=p3%*%p2
p6=p3%*%p3
p7=p5%*%p2
p8=p4%*%p4
p9=p5%*%p4
p10=p5%*%p5
p11=p10%*%p1
p12=p10%*%p2
p13=p10%*%p3
p14=p10%*%p4
p15=p10%*%p5
p16=p10%*%p6
p17=p10%*%p7
p18=p10%*%p8
p19=p10%*%p9
p20=p16%*%p4
p20
round(p20,2)
I=diag(4)
sum=I+p1+p2+p3+p4+p5+p6+p7+p8+p9+p10+p11+p12+p13+p14+p15+p16+p17+p18+p19
round(sum,3)

is.accessible(mcphone)

#All states are reachable from On, Off, Lock and Unlock.

#Summary
summary(mcphone)

#In this Markov model, all states communicate with each other and they form a closed class = {Off, On, Lock, Unlock}. It is irreducible as the closed class is the set of states in its state space. 

#Absorbing state is once the chain reaches the state, it stays there forever and impossible to leave. There is no absorbing state in this Markov model.

transientStates(mcphone)

recurrentClasses(mcphone)

#Recurrent class is if a state is recurrent, then all states in the class must be recurrent. So, the recurrent class in this Markov model is {Off, On, Lock, Unlock}.

#Periodity
period(mcphone)

#A state is said to have a period d when the chain can only revisit it after a multiple of d steps. Since d=1, the state is aperiodic.

#Limiting distribution
#Long term visiting rates
e=matrix(c(1,1,1,1),nrow = 1, byrow = T)
I=diag(4)
E=matrix(1,4,4)
e%*%solve(I+E-phoneMatrix)

# simulate discrete Markov chains
run.mc.sim <- function( P,   # probability transition matrix
                        num.iters=50, 
                        num.chains=150 )
{
  
  # number of possible states
  num.states <- nrow(P)
  
  # states X_t for all chains
  states     <- matrix(NA, ncol=num.chains, nrow=num.iters)
  
  # probability vectors pi^n through time
  all_probs  <- matrix(NA, nrow=num.iters, ncol=num.states)
  
  # forces chains to start in state 1
  pi_0      <- c(1, rep(0, num.states-1))
  
  # initialize variables for first state 
  P_n           <- P
  all_probs[1,] <- pi_0
  states[1,]    <- 1
  
  for(t in 2:num.iters) {
    
    # pi^n for this iteration
    pi_n           <- pi_0 %*% P_n
    all_probs[t,]  <- pi_n
    
    for(chain_num in seq_len(num.chains)) {
      # probability vector to simulating next state 
      p                     <- P[ states[t-1,chain_num], ]
      states[t,chain_num]   <- which(rmultinom(1, 1, p) == 1)
    }
    
    # update probability transition matrix
    P_n           <- P_n %*% P
  }
  return(list(all.probs=all_probs, states=states))
}

P<-phoneMatrix
sim1 <- run.mc.sim(P)
states <- sim1[[2]]
matplot(states[,1:5], type='l', lty=1, col=1:5, ylim=c(0,4), ylab='state', xlab='time')
abline(h=1, lty=4)
abline(h=4, lty=4)

all.probs <- sim1[[1]]
matplot(all.probs, type='l', col=1:4, lty=1, ylab='probability', xlab='time')
legend('topright', c('Off', 'On', 'Lock','Unlock'), lty=1, col=1:4)

#The probabilities quickly converge to the stationary distribution, ℼ as we found above.

#First reaching time
meanFirstPassageTime(mcphone)

#Expected return time
#for Off
Eoff=1/0.3077797
Eoff
#for On
Eon=1/0.3099211
Eon
#for Lock
Elock=1/0.1828211
Elock
#for Unlock
Eunlock=1/0.1994782
Eunlock
