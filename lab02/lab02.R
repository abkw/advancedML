library(HMM)
RNGversion(vstr = 3.6)
set.seed(12345)

#Question 1: Build a hidden Markov model (HMM) for the scenario described above.

# Initialise HMM
states_vec = c(1:10)
symbols_vec = LETTERS[seq( from = 1, to = 10 )]
start_prob = rep(0.1,10)

#The rebot can stay in the current sector or move to the next sector with equal probability.
trans_probs=matrix(data = 0, nrow = length(states_vec), ncol = length(states_vec))
diag(trans_probs) = 0.5
for (i in 1:dim(trans_probs)[2]){
  if (i == dim(trans_probs)[2]){
    trans_probs[dim(trans_probs)[1],1] = 0.5
  }
  else
  trans_probs[i,i+1] = 0.5
}
trans_probs

#the device will report that the robot is in the sectors [i−2, i+2] with equal probability
#initializing the emission probabilities, NAs are used in order to forbid the transition between the states.

accuracy = matrix(c(0.2, 0.2, 0.2, NA, NA, NA, NA, NA, 0.2, 0.2,
                    0.2, 0.2, 0.2, 0.2, NA, NA, NA, NA, NA, 0.2,
                    0.2, 0.2, 0.2, 0.2, 0.2, NA, NA, NA, NA, NA,
                    NA, 0.2, 0.2, 0.2, 0.2, 0.2, NA, NA, NA, NA,
                    NA, NA, 0.2, 0.2, 0.2, 0.2, 0.2, NA, NA, NA,
                    NA, NA, NA, 0.2, 0.2, 0.2, 0.2, 0.2, NA, NA,
                    NA, NA, NA, NA, 0.2, 0.2, 0.2, 0.2, 0.2, NA,
                    NA, NA, NA, NA, NA, 0.2, 0.2, 0.2, 0.2, 0.2,
                    0.2, NA, NA, NA, NA, NA, 0.2, 0.2, 0.2, 0.2,
                    0.2, 0.2, NA, NA, NA, NA, NA, 0.2, 0.2, 0.2), nrow = 10, ncol = 10)


hmm = initHMM(states_vec, symbols_vec, start_prob, trans_probs, accuracy)
cat('The states of HMM are: \n',hmm$States)
cat('The symbols of HMM are: \n',hmm$Symbols)
cat('The transition probabilities of HMM are: \n')
print(hmm$transProbs)

#Question 2: Simulate the HMM for 100 time steps.
simulate_hmm = function(no_of_simulations){
  
  
  #Question 2: Simulate the HMM for 100 time steps.
  length = no_of_simulations
  
  cat('>>>>>>>>>>>>>>>>>>>>>>>>>>>>Simulating for ', length , 'time steps >>>>>>>>>>>>>>>>>>>>>>>>>>>> \n')
  simulated_hmm = simHMM(hmm, length)
  
  cat('The simulated states are: \n')
  #print(simulated_hmm$states)
  cat('The simulated observations are: \n')
  #print(simulated_hmm$observation)
  
  
  #Question 3: Discard the hidden states from the sample obtained above. Use the remaining obser- vations to compute the filtered and smoothed probability distributions for each of the 100 time points. Compute also the most probable path.
  
  
  #The forward-function computes the forward probabilities. The forward probability for state X up to observation at time k is defined as the probability of observing the sequence of observations e_1, ... ,e_k and that the state at time k is X
  filtered = forward(hmm, simulated_hmm$observation)
  
  
  cat('The filtered probabilities are: \n')
  #print(filtered)
  
  #The posterior function computes the posterior probabilities of being in state X at time k for a given sequence of observations and a given Hidden Markov Model.
  smoothed = posterior(hmm, simulated_hmm$observation)
  cat('The posterior probabilities of states given the observations are: \n')
  #print(smoothed)
  
  #The values need to be normalized.
  
  #Calculating the most probable path: 
  probable_path = viterbi(hmm,simulated_hmm$observation)
  #probable_path[1:100]
  
  #Plotting the most probable path: 
  plot(c(1:length), 
       probable_path, 
       main = paste("The most probable path with sample size",length),
       xlab = "Time Points", 
       ylab = "states", 
       type = 'l', 
       col = 'blue',
       lwd = 3)
  mtext('The robot movement between different time points')
  abline(h = pretty(probable_path, 10),
         v = pretty(c(1:length), length), 
         col = "black")
  
  #Question 4: Compute the accuracy of the filtered and smoothed probability distributions, and of the most probable path. That is, compute the percentage of the true hidden states that are guessed by each method.
  
  #Normalizing the filtered states: 
  filtered = exp(filtered)
  normalized_filter = prop.table(filtered,2)
  
  #The accuracy of the filtered probabilities############################################################
  
  #Taking the maximum of the normalized filter of each time point 
  most_probable_filter = apply(normalized_filter, 2, which.max)
  most_probable_filter = as.vector(most_probable_filter)
  #Taking the maximum of the actual simulated states
  actual_states = simulated_hmm$states
  
  #Comparing the actual states with the filtered probable states:
  
  filter_accuracy = which(most_probable_filter == actual_states)
  cat('The states that are equal to the filter are: \n')
  print(length(filter_accuracy))
  
  cat('The accuracy of filtered probabilities is: \n')
  cat(length(filter_accuracy)/length(actual_states), '\n')
  
  #The accuracy of the smoothed probabilities############################################################
  most_probable_smoothed = apply(smoothed, 2, which.max)
  most_probable_smoothed = as.vector(most_probable_smoothed)
  
  #Comparing the actual states with the smoothed probable states:
  
  smooth_accuracy = which(most_probable_smoothed == actual_states)
  cat('The states that are equal to the smoothed are: \n')
  print(length(smooth_accuracy))
  cat('The accuracy of smoothed probabilities is: \n')
  cat(length(smooth_accuracy)/length(actual_states),'\n')
  
  #Calculating the probability of the most probable path ##################################################
  plot(c(1:length), 
       probable_path, 
       main = paste("The actual probable path vs the simulated sample",length),
       xlab = "Time Points", 
       ylab = "states", 
       type = 'l', 
       col = 'blue',
       lwd = 1)
  lines(c(1:length),simulated_hmm$states, col = 'red')
  legend('topright',
         legend = c("simulated", "actual"),
         col = c('blue','red'), 
         lty=1:5, 
         cex=0.6)
  
  cat('The accuracy of the most probable path is: \n')
  acuracy_most_pp = length(which(probable_path==simulated_hmm$states))/length(probable_path)
  cat(acuracy_most_pp, '\n')
  filtered_entropy = apply(normalized_filter, 2, entropy.empirical)
  plot(1:length,filtered_entropy, 
       type = 'l',
       col = 'red',
       xlab = 'Time Step',
       ylab = 'Empirical Entropy',
       main = paste('The Empirical Entropy for ',length, 'Sample Size'))
#
}

for(sample in c(100,200,300,400, 1000)){
  simulate_hmm(sample)
}

## Question 5:
#Repeat the previous exercise with different simulated samples. In general, the smoothed distributions should be more accurate than the filtered distributions. Why ? In general, the smoothed distributions should be more accurate than the most probable paths, too. Why ?

#The smoothed distribution is more accurate than the filtered one because it takes all the observations up
#to time T for zt, whereas the filtered distribution takes the observation up to *t* for the state zt, so, by
#considerting all the observation and not subset of the observation the smoothed distribution is more accurate
#the the filtered distribution.

#In the most probable path we take the state that maximizes the expression up to time z(t+1). For instance, to 
#maximize over z_0 we have to maximize for all the possible values of z_1 because we don't know the value of z_1
#thus, in most probable path we don't take all the data into consideration (only current and next states and observations), and that's why smoothed distribution
#is more accurate than the most probable path.

# Question 6:  

#From the graphs, the entropy appears to be random no matter the sample size, and that's because
#it depends only on the previous observation, so, that means the increase in the number of observations
#doesn't give us a better knowledge of where the robot is.


#Question 7:
#Consider any of the samples above of length 100. Compute the probabilities of the
#hidden states for the time step 101.


## References:
#https://cran.r-project.org/web/packages/HMM/HMM.pdf.

