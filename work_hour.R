#BH procedure
BH <- function(p_values, alpha) 
{
  ##inputs
  #p_values : vector of p values
  #alpha : target level
  
  ##output
  #rejection threshold
  
  p_sorted <- sort(p_values)
  m <- length(p_values)
  if(min(p_sorted*m/(1:m)) <= alpha)
  {
    threshold <- p_sorted[max(which(p_sorted <= (1:m) / m * alpha))]
  }
  else
  {
    threshold <- 0
  }
  return(threshold)
  
}

##batch conformal p-value
batch_conformal <- function(S_cal, S_test, eta)
{
  n <- length(S_cal)
  m <- length(S_test)
  probs <- sapply(1:(n+1), function(i) choose(i+eta-2,eta-1)*choose(n+m-i-eta+1,m-eta)/choose(n+m,m))
  t <- sum(S_cal < sort(S_test)[eta])
  return(sum(probs[(t+1):(n+1)]))
}

##------------------------------------------------------------------------------------------------

#load data
load("cpshours.RData")

set.seed(12)

#construct age groups
wdat$AGE_GROUP <- cut(
  wdat$AGE,
  breaks = c(24, 39, 49, 59, Inf),
  labels = c("25–39", "40–49", "50–59", "60+"),
  right = TRUE
)

#construct groups based on age, sex, race, and education
group_vars <- c("AGE_GROUP", "SEX", "RACE", "EDUCNUM")
subgroup_list <- split(wdat, wdat[group_vars])


#set the reference group
print(which.max(sapply(subgroup_list,function(subdata) dim(subdata)[1]))) #73
ref_group_index <- 73

comp_group_index <- setdiff(1:length(subgroup_list),c(ref_group_index,which(sapply(subgroup_list,function(data) dim(data)[1])<=4)))

#compute reference scores
score_ref <- subgroup_list[[ref_group_index]]$HOURS + runif(dim(subgroup_list[[ref_group_index]])[1],-0.01,0.01)
score_ref <- sample(score_ref,100,replace=F) #experiment under small-sample setting

comp_scores <- list()

#compute scores

for(k in 1:length(comp_group_index))
{
  sub_k <- subgroup_list[[comp_group_index[k]]]
  score_k <- sub_k$HOURS + runif(dim(sub_k)[1],-0.01,0.01)
  if(length(score_k) > 50)
  {
    score_k <- sample(score_k,50,replace=F)
  }
  comp_scores[[k]] <- score_k
}

#compute batch conformal p-values
p_values_1 <- sapply(comp_scores,function(v) batch_conformal(score_ref, v, ceiling(length(v)*0.25)))
p_values_2 <- sapply(comp_scores,function(v) batch_conformal(score_ref, v, ceiling(length(v)*0.5)))
p_values_3 <- sapply(comp_scores,function(v) batch_conformal(score_ref, v, floor(length(v)*0.75)))

rej_1_1 <- which(p_values_1 <= BH(p_values_1,0.01))
rej_1_2 <- which(p_values_1 <= BH(p_values_1,0.02))
rej_2_1 <- which(p_values_2 <= BH(p_values_2,0.01))
rej_2_2 <- which(p_values_2 <= BH(p_values_2,0.02))
rej_3_1 <- which(p_values_3 <= BH(p_values_3,0.01))
rej_3_2 <- which(p_values_3 <= BH(p_values_3,0.02))


#results
t(sapply(c(comp_group_index[rej_1_1]), function(k) subgroup_list[[k]][1,c(2,3,5,8)]))
t(sapply(c(comp_group_index[rej_1_2]), function(k) subgroup_list[[k]][1,c(2,3,5,8)]))
t(sapply(c(comp_group_index[rej_2_1]), function(k) subgroup_list[[k]][1,c(2,3,5,8)]))
t(sapply(c(comp_group_index[rej_2_2]), function(k) subgroup_list[[k]][1,c(2,3,5,8)]))
t(sapply(c(comp_group_index[rej_3_1]), function(k) subgroup_list[[k]][1,c(2,3,5,8)]))
t(sapply(c(comp_group_index[rej_3_2]), function(k) subgroup_list[[k]][1,c(2,3,5,8)]))


