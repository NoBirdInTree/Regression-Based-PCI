options(rstudio.help.showDataPreview = FALSE)

# The following scripts:

# 1. Generate (Y,W) jointly from f(Y,W|A,Z,U) via a multinomial distribution 
# (which encodes the corresponding bivariate Bernouli distribution using 4 levels).
# 2. Evaluate the performance of 1.Two-Stage 2.Naive 3.Orale estimators respectively
# 3. Plot the graph and print the result table to report their performance
# Save it

# All parameters specified below match the ones in the Section A.8 of the Appendix.

library(ggplot2)
library(reshape2)
library(latex2exp)

multivariate_bernouli_generator = function(sample_space, prob){
  random_rv = runif(1)
  p = prob[1]
  for(i in seq(1,length(prob))){
    if(random_rv < p){
      return (sample_space[i])
    }
    else{
      p = p+prob[i+1]
    }
  }
}

# set parameters
beta_0 = -1.4
beta_a = 1.2
beta_u = -0.7
beta_w = 0.5
alpha_0 = -0.8
alpha_u = 0.5
scale= 0.3


N_Choice = c(250,500,1000,1500)
result_df = data.frame()
replication_time = 500
epoch = 300

# Load the previously generated large dataset of i.i.d. (A,Z,U)
load("azu.Rda")

for(N in N_Choice){
  print(paste0("The current N: ", N))
  set.seed(821743*N)
  idx = sample(1:dim(df)[1],N)
  # Draw (A,Z,U) from our previously generated large dataset (>5 million pairs)
  a = df$a[idx]
  z = df$z[idx]
  u = df$u[idx]
  
  # Calculate the probability of f(Y,W|A,Z,U)
  y_0_w_0 = 1
  y_1_w_0 = exp(beta_0 + beta_a*a + beta_u*u)
  y_0_w_1 = exp(alpha_0 + alpha_u*u)
  y_1_w_1 = exp(beta_0 + beta_a*a + beta_u*u) * exp(alpha_0 + alpha_u*u) * exp(beta_w)
  mat = cbind(y_0_w_0,y_1_w_0,y_0_w_1,y_1_w_1)
  mat = t(apply(mat, 1, function(x)(x/sum(x))))
  
  sample_space = list(c(0,0),c(1,0),c(0,1),c(1,1))
  for(r in c(1:replication_time)){
    if (r%%20== 0) print(paste0("The current R: ", r))
    # Generate (Y,W | U, A,Z) conditionally on (A,Z)
    y_w_joint = t(apply(mat, 1, function(x)(multivariate_bernouli_generator(sample_space,x))))
    y_w_joint = matrix(unlist(y_w_joint), ncol = 2, byrow = TRUE)
    y = y_w_joint[,1]
    w = y_w_joint[,2]
    tmp_df = data.frame(a,z,u,y,w)
    
    # Two-Stage Estimator
    
    # This vector saves the estimates for beta_a from bootstrap with epoch number prespecified
    hat_beta_usubstitute = rep(-999,epoch)
    for(e in seq(1:epoch)){
      indicies = sample(c(1:N),replace = TRUE)
      df_bootstrap = tmp_df[indicies,]
      m_sub_bootstrap = glm(w~a*z+y,family = binomial,data=df_bootstrap)
      # coef 1: the intercept
      # coef 2: the slope for A
      # coef 3: the slope for Z
      # coef 4: the slope for Y
      # coef 5: the slope for AZ
      df_bootstrap['u_sub'] =
        m_sub_bootstrap$coefficients[1] +
        m_sub_bootstrap$coefficients[2]*df_bootstrap$a +
        m_sub_bootstrap$coefficients[3]*df_bootstrap$z +
        m_sub_bootstrap$coefficients[4] +
        m_sub_bootstrap$coefficients[5]*df_bootstrap$a*df_bootstrap$z
      
      # 
      m_y_a_usubstitute_bootstrap = glm(y~a+u_sub+w,family = binomial, data=df_bootstrap)
      hat_beta_usubstitute[e] = m_y_a_usubstitute_bootstrap$coefficients[[2]]
    }
    
    # Obtain the estimate from Two-Stage estimator (S.E. is from bootstrap)
    m_sub = glm(w~a*z+y,family = binomial,data=tmp_df)
    summary(m_sub)
    # coef 1: the intercept
    # coef 2: the slope for A
    # coef 3: the slope for Z
    # coef 4: the slope for Y
    # coef 5: the slope for AZ
    tmp_df['u_sub'] =
      m_sub$coefficients[1] +
      m_sub$coefficients[2]*a +
      m_sub$coefficients[3]*z +
      m_sub$coefficients[4] +
      m_sub$coefficients[5]*a*z
    beta_a_estimate_star = glm(y~a+u_sub+w,family = binomial,data=tmp_df)$coefficients[[2]]
    
    # Oracle Estimator
    
    # obtain the estimate from Oracle estimator (S.E. from information matrix)
    tmp = glm(y~a+u+w,family = binomial,data=tmp_df)
    beta_a_estimate = tmp$coefficients[[2]]
    beta_a_estimate_se = sqrt(diag(vcov(tmp)))[2]
    
    # Naive Estimator
  
    # obtain the estimate from Naive estimator (S.E. from information matrix)
    tmp = glm(y~a+w,family = binomial,data=tmp_df)
    beta_a_estimate_ignore = tmp$coefficients[[2]]
    beta_a_estimate_ignore_se = sqrt(diag(vcov(tmp)))[2]
    
    # Save all these results

    BCI = quantile(hat_beta_usubstitute,probs=c(.025,.975))
    
    tmp = data.frame(N=N,
                     rep_time=r,
                     beta_a_estimate = beta_a_estimate, # the estimate from Oracle Estimator
                     beta_a_estimate_se = beta_a_estimate_se, # the S.E. of the estimate from Oracle Estimator
                     beta_a_estimate_star = beta_a_estimate_star, # the estimate from Two-Stage Estimator
                     BSE = sd(hat_beta_usubstitute), # the bootstrap S.E. of Two-Stage Estimator
                     BM = mean(hat_beta_usubstitute), # the bootstrap mean of Two-Stage Estimator
                     beta_a_estimate_ignore = beta_a_estimate_ignore, # the estimate from Naive Estimator
                     beta_a_estimate_ignore_se = beta_a_estimate_se, # the S.E. of the estimate from Naive Estimator
                     BCI_l = BCI[1], # the 2.5% quantile of the estimate from boostrap for Two-Stage Estimator
                     BCI_u = BCI[2]. # the 97.5% quantile of the estimate from boostrap for Two-Stage Estimator
                     )
    result_df = rbind(result_df,tmp)
  }
}

# Save the simulation result
save(result_df,file="performace_evaluation.Rda")

# Load the simulation result
load("performace_evaluation.Rda")
colnames(result_df)[c(8,9)] = c('beta_a_estimate_ignore','beta_a_estimate_ignore_se')
result_df$N = as.factor(result_df$N)
result_df$Bias_beta_a_ignore = result_df$beta_a_estimate_ignore - beta_a
result_df$Bias_beta_a_star = result_df$beta_a_estimate_star - beta_a
result_df$Bias_beta_a = result_df$beta_a_estimate - beta_a
formatted_result_df = result_df[,c('N','Bias_beta_a_ignore','Bias_beta_a_star','Bias_beta_a')]
colnames(formatted_result_df) =c('N','Naive','Two-Stage','Oracle')
formatted_result_df = melt(formatted_result_df, id='N')

# Evaluate Two-Stage estimator
lower = result_df$BCI_l
upper = result_df$BCI_u
result_df$Covered_beta_a_estimate_star = (lower < beta_a & beta_a < upper)
Bias = rep(-999,length(N_Choice))
BSE = rep(-999,length(N_Choice))
ESE = rep(-999,length(N_Choice))
Coverage_BSE = rep(-999,length(N_Choice))
Replication = rep(replication_time,length(N_Choice))
Sanity_Check = rep(-999,length(N_Choice))
v = 0
for(N in N_Choice){
  v = v+1
  tmp = result_df[result_df$N==N,]
  Bias[v] = mean(tmp$Bias_beta_a_star)
  ESE[v] = sd(tmp$beta_a_estimate_star)
  BSE[v] = mean(tmp$BSE)
  Coverage_BSE[v] = sum(tmp$Covered_beta_a_estimate_star)/nrow(tmp)
  Sanity_Check[v] = sum(tmp$beta_a_estimate_star - tmp$BM)/nrow(tmp)
}
quality_df_beta_a_estimate_star = data.frame(N_Choice,Bias,ESE,BSE,Coverage_BSE,Sanity_Check,Replication)
quality_df_beta_a_estimate_star
present = quality_df_beta_a_estimate_star
present = present[,c("N_Choice","Bias","ESE","BSE","Coverage_BSE","Replication")]
library(xtable)
print(xtable(present,digits=c(0,0,2,2,2,2,0)),include.rownames=FALSE)

# Graphical Evaluation
colnames(formatted_result_df)[2] = 'Estimator'
p<-ggplot(formatted_result_df, aes(x=N, y=value, fill=Estimator)) +
  geom_boxplot(position = "dodge2") +
  ylim(-2.5,2.5) +
  xlim("250", "500", "1000", "1500","-") + 
  ylab(TeX(r'($Bias$)')) + 
  theme(axis.title.y = element_text(angle = 0)) +
  xlab("N") + 
  geom_hline(yintercept = 0, color="red", linetype="dashed") +
  scale_fill_manual(
    values=c('orchid2','#FF9999','darkslategray1')
  ) +
  theme_bw() +
  theme(legend.key.size = unit(1.1, 'cm')) +
  theme(legend.position = c(0.90, 0.5)) + 
  theme(legend.text.align = 0) +
  theme(legend.title = element_text(size=13)) +
  theme(legend.text = element_text(size=10)) +
  coord_cartesian(clip = 'off')
p

# ------

# Evaluate Naive estimator
# alpha = 0.05
# lower = result_df$beta_a_estimate - qnorm(1-alpha/2)*result_df$beta_a_estimate_se
# upper = result_df$beta_a_estimate + qnorm(1-alpha/2)*result_df$beta_a_estimate_se
# result_df$Covered_beta_a_estimate = (lower < beta_a & beta_a < upper)
# 
# Bias = rep(-999,length(N_Choice))
# ESE = rep(-999,length(N_Choice))
# SE_Wald = rep(-999,length(N_Choice))
# Coverage_Wald = rep(-999,length(N_Choice))
# Replication = rep(replication_time,length(N_Choice))
# v = 0
# for(N in N_Choice){
#   v = v+1
#   tmp = result_df[result_df$N==N,]
#   Bias[v] = mean(tmp)
#   ESE[v] = sd(tmp$beta_a_estimate)
#   SE_Wald[v] = mean(tmp$beta_a_estimate_se)
#   Coverage_Wald[v] = sum(tmp$Covered_beta_a_estimate)/nrow(tmp)
# }
# quality_df_beta_a_estimate = data.frame(N_Choice,Bias,ESE,SE_Wald,Coverage_Wald,Replication)
# quality_df_beta_a_estimate
# 
# library(xtable)
# xtable(quality_df_beta_a_estimate,digits=c(0,0,6,6,6,6,0))