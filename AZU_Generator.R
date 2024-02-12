options(rstudio.help.showDataPreview = FALSE)

# The following scripts:

# 1. Generate A and Z independently from N(0,0.5)
# 2. Generate U from f(U|A,Z) via the accept-reject sampling algorithm.
# Repeat step 1 and 2 many times and generate a large set of data pair (A,Z,U)
# Save it

# All parameters specified below match the ones in the Section A.8 of the Appendix.

set.seed(666)

# m_y0 is E[U|A,Z,Y=0,W=0]
m_y0 = function(A,Z){return (-0.4 + 0.8*A + 1.2*Z - 1*A*Z)}

# set parameters
beta_0 = -1.4
beta_a = 1.2
beta_u = -0.7
beta_w = 0.5
alpha_0 = -0.8
alpha_u = 0.5
scale= 0.3

# The trail size
N = 20000000
df = data.frame()

for(i in (1:8)){
  # Generate continuous A,Z
  a = rnorm(N,0,0.5)
  z = rnorm(N,0,0.5)

  # ACCEPT-REJECT ALGORITHM. Generate U from f(U|A,Z).
  trail_size = N
  u_lower_bound = -7
  u_upper_bound = 4
  u = runif(N,u_lower_bound,u_upper_bound)
  eps = u-m_y0(a,z)
  merged_df = data.frame(a_sampled = a, z_sampled = z, u_sampled = u, eps_sampled=eps)
  merged_df$density = dlogis(eps,scale=scale)*
    ( 1+exp(beta_0+beta_a*a+beta_u*u) + (1+exp(beta_0+beta_a*a+beta_u*u))*exp(alpha_0+alpha_u*u+beta_w))/
    ( (1+exp(beta_0+beta_a*a)) * (1+exp(alpha_0)) )
  
  # Normalize the density
  merged_df$density = merged_df$density/max(merged_df$density)

  # Only accept the pair when runif gives a lower value
  merged_df$accpted = runif(trail_size,0,1) < merged_df$density
  df1 = merged_df[merged_df$accpted==1,colnames(merged_df)[1:3]]
  
  df = rbind(df,df1)
}

save(df, file="azu.Rda")

# summary(m_y0(a,z))
# summary(merged_df$u_sampled[merged_df$accpted])
# hist(merged_df$u_sampled[merged_df$accpted], breaks = 200)