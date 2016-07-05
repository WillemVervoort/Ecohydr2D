# growth rate functions
# Ursino 2007, Table 1

n = seq(0,1,by=0.01)
w = seq(0,1,by=0.01)

m = 0.45
r = 0.36

# hypothesis 1 (linear)
foo = function(w,n,m,r) {- m*(1-(r*w/m)*n)}

plot(w,foo(w,0.5,m,r))

# hypothesis 2
g = 5 # acts as Emax
k1 = 0.075 # curvature, smaller is more curve
m = 2.74 # scales on y axis
foo = function(w,m,k1,g) {g*(w/(k1+w))-m}
# not a function of n (biomass itself, no feedback)
plot(w,foo(w,m,k1,g))

# hypothesis 3a
b = 2 # Emax (similar)
k_p_n = 10 # scales y-axis

foo = function(w,n,b,k_p_n) {b*(1-n/(k_p_n*w))}
par(mfrow=c(2,1))
plot(w,foo(w,0.5,b,k_p_n))
plot(n,foo(0.5,n,b,k_p_n))
par(mfrow=c(1,1))

# hypothesis 3b (both linear)
b_p = 2 # 
k_n = 5 # scale

foo = function(w,n,b_p,k_n) {b_p*w*(1-n/(k_n))}
par(mfrow=c(2,1))
plot(w,foo(w,0.5,b_p,k_n))
plot(n,foo(0.5,n,b_p,k_n))
par(mfrow=c(1,1))


# from slides
g_max = 5 # 
k1 = 0.075 # curvature

# no negative feedback with growth
foo = function(w,n,g_max,k1) {g_max*n*(w/(w+k1))}
par(mfrow=c(2,1))
plot(w,foo(w,0.5,g_max,k1))
plot(n,foo(0.5,n,g_max,k1))
par(mfrow=c(1,1))

