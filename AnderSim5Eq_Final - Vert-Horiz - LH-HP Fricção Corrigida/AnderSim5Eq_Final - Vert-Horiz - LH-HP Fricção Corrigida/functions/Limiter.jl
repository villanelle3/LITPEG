function Limiter(r)
  #limiter = (r+abs(r))/(r+1) # Van Leer
  #limiter = maximum([0.0, minimum([2*r, (1+r)/2, 2])])  #VANLEER 77 (MUSCL)
  #limiter = maximum([0.0, minimum([2*r, 1]),minimum([r, 2])])   #Superbee
  #limiter = maximum([0.0, minimum([1,r])]) # Minmod
  #limiter = 1 # SOU
 limiter = 0 # FOU
  return limiter
end
