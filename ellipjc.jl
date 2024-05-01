using Polynomials

function ellipkkp(L)
  if L > 10
    K = pi/2;
    Kp = pi*L + log(4);
    return K,Kp
  end
  m = exp(-2*pi*L);
  a0 = 1;
  b0 = sqrt(1-m);
  s0 = m;
  i1 = 0; 
  mm = 1;
  a1 = (a0+b0)/2;
  while mm > eps() # eps funziona?
    a1 = (a0+b0)/2;
    b1 = sqrt(a0.*b0);
    c1 = (a0-b0)/2;
    i1 = i1 + 1;
    w1 = 2^i1*c1.^2;
    mm = max(max(w1));
    s0 = s0 + w1;
    a0 = a1;
    b0 = b1;
  end
  K = pi./(2*a1);
  # im = find(m==1);
  # if isempty(im)
  #   K[im] = K[im]*inf;
  # end
    a0 = 1;
    b0 = sqrt(m);
    s0 = 1-m;
    i1 = 0;
    mm = 1;
    while mm > eps() # eps() funziona?
      a1 = (a0+b0)/2;
      b1 = sqrt(a0.*b0);
      c1 = (a0-b0)/2;
      i1 = i1 + 1;
      w1 = 2^i1*c1.^2;
      mm = max(max(w1));
      s0 = s0 + w1;
      a0 = a1;
      b0 = b1;
    end
    Kp = pi./(2*a1);
    # im = find(m==0);
    # if isempty(im)
    #   Kp[im] = Kp[im]*inf;
    # end
  return K,Kp  
end

function ellipjc(u,L,flag)
  if flag == 1
    high = zeros(size(u));
    m = L;
  else
    K, Kp = ellipkkp(L);
    high = imag(u) .> Kp/2;
    u[high] = im*Kp .- u[high];
    m = exp(-2*pi*L);
  end
if m < 4*eps()
  sinu = sin.(u);
  cosu = cos.(u);
  sn = sinu + m/4*(sinu .* cosu - u).*cosu;
  cn = cosu + m/4*(-sinu .* cosu + u).*sinu;
  dn = 1 .+ m/4*(cosu.^2 .- sinu.^2 .- 1);
else
  if m > 1e-3
    kappa = (1-sqrt(1-m))/(1+sqrt(1-m));
  else
    pol = Polynomial([0,1,2,5,14,42,132]);
    kappa = pol(m/4);
  end
  mu = kappa^2;
  println(mu);
  v = u/(1+kappa);
  sn1,cn1,dn1 = ellipjc(v,mu,1);
  denom = (1 .+ kappa*sn1.^2);
  sn = (1 .+ kappa)*sn1 ./ denom;
  cn = cn1.*dn1 ./ denom;
  dn = (1 .- kappa .* sn1.^2) ./ denom;
end
if any(high[:] .!= 0) 
  snh = sn[high];
  cnh = cn[high];
  dnh = dn[high];
  sn[high] = -1 ./(sqrt(m)*snh);
  cn[high] = im*dnh./(sqrt(m)*snh);
  dn[high] = im*cnh./snh;
end
return sn,cn,dn
end

a,b,c=ellipjc([0.3],0.3,0);
println(a," ",b," ",c)

