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
  high = zeros(size(u));
  m = L;
end
if m < 4*eps()
  sinu = sin(u);
  cosu = cos(u);
  sn = sinu + m/4*(sinu.*cosu-u).*cosu;
  cn = cosu + m/4*(-sinu.*cosu+u).*sinu;
  dn = 1 + m/4*(cosu.^2-sinu.^2-1);
else
  if m > 1e-3
    kappa = (1-sqrt(1-m))/(1+sqrt(1-m));
  else
    kappa = polyval([132,42,14,5,2,1,0],m/4);
  end
  mu = kappa^2;
  v = u/(1+kappa);
  [sn1,cn1,dn1] = ellipjc(v,mu,1);
  denom = (1+kappa*JC[1].^2);
  sn = (1+kappa)*JC[1] ./ denom;
  cn = JC[2].*JC[3] ./ denom;
  dn = (1-kappa*JC[1].^2) ./ denom;
end
if any(high((:)) %rivedere
  snh = sn(high);
  cnh = cn(high);
  dnh = dn(high);
  sn(high) = -1./(sqrt(m)*snh);
  cn(high) = i*dnh./(sqrt(m)*snh);
  dn(high) = i*cnh./snh;
end
end

function ellipjc(u,L)
    KK = ellipkkp(L); 
   high = imag(u) > KK[2]/2;
  u[high] = i*KK[2] - u[high];
     m = exp(-2*pi*L); 
