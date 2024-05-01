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


L=0.5;
K,Kp=ellipkkp.([0.5,0.75])
println(K," ",Kp)
