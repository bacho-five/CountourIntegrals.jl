using Polynomials
using LinearAlgebra
using Elliptic
function ellipkkp(L)
    if L > 10
        K = pi / 2
        Kp = pi * L + log(4)
        return K, Kp
    end
    m = exp(-2 * pi * L)
    a0 = 1
    b0 = sqrt(1 - m)
    s0 = m
    i1 = 0
    mm = 1
    a1 = (a0 + b0) / 2
    while mm > eps() # eps funziona?
        a1 = (a0 + b0) / 2
        b1 = sqrt(a0 .* b0)
        c1 = (a0 - b0) / 2
        i1 = i1 + 1
        w1 = 2^i1 * c1 .^ 2
        mm = max(max(w1))
        s0 = s0 + w1
        a0 = a1
        b0 = b1
    end
    K = pi ./ (2 * a1)
    # im = find(m==1);
    # if isempty(im)
    #   K[im] = K[im]*inf;
    # end
    a0 = 1
    b0 = sqrt(m)
    s0 = 1 - m
    i1 = 0
    mm = 1
    while mm > eps() # eps() funziona?
        a1 = (a0 + b0) / 2
        b1 = sqrt(a0 .* b0)
        c1 = (a0 - b0) / 2
        i1 = i1 + 1
        w1 = 2^i1 * c1 .^ 2
        mm = max(max(w1))
        s0 = s0 + w1
        a0 = a1
        b0 = b1
    end
    Kp = pi ./ (2 * a1)
    # im = find(m==0);
    # if isempty(im)
    #   Kp[im] = Kp[im]*inf;
    # end
    return K, Kp
end
function utile2(Kp, K, N, L, m, M, k, Id, A, b, f, X)
    t = 0.5im .* Kp .- K .+ (0.5:N) .* 2 .* K ./ N
    u, cn, dn = Elliptic.ellipj(t, L)
    w = (m * M)^(1 / 4) * ((1 / k .+ u) ./ (1 / k .- u))
    dzdt = cn .* dn ./ (1 / k .- u) .^ 2
    S = 0
    for j = 1:N
        po = (b' * A * (((w[j]^2) * Id - A) \ b))
        S = S + (f(w[j]^2) / w[j]) * po * dzdt[j]
        S = -8 * K * (m * M)^(1 / 4) * imag(S) / (k * pi * N)
    end
    error = norm(S - X) / norm(X)
    #X = X / N
    #println(S ./ X)
    #println(w)
    println(N, " ", error)
end

function progetto(b)
    f = sqrt
    n = size(b, 1)
    a = ones(n - 1, 1)
    a = -a
    a = vec(a)
    c = ones(n, 1)
    c = vec(c)
    Id = Diagonal(c)
    c = 2 * c
    b = vec(b)
    A = Tridiagonal(a, c, a)
    e = eigvals(A)
    m = minimum(e)
    M = maximum(e)
    # M = 4
    # m = 1 / n^2
    X = sqrt(A)
    X = (b' * X * b) / (b' * b)
    #display(X)
    # X = 0
    k = ((M / m)^(1 / 4) - 1) / ((M / m)^(1 / 4) + 1)
    L = -log(k) / pi
    K, Kp = ellipkkp(L)
    utile2(Kp, K, 5, L, m, M, k, Id, A, b, f, X)
    for N = 5:5:20
        o = Int64(N / 5)
        Z[o] = @elapsed utile2(Kp, K, N, L, m, M, k, Id, A, b, f, X)
    end
end
n = 30
b = zeros(n)
b[1] = 1
Z = zeros(40)
progetto(vec(b))