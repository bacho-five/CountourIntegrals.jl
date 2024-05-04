using Distributed
addprocs(10)
using Plots

@everywhere begin
    using MAT
    A = matread("matrice1.mat")["Problem"]["A"]
    A = Matrix(A)
    using Polynomials
    using LinearAlgebra
    using SharedArrays

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

    function ellipjc(u, L, flag)
        if flag == 1
            high = zeros(size(u))
            m = L
        else
            K, Kp = ellipkkp(L)
            high = imag(u) .> Kp / 2
            u[high] = im * Kp .- u[high]
            m = exp(-2 * pi * L)
        end
        if m < 4 * eps()
            sinu = sin.(u)
            cosu = cos.(u)
            sn = sinu + m / 4 * (sinu .* cosu - u) .* cosu
            cn = cosu + m / 4 * (-sinu .* cosu + u) .* sinu
            dn = 1 .+ m / 4 * (cosu .^ 2 .- sinu .^ 2 .- 1)
        else
            if m > 1e-3
                kappa = (1 - sqrt(1 - m)) / (1 + sqrt(1 - m))
            else
                pol = Polynomial([0, 1, 2, 5, 14, 42, 132])
                kappa = pol(m / 4)
            end
            mu = kappa^2
            v = u / (1 + kappa)
            sn1, cn1, dn1 = ellipjc(v, mu, 1)
            denom = (1 .+ kappa * sn1 .^ 2)
            sn = (1 .+ kappa) * sn1 ./ denom
            cn = cn1 .* dn1 ./ denom
            dn = (1 .- kappa .* sn1 .^ 2) ./ denom
        end
        if any(high[:] .!= 0)
            snh = sn[high]
            cnh = cn[high]
            dnh = dn[high]
            sn[high] = -1 ./ (sqrt(m) * snh)
            cn[high] = im * dnh ./ (sqrt(m) * snh)
            dn[high] = im * cnh ./ snh
        end
        return sn, cn, dn
    end
end
function utile(X, Y, b, A, w, Id, f, dzdt, K, m, M, N, k)
    @sync @async for j = 1:N
        po = (b' * A * (((w[j]^2) * Id - A) \ b))
        Y[j] += (f(w[j]^2) / w[j]) * po * dzdt[j]
    end
    S = sum(Y)
    S = -8 * K * (m * M)^(1 / 4) * imag(S) / (k * pi * N)
    error = norm(X .- S) / norm(X)
    println(N, " ", error)
end
function progettoparallelo(A, b)
    @everywhere begin
        f = sqrt
        n = size(b, 1)
        c = ones(n, 1)
        Id = Diagonal(vec(c))
        e = eigvals(A)
        m = minimum(e)
        M = maximum(e)
        #M = 4
        #m = 1 / n^2
        X = sqrt(A)
        X = (b' * X * b)
        #X = (b' * X * b) / (b' * b)
        #display(X)
        # X = 0
        k = ((M / m)^(1 / 4) - 1) / ((M / m)^(1 / 4) + 1)
        L = -log(k) / pi
        K, Kp = ellipkkp(L)
    end
    for N = 5:5:20
        Y = SharedArray{ComplexF64}(20)
        @everywhere begin
            N1 = $N
            t = 0.5im .* Kp .- K .+ (0.5:N1) .* 2 .* K ./ N1
            u, cn, dn = ellipjc(t, L, 0)
            w = (m * M)^(1 / 4) * ((1 / k .+ u) ./ (1 / k .- u))
            dzdt = cn .* dn ./ (1 / k .- u) .^ 2
        end
       # utile(X, Y, b, A, w, Id, f, dzdt, K, m, M, N, k)
        o = Int64(N / 5)
        Z2[o] = @elapsed utile(X, Y, b, A, w, Id, f, dzdt, K, m, M, N, k)
    end
    return Z
end
function utile(X, S, b, A, w, Id, f, dzdt, K, m, M, N, k)
    for j = 1:N
        po = (b' * A * (((w[j]^2) * Id - A) \ b))
        S = S + (f(w[j]^2) / w[j]) * po * dzdt[j]
    end
    S = -8 * K * (m * M)^(1 / 4) * imag(S) / (k * pi * N)
    error = norm(S - X) / norm(X)
    X = X / N
    # println(S ./ X)
    println(N, " ", error)
end

function progetto(A, b)
    f = sqrt
    n = size(b, 1)
    c = ones(n, 1)
    c = vec(c)
    Id = Diagonal(c)
    e, _ = eigen(A)
    m = minimum(e)
    M = maximum(e)
    #M = 4
    #m = 1 / n^2
    X = sqrt(A)
    X = (b' * X * b)
    #X = (b' * X * b) / (b' * b)
    #display(X)
    # X = 0
    k = ((M / m)^(1 / 4) - 1) / ((M / m)^(1 / 4) + 1)
    L = -log(k) / pi
    K, Kp = ellipkkp(L)
    #utile2(Kp, K, 5, L, m, M, k, Id, A, b, f, X)
    for N = 5:5:20
        t = 0.5im .* Kp .- K .+ (0.5:N) .* 2 .* K ./ N
        u, cn, dn = ellipjc(t, L, 0)
        w = (m * M)^(1 / 4) * ((1 / k .+ u) ./ (1 / k .- u))
        dzdt = cn .* dn ./ (1 / k .- u) .^ 2
        S = 0
        o = Int64(N / 5)
        Z1[o] = @elapsed utile(X, S, b, A, w, Id, f, dzdt, K, m, M, N, k)
    end
end
@everywhere begin
    n = size(A, 1)
    b = zeros(n, 1)
    b[1] = 1
end
Z1 = zeros(4, 1)
Z2 = zeros(4, 1)
progetto(A,b)
progettoparallelo(A, vec(b))
x = 1:4
plot(x, [Z1, Z2], title="Comparazione dei tempi", label=["Tempi con un nodo" "Tempi in parallelo"])
png("s1.6")
