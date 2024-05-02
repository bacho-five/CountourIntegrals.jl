using Distributed
using SlurmClusterManager
addprocs(SlurmManager())
println(workers())
using Plots
using MAT
A = matread("bcsstk01.mat")["Problem"]["A"]
A = Matrix(A)

@everywhere begin
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
function utile(N, Kp, K, L, m, M, Id, A, b, f)
    Y = SharedArray{ComplexF64}(40)
    @everywhere begin
        N1 = $N
        t = 0.5im .* Kp .- K .+ (0.5:N1) .* 2 .* K ./ N1
        u, cn, dn = ellipjc(t, L, 0)
        w = (m * M)^(1 / 4) * ((1 / k .+ u) ./ (1 / k .- u))
        dzdt = cn .* dn ./ (1 / k .- u) .^ 2
    end
    @sync @distributed for j = 1:N1
        po = (b' * A * (((w[j]^2) * Id - A) \ b))
        Y[j] += (f(w[j]^2) / w[j]) * po * dzdt[j]
    end
    S = sum(Y)
    S = -8 * K * (m * M)^(1 / 4) * imag(S) / (k * pi * N)
    error = norm(X .- S) / norm(X)
    println(N1, " ", error)
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
        # X = (b' * X * b) / (b' * b)
        #display(X)
        # X = 0
        k = ((M / m)^(1 / 4) - 1) / ((M / m)^(1 / 4) + 1)
        L = -log(k) / pi
        K, Kp = ellipkkp(L)
    end
    #Y = SharedArray{ComplexF64}(20);
    for N = 5:5:15
        o = Int64(N / 5)
        Z[o] = @elapsed utile(N, Kp, K, L, m, M, Id, A, b, f)
    end
    return Z
end
@everywhere begin
    a = [1 2 3 4 5 6 8]
    A = Diagonal(vec(a))
    n = 7
    b = zeros(n, 1)
    b[1] = 1
end
Z = zeros(4)
progettoparallelo(A, vec(b))
x = range(5, 15, length=3);
plot(x, Z)
png("solutions")
