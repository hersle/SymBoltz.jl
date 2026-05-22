nodes_trapezoidal(N) = collect(range(-1, 1, length = N))
weights_trapezoidal(N) = [i == 1 || i == N ? 1/(N-1) : 2/(N-1) for i in 1:N]
transform_quadrature_domain(x, a, b) = (b+a)/2 .+ (b-a)/2 .* x