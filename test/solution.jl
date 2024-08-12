using Test

idxs = [M.g.a, M.g.a, M.g.a, M.g.a]
ts = [1.0, 2.0, 3.0]
ks = [1.0, 10.0] ./ SymBoltz.k0

has_size(arr::Array, siz) = siz == size(arr)
has_size(arr::Number, siz) = siz == ()

# size of solution output should match input arguments in order
#@test has_size(sol(ts[1], idxs[1]), ())
@test has_size(sol(ts[1], idxs[1]), ())
@test has_size(sol(ts[1], idxs), (length(idxs),))
@test has_size(sol(ts, idxs[1]), (length(ts),))
@test has_size(sol(ts, idxs), (length(ts), length(idxs)))

@test has_size(sol(ks[1], ts[1], idxs[1]), ())
@test has_size(sol(ks[1], ts[1], idxs), (length(idxs),))
@test has_size(sol(ks[1], ts, idxs[1]), (length(ts),))
@test has_size(sol(ks[1], ts, idxs), (length(ts), length(idxs)))
@test has_size(sol(ks, ts[1], idxs[1]), (length(ks),))
@test has_size(sol(ks, ts[1], idxs), (length(ks), length(idxs)))
@test has_size(sol(ks, ts, idxs[1]), (length(ks), length(ts)))
@test has_size(sol(ks, ts, idxs), (length(ks), length(ts), length(idxs)))