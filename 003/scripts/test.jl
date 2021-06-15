using Plots; gr()
size = 100
x = range(-2*pi, stop=2*pi, length=size)
y = range(-2*pi, stop=2*pi, length=size)

r2 = Array{Float64}(undef,size,size)
for i = 1:size
  for j = 1:size
    r2[i,j] = (x[i]^2 + y[j]^2)
        # z(i,j) = sin(x(i))*cos(y(j))*sin(r2)/log(r2+1)
  end
end


contourf(x, y, r2, framestyle=:box, minorticks=true, tick_dir=:out)

# savefig(pl, "test.tikz")