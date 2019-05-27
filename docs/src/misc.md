Miscellaneous
============
 
# Logo
-------
Code used to generate the `msJ`logo. 

```julia
using Plots
gr()

@. model(x, p) = p[4] + p[3] * exp(- ( (x-p[2])/p[1] )^2)

x1 = range(-7.5, stop=-2.5, length=100)
p1 = [1, -5, 1, 0]
y1 = model(x1,p1);
plot(x1,y1, color = RGBA(0.884, 0.2, 0.2), fillrange = 0, fillalpha = 0.9, label = "", 
	thickness_scaling = 2.0, 
    xaxis=nothing,
    yaxis=nothing,
    background_color=:transparent, foreground_color=:black,)

x2 = range(-2.5, stop=2.5, length=100)
p2 = [1, 0, 2, 0]
y2 = model(x2, p2)
plot!(x2,y2, color= RGBA(0.22, 0.596, 0.149), fillrange = 0, fillalpha = 0.9, label = "")

x3 = range(2.5, stop=7.5, length=100)
p3 = [1, 5, 1, 0]
y3 = model(x3,p3)
plot!(x3,y3,color = RGBA(0.584, 0.345, 0.608), fillrange = 0, fillalpha = 0.9, label = "")

x4 = range(-17.5, stop=-12.5, length=100)
p4 = [1, -15, 1, 0]
y4 = model(x4,p4);
plot!(x4,y4, color = RGBA(0.255, 0.412, 0.882), fillrange = 0, fillalpha = 0.9, label = "")

savefig("docs/src/assets/logo.png")
```

![logo](assets/logo.png)
