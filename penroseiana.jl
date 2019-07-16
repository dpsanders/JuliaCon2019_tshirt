# # Designing the JuliaCon 2019 T-shirt

# ### Cormullion and David P. Sanders

# We offered to design some graphics for the JuliaCon 2019 T-shirt. Several ideas were proposed, and were received variously with “meh” or “cool”, and there was much enjoyable bike-shedding. Eventually, a design derived from the [Penrose tiling](https://en.wikipedia.org/wiki/Penrose_tiling) was developed, using Julia code.

# In this notebook we explain some of the steps in arriving at the final design.

# The basic drawing was made using [Luxor.jl](https://github.com/JuliaGraphics/Luxor.jl), and the algorithms for resolving the colouring used [LightGraphs]( https://github.com/scheinerman/SimpleGraphs.jl).

# ## Penrose tilings

# A Penrose tiling is a tiling (covering) of the plane using two particular tile types.

# ![penrose tiles](penrose-tiles.png)

# The surprising thing is that the tiling turns out to be *quasiperiodic*: it never repeats itself, but does “almost repeat itself”.

# A common method to produce Penrose tilings is the method of **inflation&ndash;deflation**: we start with a given set of tiles and recursively replace each tile with a certain combination of smaller tiles.

# ![recursive division](recursion.svg)

# In this implementation, we build everything out of triangles, represented by the following `struct`. Here, `A`, `B` and `C` are the 3 points in the triangle, `Point` is a type from `Luxor` representing the coordinates of a point in 2D, and `flag` is a Boolean variable representing the “type” of triangle.

using Luxor
using LightGraphs
using Random

struct PenroseTriangle
    flag::Bool
    A::Point
    B::Point
    C::Point
end

# We start by generating $n$ triangles in a circle with center $c$ and radius $r$. The `polar` function generates a vector with given radius and angle

function initialize(center::Point, r=100, n=10)
    penrosetriangles = PenroseTriangle[]
    A = center
    for i in 1:n
        ϕ = (i - 1) * (2π / n)
        C = A + polar(r, ϕ)

        ϕ = (i) * (2π / n)
        B = A + polar(r, ϕ)

        if i % 2 == 1
            triangle = PenroseTriangle(true, B, A, C)
        else
            triangle = PenroseTriangle(true, C, A, B)
        end
        push!(penrosetriangles, triangle)
    end
    return penrosetriangles
end

# Let's draw these triangles. For now we can just apply a polygon-drawing function `poly()` to the triangles. (And we're currently working in global scope, which you'll want to be aware of!)

@svg begin
    for triangle in initialize(O, 150)
        randomhue()
        poly([triangle.A, triangle.B, triangle.C], close=true, :stroke)
    end
end 500 300 "figure1"

# Now we subdivide recursively using the following function:

function subdivide(penrosetriangles::Vector{PenroseTriangle})
    result = PenroseTriangle[]
    for triangle in penrosetriangles
        A, B, C = triangle.A, triangle.B, triangle.C

        if triangle.flag
            Q = A + (B - A) / MathConstants.golden
            R = B + (C - B) / MathConstants.golden
            push!(result, PenroseTriangle(false, R, Q, B))
            push!(result, PenroseTriangle(true,  Q, A, R))
            push!(result, PenroseTriangle(true,  C, A, R))
        else
            P = C + (A - C) / MathConstants.golden
            push!(result, PenroseTriangle(false, B, P, A))
            push!(result, PenroseTriangle(true,  P, C, B))
        end
    end
    return result
end

# Let's see the result of a single subdivision. We can use `foldl()` to run `subdivide()` a number of times.

@svg begin
    penrosetriangles = foldl((x, y) -> subdivide(x), 1:1, init=initialize(O, 150))
    setlinecap("butt")
    setlinejoin("bevel")
    setline(2)
    for triangle in penrosetriangles
        p = [triangle.A, triangle.B, triangle.C]
        ispolyclockwise(p) && reverse!(p)
        poly(p, close=true, :stroke)
    end

    # draw 2 more levels
    penrosetriangles = foldl((x, y) -> subdivide(x), 1:2, init=initialize(O, 150))
    setopacity(0.5)

    for triangle in penrosetriangles
        randomhue()
        p = [triangle.A, triangle.B, triangle.C]
        ispolyclockwise(p) && reverse!(p)
        poly(p, close=true, :fill)
    end
end 500 300 "figure2"

# The more levels you recurse to, the richer and more complex the pattern:

@svg begin
    setlinecap("butt")
    setlinejoin("bevel")
    setline(1)

    # draw first level
    penrosetriangles = foldl((x, y) -> subdivide(x), 1:1, init=initialize(O, 150))

    for triangle in penrosetriangles
        poly([triangle.A, triangle.B, triangle.C], close=true, :stroke)
    end

    # draw 4 more levels
    setopacity(0.5)
    penrosetriangles = foldl((x, y) -> subdivide(x), 1:4, init=initialize(O, 150))

    for triangle in penrosetriangles
        randomhue()
        poly([triangle.A, triangle.B, triangle.C], close=true, :fill)
    end
end 500 300 "figure3"

# This provides the basis for our Penrose-inspired design.

# # Colouring the triangles

# A Penrose tiling can be thought of as a planar map, and the [Four-colour Theorem](https://en.wikipedia.org/wiki/Four_color_theorem) says that any map of the plane can be coloured using at most 4 colours. Here, a map colouring is one in which neighbouring "countries" are coloured by different colours if they are joined by an edge (not just a point).

# Since there are 4 colours in the Julia logo, we decided to find the perfect colouring of our Penrose tiling using the Julia colours!

# ## Mapping to a graph

# To do so, we think of the collection of triangles as a **graph** or **network**, i.e. a collection of vertices joined by edges. We need to find the **adjacency matrix** of the graph, i.e. a matrix $A_{ij}$ which contains $1$ if nodes $i$ and $j$ are connected, and $0$ if they are not connected.

# The problem is that the inflation&ndash;deflation method that we are using to construct the triangles gives us no information about which triangles are neighbours (at least, we are not aware of how to directly extract this information from the construction).

# Hence, the first thing to do is to extract the adjacency matrix, for which we just use brute force: given two triangles, we calculate the distances between each of their vertices. If exactly two of these distances are (approximately) zero, then those two vertices coincide, and hence the two triangles share an edge:

function incidence_matrix(ptriangles::Vector{PenroseTriangle})
    n = length(ptriangles)
    g = SimpleGraph(n)

    for i in 1:n
        for j in i+1:n
            ti = ptriangles[i]
            tj = ptriangles[j]

            points_i = [ti.A, ti.B, ti.C]
            points_j = [tj.A, tj.B, tj.C]

            coincidental = 0
            tol = 1e-2

            for k in 1:3, l in 1:3
                if distance(points_i[k], points_j[l]) < tol
                    coincidental += 1
                end
            end

            if coincidental == 2  # share edge
                add_edge!(g, i, j)
            end
        end
    end
    return g
end

# Use this to create an incidence matrix for our triangles, and we'll draw the graph over the triangles, to show what's going on.

function drawgraph()
    penrosetriangles = foldl((x, y) -> subdivide(x), 1:3, init=initialize(O, 150))

    @svg begin
        setlinecap("butt")
        setlinejoin("bevel")
        setline(1)
        @layer begin
            setopacity(0.9)
            # Julia colors are built-in to Luxor
            juliacolors = [Luxor.julia_green, Luxor.julia_red,Luxor.julia_purple, Luxor.julia_blue]

            # remove the middle ones:
            outertriangles = PenroseTriangle[]
            for t in penrosetriangles
                t1 = [t.A, t.B, t.C]
                pc1 = polycentroid(t1)
                if distance(pc1, O) > 40
                    push!(outertriangles, t)
                end
            end

            @show length(outertriangles)
            # draw them with random colors
            for (k, t) in enumerate(outertriangles)
                sethue(juliacolors[rand(1:end)]...)
                t1 = [t.A, t.B, t.C]
                if !ispolyclockwise(t1)
                    reverse!(t1)
                end
                poly(offsetpoly(t1, -1), :fill)
            end
        end

        # create graph

        g = incidence_matrix(outertriangles)

        # draw graph over triangles

        for edge in edges(g)
            e, f = src(edge), dst(edge)

            st = outertriangles[e]
            ft = outertriangles[f]

            stflag = st.flag

            t1 = [st.A, st.B, st.C]
            t2 = [ft.A, ft.B, ft.C]

            pc1 = polycentroid([st.A, st.B, st.C])
            pc2 = polycentroid([ft.A, ft.B, ft.C])

            if !isapprox(pc1, pc2)
                sethue("white")
                circle(pc1,  1.5, :fill)
                circle(pc2,  1.75, :fill)
                sethue("black")
                arrow(between(pc1, pc2, 0.075), between(pc1, pc2, 0.925), arrowheadlength=5,
                arrowheadangle=π/6, linewidth=1)
            end
        end
    end 420 500 "fig4-network"
end

drawgraph()

# ![network](fig4-network.svg)

# ## Colouring the graph

# Finally, we want to find a **graph colouring** with the property that adjacent nodes do not share the same colour. This is accomplished by the following lines:

# ```
# colorsequence = Random.shuffle!(collect(1:length(outertriangles)))
# colouring = LightGraphs.perm_greedy_color(g, colorsequence)
# ```

# The `perm_greedy_color()` function in LightGraphs uses a “greedy heuristic” to color the graph using the given order.

# Here’s a function to draw the triangles with the “correct” map coloring:

function drawgreedy()
    PenroseTriangles = foldl((x, y) -> subdivide(x), 1:3, init=initialize(O, 150))

    outertriangles = PenroseTriangle[]
    for t in PenroseTriangles
        t1 = [t.A, t.B, t.C]
        pc1 = polycentroid(t1)
        if distance(pc1, O) > 40
            push!(outertriangles, t)
        end
    end

    g = incidence_matrix(outertriangles)

    colorsequence = Random.shuffle!(collect(1:length(outertriangles)))
    colouring = LightGraphs.perm_greedy_color(g, colorsequence)

    @svg begin
        juliacolors = [Luxor.julia_green, Luxor.julia_red,Luxor.julia_purple, Luxor.julia_blue]
        for (k, t) in enumerate(outertriangles)
            sethue(juliacolors[colouring.colors[k]]...)
            t1 = [t.A, t.B, t.C]
            if !ispolyclockwise(t1)
                reverse!(t1)
            end
            poly(offsetpoly(t1, -1), :fill)
        end
    end 500 300 "fig5-greedy"
end

drawgreedy()

# It turns out that for the Penrose tiling, there are rather few nodes that use the fourth colour.

# To make the result visually more interesting, we randomly choose approximately a quarter of the nodes to have the 4th colour, provided that none of its neighbours have that colour:

function tweak()
    PenroseTriangles = foldl((x, y) -> subdivide(x), 1:3, init=initialize(O, 150))
    outertriangles = PenroseTriangle[]
    for t in PenroseTriangles
        t1 = [t.A, t.B, t.C]
        pc1 = polycentroid(t1)
        if distance(pc1, O) > 40
            push!(outertriangles, t)
        end
    end
    g = incidence_matrix(outertriangles)

    colorsequence = Random.shuffle!(collect(1:length(outertriangles)))
    colouring = LightGraphs.perm_greedy_color(g, colorsequence)
    n = length(colouring.colors)

    for i in 1:n ÷ 4
        which = rand(1:n)
        if 4 ∉ [colouring.colors[x] for x in neighbors(g, which)]
            colouring.colors[which] = 4
        end
    end

    # map(x -> sum(colouring.colors .== x ), 1:4)

    @svg begin
        juliacolors = [Luxor.julia_green, Luxor.julia_red,Luxor.julia_purple, Luxor.julia_blue]
        for (k, t) in enumerate(outertriangles)
            sethue(juliacolors[colouring.colors[k]]...)
            t1 = [t.A, t.B, t.C]
            if !ispolyclockwise(t1)
                reverse!(t1)
            end
            poly(offsetpoly(t1, -1), :fill)
        end
    end 500 300 "fig6-tweak"
end

tweak()

using Literate #src
Literate.notebook("penroseiana.jl", ".") #src
Literate.markdown("penroseiana.jl", ".", name="penrosiana", documenter=false) #src
