# julia 1.2
using Luxor, ColorSchemes, Colors

using LightGraphs

struct PenroseTriangle
    red::Bool # not red, more of a flag
    pointA::Point
    pointB::Point
    pointC::Point
end

function PenroseTiling(centerpos::Point, radius, depth=4;
        type=:P3)
    triangles = PenroseTriangle[]
    A = centerpos
    n = 10 # a circle of triangles 360/10 -> 36°
    for i in 1:n
        phi = (i - 1) * (2π / n)
        C = A + polar(radius, phi)
        phi = (i) * (2π / n)
        B = A + polar(radius, phi)
        if type == :P3
            if i % 2 == 1
                triangle = PenroseTriangle(true, A, C, B)
            else
                triangle = PenroseTriangle(true, A, B, C)
            end
        else # P2
            if i % 2 == 1
                triangle = PenroseTriangle(true, B, A, C)
            else
                triangle = PenroseTriangle(true, C, A, B)
            end
        end
        push!(triangles, triangle)
    end
    for i in 1:depth
        triangles = subdivide(triangles, type=type)
    end
    return triangles
end

function subdivide(triangles;
        type=:P3)
    result = []
    for triangle in triangles
        A, B, C = triangle.pointA, triangle.pointB, triangle.pointC
        if triangle.red == true
            if type == :P3
                #  P3 rhombus
                P = A + (B - A) / MathConstants.golden
                push!(result, PenroseTriangle(true, C, P, B))
                push!(result, PenroseTriangle(false, P, C, A))
            else # P2 half kite
                Q = A + (B - A) / MathConstants.golden
                R = B + (C - B) / MathConstants.golden
                push!(result, PenroseTriangle(false, R, Q, B))
                push!(result, PenroseTriangle(true,  Q, A, R))
                push!(result, PenroseTriangle(true,  C, A, R))
            end
        else
            if type == :P3
                # P3 rhombus
                Q = B + (A - B) / MathConstants.golden
                R = B + (C - B) / MathConstants.golden
                push!(result, PenroseTriangle(true, R, Q, A))
                push!(result, PenroseTriangle(false, R, A, C))
                push!(result, PenroseTriangle(false, Q, R, B))
            else # P2 kite/dart
                P = C + (A - C) / MathConstants.golden
                push!(result, PenroseTriangle(false, B, P, A))
                push!(result, PenroseTriangle(true,  P, C, B))
            end
        end
    end
    return result
end

function Base.convert(::Type{Vector{Point}},
    pt::PenroseTriangle)
    return [pt.pointA,pt.pointB,pt.pointC]
end

function drawtiles(radius, depth, type, foregroundcolor;
        coloringstyle=1)

    @layer begin
        rotate(π/2)
        juliacolors = [Luxor.julia_blue, Luxor.julia_red,Luxor.julia_purple,Luxor.julia_green]
        triangles = PenroseTiling(O, radius, depth, type=type)

        ### add graph-based coloring style
        if coloringstyle==5
            g = incidence_graph(triangles)
            coloring = LightGraphs.perm_greedy_color(g, [1:length(triangles);])

            # make coloring more interesting by adding random greens:
            n = length(triangles)
            for i in 1:n ÷ 4
                which = rand(1:n)
                if 4 ∉ [coloring.colors[x] for x in neighbors(g, which)]
                    coloring.colors[which] = 4
                end
            end
        end

        for (n, triangle) in enumerate(triangles)
            pgon = [triangle.pointA, triangle.pointB, triangle.pointC]
            # we need clockwise polys
            if !ispolyclockwise(pgon)
                pgon = reverse(pgon)
            end
            # find distance from center to the innermost point of this triangle
            d = minimum([distance(O, pgon[1]), distance(O, pgon[2]), distance(O, pgon[3])])
            if d > 148 # <- !
                if coloringstyle == 1
                    sethue(LCHab(80, 80, rand(0:359)))
                elseif coloringstyle == 2
                    sl = slope(O, polycentroid(pgon))
                    sethue(LCHab(60, 100, rescale(sl, 0, π, 0, 360)))
                elseif coloringstyle == 3
                    sethue(juliacolors[rand(1:end)])
                elseif coloringstyle == 4
                    sethue(juliacolors[mod1(n, end)])
                elseif coloringstyle == 5
                    sethue(juliacolors[coloring.colors[n]])
                else
                    sethue(foregroundcolor)
                end
                polysmooth(offsetpoly(pgon, rescale(depth, 3, 6, -2, -1)), 0.5,  :fillstroke)
            end
        end
    end
end

function incidence_graph(triangles)
    n = length(triangles)
    g = SimpleGraph(n)
    for i in 1:n
        for j in i+1:n
            ti = triangles[i]
            tj = triangles[j]

            points_i = [ti.pointA, ti.pointB, ti.pointC]
            points_j = [tj.pointA, tj.pointB, tj.pointC]

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

function makeone(w, h, fname;
        type=:P3,
        depth=3,
        backgroundcolor = "",
        foregroundcolor = "red",
        style=1)
        Juno.clearconsole()

    Drawing(w, h, fname)
    origin()
    if backgroundcolor != ""
        background(backgroundcolor)
    end
    thickness = 1.5
    radius=w/2
    setline(thickness)

    # create clipping mask
    box(O, w, h, :path)
    newsubpath()
    circlepath(O, radius - 115, :path)
    clip()
    drawtiles(radius, depth, type, foregroundcolor, coloringstyle=style)
    clipreset()

    @layer begin
        sethue(foregroundcolor)
        toptext, bottomtext =  "JULIACON 2019", "BALTIMORE • USA"
        fontsize(26)
        fontface("ChunkFive")
        # these two radius values usually need adjusting by eye because fonts ...
        textcurvecentered(toptext,   -π/2, 119, O, letter_spacing=6)
        textcurvecentered(bottomtext, π/2, 136, O, clockwise=false, letter_spacing=5)
    end

    circlepath(O, 105, :path)
    clip()
    @layer begin
        scale(0.55)
        translate(-178, -125)
        julialogo(color=true, bodycolor=foregroundcolor)
    end
    clipreset()
    finish()
    preview()
end

function makeall()
    for fc_bc_pair in (("red", "white"), ("black", "white"), ("orange", "black"), ("white", ""))
        for i in 2:6
            for style in 1:5
                makeone(512, 512, "/tmp/design-typeP2-depth$(i)-$(join(fc_bc_pair))-$(style).svg", type=:P2, depth=i, foregroundcolor=first(fc_bc_pair), backgroundcolor=last(fc_bc_pair), style=style)
            end
        end
    end
end

# makeall()

makeone(512, 512, "/tmp/design.svg", type=:P2, depth=6, foregroundcolor="white", backgroundcolor="black", style=5)

makeone(512, 512, "/tmp/design.pdf", type=:P2, depth=6, foregroundcolor="white", backgroundcolor="black", style=5)

makeone(512, 512, "/tmp/penrose_level_depth4.pdf", type=:P2, depth=4, foregroundcolor="white", backgroundcolor="", style=5)
