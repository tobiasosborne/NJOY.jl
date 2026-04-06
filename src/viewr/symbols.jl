#=
    Symbol and circle drawing routines for the NJOY graph engine.
    Faithful port of graph.f90: dsym (25 symbol types), circle.
=#

using Printf

"""
    circle!(gs, io, x, y, r, w)

Draw a circle via PostScript arc. Coordinates in page units.
Matches graph.f90 `circle` (lines 2763-2793).
"""
function circle!(gs::GraphState, io::IO, x::Float64, y::Float64,
                 r::Float64, w::Float64)
    if gs.land == 1
        u1 = gs.uwidth - 72 * y + gs.ushift
        v1 = 72 * x + gs.vshift
    else
        u1 = 72 * x + gs.ushift
        v1 = 72 * y + gs.vshift
    end
    rr = 72 * r

    println(io, "stroke")
    println(io, "newpath")
    u1 = clamp(u1, -1000.0, 2000.0)
    v1 = clamp(v1, -1000.0, 2000.0)
    @printf(io, " %9.2f%9.2f%9.2f 0 360 arc\n", u1, v1, rr)
end

"""
    dsym!(gs, io, x, y, isym, size)

Draw one of 25 plotting symbols at position (x, y) in page coordinates.
Matches graph.f90 `dsym` (lines 2376-2761).
"""
function dsym!(gs::GraphState, io::IO, x::Float64, y::Float64,
               isym::Int, size::Float64)
    w = 0.01      # stroke width
    rr = 1.128    # circle radius multiplier
    c8 = 0.414    # octagon chamfer
    cd = 1.4      # diamond scale
    cx = 0.82     # X-stroke on circle
    t1 = 1.75     # triangle proportions
    t2 = 1.52
    t3 = 0.88
    sz = 6 * size / 10
    foregr = 0.0
    backgr = 1.0

    if isym == 0  # square
        xx = [x-sz, x+sz, x+sz, x-sz, x-sz]
        yy = [y+sz, y+sz, y-sz, y-sz, y+sz]
        poly2!(gs, io, xx, yy, 5, w, backgr)
    elseif isym == 1  # octagon
        xx = [x-c8*sz, x+c8*sz, x+sz, x+sz, x+c8*sz, x-c8*sz, x-sz, x-sz, x-c8*sz]
        yy = [y+sz, y+sz, y+c8*sz, y-c8*sz, y-sz, y-sz, y-c8*sz, y+c8*sz, y+sz]
        poly2!(gs, io, xx, yy, 9, w, backgr)
    elseif isym == 2  # triangle
        xx = [x-t2*sz, x, x+t2*sz, x-t2*sz]
        yy = [y-t3*sz, y+t1*sz, y-t3*sz, y-t3*sz]
        poly2!(gs, io, xx, yy, 4, w, backgr)
    elseif isym == 3  # cross (+)
        moveh!(gs, io, x, y+sz)
        drawh!(gs, io, x, y-sz, w, 0)
        moveh!(gs, io, x-sz, y)
        drawh!(gs, io, x+sz, y, w, 0)
    elseif isym == 4  # ex (X)
        moveh!(gs, io, x-sz, y+sz)
        drawh!(gs, io, x+sz, y-sz, w, 0)
        moveh!(gs, io, x-sz, y-sz)
        drawh!(gs, io, x+sz, y+sz, w, 0)
    elseif isym == 5  # diamond
        xx = [x-cd*sz, x, x+cd*sz, x, x-cd*sz]
        yy = [y, y+cd*sz, y, y-cd*sz, y]
        poly2!(gs, io, xx, yy, 5, w, backgr)
    elseif isym == 6  # inverted triangle
        xx = [x-t2*sz, x, x+t2*sz, x-t2*sz]
        yy = [y+t3*sz, y-t1*sz, y+t3*sz, y+t3*sz]
        poly2!(gs, io, xx, yy, 4, w, backgr)
    elseif isym == 7  # exed square
        xx = [x-sz, x+sz, x+sz, x-sz, x-sz]
        yy = [y+sz, y+sz, y-sz, y-sz, y+sz]
        poly2!(gs, io, xx, yy, 5, w, backgr)
        moveh!(gs, io, x-sz, y+sz)
        drawh!(gs, io, x+sz, y-sz, w, 0)
        moveh!(gs, io, x-sz, y-sz)
        drawh!(gs, io, x+sz, y+sz, w, 0)
    elseif isym == 8  # crossed ex
        moveh!(gs, io, x-sz, y+sz)
        drawh!(gs, io, x+sz, y-sz, w, 0)
        moveh!(gs, io, x-sz, y-sz)
        drawh!(gs, io, x+sz, y+sz, w, 0)
        moveh!(gs, io, x-sz, y)
        drawh!(gs, io, x+sz, y, w, 0)
        moveh!(gs, io, x, y+sz)
        drawh!(gs, io, x, y-sz, w, 0)
    elseif isym == 9  # crossed diamond
        xx = [x-cd*sz, x, x+cd*sz, x, x-cd*sz]
        yy = [y, y+cd*sz, y, y-cd*sz, y]
        poly2!(gs, io, xx, yy, 5, w, backgr)
        moveh!(gs, io, x, y+cd*sz)
        drawh!(gs, io, x, y-cd*sz, w, 0)
        moveh!(gs, io, x-cd*sz, y)
        drawh!(gs, io, x+cd*sz, y, w, 0)
    elseif isym == 10  # crossed octagon
        xx = [x-c8*sz, x+c8*sz, x+sz, x+sz, x+c8*sz, x-c8*sz, x-sz, x-sz, x-c8*sz]
        yy = [y+sz, y+sz, y+c8*sz, y-c8*sz, y-sz, y-sz, y-c8*sz, y+c8*sz, y+sz]
        poly2!(gs, io, xx, yy, 9, w, backgr)
        moveh!(gs, io, x, y+sz)
        drawh!(gs, io, x, y-sz, w, 0)
        moveh!(gs, io, x-sz, y)
        drawh!(gs, io, x+sz, y, w, 0)
    elseif isym == 11  # double triangle
        xx = [x-t2*sz, x, x+t2*sz, x-t2*sz]
        yy = [y-t3*sz, y+t1*sz, y-t3*sz, y-t3*sz]
        poly2!(gs, io, xx, yy, 4, w, backgr)
        xx2 = [x-t2*sz, x, x+t2*sz, x-t2*sz]
        yy2 = [y+t3*sz, y-t1*sz, y+t3*sz, y+t3*sz]
        poly2!(gs, io, xx2, yy2, 4, w, backgr)
    elseif isym == 12  # crossed square
        xx = [x-sz, x+sz, x+sz, x-sz, x-sz]
        yy = [y+sz, y+sz, y-sz, y-sz, y+sz]
        poly2!(gs, io, xx, yy, 5, w, backgr)
        moveh!(gs, io, x, y+sz)
        drawh!(gs, io, x, y-sz, w, 0)
        moveh!(gs, io, x-sz, y)
        drawh!(gs, io, x+sz, y, w, 0)
    elseif isym == 13  # exed octagon
        xx = [x-c8*sz, x+c8*sz, x+sz, x+sz, x+c8*sz, x-c8*sz, x-sz, x-sz, x-c8*sz]
        yy = [y+sz, y+sz, y+c8*sz, y-c8*sz, y-sz, y-sz, y-c8*sz, y+c8*sz, y+sz]
        poly2!(gs, io, xx, yy, 9, w, backgr)
        moveh!(gs, io, x-sz, y+sz)
        drawh!(gs, io, x+sz, y-sz, w, 0)
        moveh!(gs, io, x-sz, y-sz)
        drawh!(gs, io, x+sz, y+sz, w, 0)
    elseif isym == 14  # triangle + square
        xx = [x-sz, x+sz, x+sz, x-sz, x-sz]
        yy = [y+sz, y+sz, y-sz, y-sz, y+sz]
        poly2!(gs, io, xx, yy, 5, w, backgr)
        moveh!(gs, io, x-sz, y-sz)
        drawh!(gs, io, x, y+sz, w, 0)
        drawh!(gs, io, x+sz, y-sz, w, 0)
    elseif isym == 15  # filled circle
        circle!(gs, io, x, y, rr*sz, w)
        fillh!(gs, io, foregr)
    elseif isym == 16  # open circle
        circle!(gs, io, x, y, rr*sz, w)
        fillh!(gs, io, backgr)
    elseif isym == 17  # open square
        xx = [x-sz, x+sz, x+sz, x-sz, x-sz]
        yy = [y+sz, y+sz, y-sz, y-sz, y+sz]
        poly2!(gs, io, xx, yy, 5, w, backgr)
    elseif isym == 18  # filled square
        xx = [x-sz, x+sz, x+sz, x-sz, x-sz]
        yy = [y+sz, y+sz, y-sz, y-sz, y+sz]
        poly2!(gs, io, xx, yy, 5, w, foregr)
    elseif isym == 19  # filled diamond
        xx = [x-cd*sz, x, x+cd*sz, x, x-cd*sz]
        yy = [y, y+cd*sz, y, y-cd*sz, y]
        poly2!(gs, io, xx, yy, 5, w, foregr)
    elseif isym == 20  # filled triangle
        xx = [x-t2*sz, x, x+t2*sz, x-t2*sz]
        yy = [y-t3*sz, y+t1*sz, y-t3*sz, y-t3*sz]
        poly2!(gs, io, xx, yy, 4, w, foregr)
    elseif isym == 21  # filled inverted triangle
        xx = [x-t2*sz, x, x+t2*sz, x-t2*sz]
        yy = [y+t3*sz, y-t1*sz, y+t3*sz, y+t3*sz]
        poly2!(gs, io, xx, yy, 4, w, foregr)
    elseif isym == 22  # crossed circle
        circle!(gs, io, x, y, rr*sz, w)
        fillh!(gs, io, backgr)
        moveh!(gs, io, x, y+rr*sz)
        drawh!(gs, io, x, y-rr*sz, w, 0)
        moveh!(gs, io, x-rr*sz, y)
        drawh!(gs, io, x+rr*sz, y, w, 0)
    elseif isym == 23  # exed circle
        circle!(gs, io, x, y, rr*sz, w)
        fillh!(gs, io, backgr)
        moveh!(gs, io, x-cx*sz, y+cx*sz)
        drawh!(gs, io, x+cx*sz, y-cx*sz, w, 0)
        moveh!(gs, io, x-cx*sz, y-cx*sz)
        drawh!(gs, io, x+cx*sz, y+cx*sz, w, 0)
    elseif isym == 24  # exed diamond
        xx = [x-cd*sz, x, x+cd*sz, x, x-cd*sz]
        yy = [y, y+cd*sz, y, y-cd*sz, y]
        poly2!(gs, io, xx, yy, 5, w, backgr)
        moveh!(gs, io, x-cd*sz/2, y+cd*sz/2)
        drawh!(gs, io, x+cd*sz/2, y-cd*sz/2, w, 0)
        moveh!(gs, io, x-cd*sz/2, y-cd*sz/2)
        drawh!(gs, io, x+cd*sz/2, y+cd*sz/2, w, 0)
    else  # default: circle
        circle!(gs, io, x, y, rr*sz, w)
        fillh!(gs, io, backgr)
    end
end
