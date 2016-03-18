# function ring()
#     nodes = [PointMass(0.1) for i in 1:10]
#     edges = []
#     for j = 1:10
#         for k = j+1:10
#             push!(edges, [j,k])
#         end
#     end
#     # edges = Any[[i, i+1] for i in 1:9]
#     # push!(edges, [10,1])
#     edges = [DampedSpring(e[1], e[2], 2000.0, 4.0, 0.05) for e in edges]
#     gravity = [0, -9.81]
#     r = SoftRobot(nodes, edges)

#     positions = [[i/40 + 0.75; 0.8 + rand(1)*0.01 + 0.05 * mod(i,2)] for (i,n) in enumerate(r.nodes)]
#     velocities = [zeros(2) for n in r.nodes]
#     state = SoftRobotState(positions, velocities)

#     r, state
# end

#   2 5 8
# 1 3 6 9
#   4 7 10


function tetrahedron()
    # k = 4000.0
    k = 4000.0
    b = 4.0
    l = 0.05
    m = 0.1

    nodes = [PointMass(m) for i in 1:4]
    edges = [DampedSpring(x, y, k, b, l) for (x, y) in Any[[1,2], [1,3], [1,4], [2,3], [2,4], [3, 4]]]
    positions = Point{3, Float64}[[0.; 0; 0.2], [0.05; 0; 0.2], [0.; 0.05; 0.2], [0.02; 0.02; 0.25]]
    faces = convex_hull(positions)
    robot = SoftRobot(nodes, edges, faces)
    velocities = [Point{3, Float64}(0) for n in robot.nodes]
    state = SoftRobotState(positions, velocities)
    robot, state
end

function blob(;k=4000., b=4.0, m=0.1, num_nodes=20)
    nodes = [PointMass(m) for i in 1:num_nodes]
    positions = Point{3, Float64}[rand(3)+[0.; 0; 0.1] for i in 1:num_nodes]
    edges = DampedSpring[]
    for i = 1:num_nodes
        for j = (i+1):num_nodes
            rest_length = norm(positions[i] - positions[j])
            push!(edges, DampedSpring(i, j, k, b, rest_length))
        end
    end

    faces = convex_hull(positions)
    robot = SoftRobot(nodes, edges, faces)
    velocities = [Point{3, Float64}(0) for n in robot.nodes]
    state = SoftRobotState(positions, velocities)
    robot, state
end

function cube(;k=4000., b=4.0, m=0.1)
    num_nodes = 8
    nodes = [PointMass(m) for i in 1:num_nodes]
    positions = Point{3, Float64}[product([[0.; 1.0] for i = 1:3]...)...]
    edges = DampedSpring[]
    for i = 1:num_nodes
        for j = (i+1):num_nodes
            rest_length = norm(positions[i] - positions[j])
            push!(edges, DampedSpring(i, j, k, b, rest_length))
        end
    end

    faces = convex_hull(positions)
    robot = SoftRobot(nodes, edges, faces)
    velocities = [Point{3, Float64}(0) for n in robot.nodes]
    state = SoftRobotState(positions, velocities)
    robot, state
end


function snake()
    k = 4000.0
    b = 4.0
    l = 0.05
    m = 0.1
    nodes = [PointMass(m)]
    edges = [DampedSpring(x, y, k, b, l) for (x, y) in Any[[1,2], [1,3], [1,4]]]
    x = 0.4
    y = 0.8
    positions = Point{3, Float64}[[x, y, 0]]
    rows = 8
    for row = 1:8
        x += l
        for i = 1:3
            push!(nodes, PointMass(m))
            push!(positions, Point(x, y + (i - 2) * l, 0))
        end
        for i = 1:2
            push!(edges, DampedSpring(1 + i + (row-1) * 3, i + 2 + (row-1) * 3, k, b, l))
        end
        if row != rows
            for i = 1:3
                push!(edges, DampedSpring(1 + i + (row-1) * 3, 1 + i + (row) * 3, k, b, l))
            end
            for i = 1:2
                push!(edges, DampedSpring(i + 2 + (row-1) * 3, 1 + i + row * 3, k, b, sqrt(2)*l))
                push!(edges, DampedSpring(1 + i + (row-1) * 3, i + 2 + row * 3, k, b, sqrt(2)*l))
            end
        end
    end
    x += l
    push!(nodes, PointMass(m))
    push!(positions, Point(x, y, 0))
    for i = length(nodes)-3:length(nodes)-1
        push!(edges, DampedSpring(i, length(nodes), k, b, l))
    end

    faces = convex_hull(positions)
    robot = SoftRobot(nodes, edges, faces)

    velocities = [Point{3, Float64}(0) for n in r.nodes]
    state = SoftRobotState(positions, velocities)
    robot, state
end
