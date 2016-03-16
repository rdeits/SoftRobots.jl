VERSION >= v"0.4" && __precompile__()

module SoftRobots

using SpatialFields
import MultiPoly: MPoly
using GeometryTypes
import DataStructures: OrderedDict
using PyCall
global spatial
using LCMGL

function __init__()
    global spatial = pyimport("scipy.spatial")
end

abstract Object
abstract StaticObject <: Object
abstract DynamicObject <: Object

type DampedSpring
    parent::Int
    child::Int
    stiffness::Float64
    damping::Float64
    rest_length::Float64
end

type PointMass
    mass::Float64
end

type SoftRobot <: DynamicObject
    nodes::Vector{PointMass}
    edges::Vector{DampedSpring}
    faces::Vector{Face{3}}
end

immutable FixedObject <: StaticObject
end

type World
    objects::Array{Object}
    gravity::Vector
    friction::Number

    World(gravity::Vector) = new([], gravity, 1.0)
end

function World3D()
    World([0, 0, -9.81])
end

abstract ObjectState
abstract DynamicObjectState <: ObjectState
abstract StaticObjectState <: ObjectState

typealias WorldState Dict{Object, ObjectState}

type SoftRobotState{T} <: DynamicObjectState
    positions::Vector{Point{3, T}}
    velocities::Vector{Point{3, T}}
    barrier::ScalarField
    velocity_field::VectorField

    SoftRobotState(positions, velocities) = new(positions, velocities)
end

function SoftRobotState{T}(positions::Vector{Vector{T}}, velocities::Vector{Vector{T}})
	SoftRobotState{T}(map(x -> convert(Point{3, T}, x), positions), map(x -> convert(Point{3, T}, x), velocities))
end

function SoftRobotState{T}(positions::Vector{Point{3, T}}, velocities::Vector{Point{3, T}})
    SoftRobotState{T}(map(x -> convert(Point{3, T}, x), positions), map(x -> convert(Point{3, T}, x), velocities))
end

type FixedObjectState <: StaticObjectState
    barrier::ScalarField
    velocity_field::VectorField

    FixedObjectState(barrier) = new(barrier, [MPoly{Float64}(OrderedDict(zeros(Int64, length(barrier.vars)) => 0.0), barrier.vars) for i in 1:length(barrier.vars)])
end

function update!(world::World, states::Dict{Object, ObjectState}, dt::Number)
    for (object, state) in states
        update_barrier!(object, state)
        update_velocity_field!(object, state)
    end
    for object in keys(states)
        collisions!(object, states, world.friction, dt)
    end
    for (object, state) in states
        dynamics!(object, state, world.gravity, dt)
    end
end

function update_barrier!(object::FixedObject, state::FixedObjectState)
    # Nothing to do here
end

function update_velocity_field!(object::FixedObject, args...)
    # Nothing to do here
end

function collisions!(object::FixedObject, args...)
    # Nothing to do here
end

function dynamics!(object::FixedObject, args...)
    # Nothing to do here
end

function normalize!(x)
    n = norm(x)
    @inbounds @simd for i = 1:length(x)
        x[i] /= n
    end
end

function update_barrier!(robot::SoftRobot, state::SoftRobotState)
    # points = hcat([state.positions[s.parent] for s in robot.spines]...)
    # normals = hcat([state.positions[s.parent] - state.positions[s.child] for s in robot.spines]...)
    # for i = 1:size(normals, 1)
    #     normalize!(normals[i,:])
    # end
    # state.barrier = HermiteRadialField(points, normals)
end

function update_velocity_field!{T}(robot::SoftRobot, state::SoftRobotState{T})
    mean_velocity = mean(state.velocities)
    # TODO: linear regression on velocity as a function of [x, y]
    dimension = length(state.positions[1])
    state.velocity_field = PolynomialVectorField([MPoly{T}(OrderedDict(zeros(dimension) => v), [:x, :y, :z][1:dimension]) for v in mean_velocity])
end

function collisions!(robot::SoftRobot, states::Dict{Object, ObjectState}, friction::Number, dt::Number)
    my_state = states[robot]
    for (obj, other_state) in states
        if obj == robot
            # Don't collide with yourself
            continue
        end
        barrier_collisions!(robot, my_state, other_state, friction, dt)
    end
end

function dynamics!(robot::SoftRobot, state::SoftRobotState, gravity::Vector, dt::Number)
    dim = length(state.positions[1])
    total_forces = [robot.nodes[i].mass * gravity for i = 1:length(robot.nodes)]
    update_internal_forces!(robot, state, total_forces)
    euler_integrate!(robot, state, total_forces, dt)
end

function rotmat(theta)
    [cos(theta) -sin(theta);
     sin(theta) cos(theta)]
end

function barrier_collisions!{T}(robot::SoftRobot, my_state::SoftRobotState{T}, other_state::ObjectState, friction::Number, dt::Number)
	barrier_graident = grad(other_state.barrier)
    for i = 1:length(robot.nodes)
        barrier_velocity = Point(evaluate(other_state.velocity_field, my_state.positions[i])...)
        relative_velocity = my_state.velocities[i] - barrier_velocity
        predicted_relative_position = my_state.positions[i] + relative_velocity * dt
        if evaluate(other_state.barrier, predicted_relative_position) <= 0
            # println("collision with barrier: " , other_state.barrier)
            # @show my_state.velocities[i]
            # @show my_state.positions[i]
            # @show relative_velocity
            gradient = Point(evaluate(barrier_graident, predicted_relative_position))
            normal = gradient / norm(gradient)
            penetration_velocity = -dot(relative_velocity, normal)
            # @show penetration_velocity
            if penetration_velocity > 0
                # my_state.velocities[i] += normal * penetration_velocity
                my_state.velocities[i] = barrier_velocity
            end
            # @show my_state.velocities[i]
            # gradient = rotmat(friction * pi*(rand(1) - 0.5)) * gradient
            # normal = gradient / norm(gradient)
            # penetration_velocity = -dot(relative_velocity, normal)
            # if penetration_velocity > 0
            #     my_state.velocities[i] += normal * penetration_velocity
            # end
        end
    end
end

function update_internal_forces!{T}(robot::SoftRobot, state::SoftRobotState{T}, total_forces)
    dim = length(total_forces[1])
    displacement = Array{T}(dim)
    direction = Array{T}(dim)
    relative_velocity = Array{T}(dim)
    force = Array{T}(dim)
    distance = 0.
    edge_velocity = 0.
    force = 0.
    for edge in robot.edges
        for j = 1:dim
            displacement[j] = state.positions[edge.child][j] - state.positions[edge.parent][j]
        end
        distance = norm(displacement)
        for j = 1:dim
            direction[j] = displacement[j] ./ distance
            relative_velocity[j] = state.velocities[edge.child][j] - state.velocities[edge.parent][j]
        end
        edge_velocity = dot(relative_velocity, direction)
        for j = 1:dim
            force = direction[j] * (edge.stiffness * (edge.rest_length - distance) - edge.damping * edge_velocity)
            total_forces[edge.child][j] += force
            total_forces[edge.parent][j] -= force
        end
    end
end

function euler_integrate!(robot::SoftRobot, state::SoftRobotState, total_forces, dt::Real)
    dim = length(state.positions[1])
    for i = 1:length(robot.nodes)
        half_delta_v = Point(total_forces[i]) / robot.nodes[i].mass * dt / 2.0
        state.velocities[i] += half_delta_v
        # state.velocities[i] += Point(total_forces[i]) / robot.nodes[i].mass * dt
        state.positions[i] += state.velocities[i] * dt
        state.velocities[i] += half_delta_v
    end
end

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

function convex_hull(nodes)
    hull = spatial[:ConvexHull](hcat(map(x -> convert(Vector, x), nodes)...)')
    simplices = hull[:simplices]
    return [Face((simplices[i,:]+1)...) for i in 1:size(simplices, 1)]
end

function tetrahedron()
    # k = 4000.0
    k = 4000.0
    b = 4.0
    l = 0.05
    m = 0.1

    nodes = [PointMass(m) for i in 1:4]
    edges = [DampedSpring(x, y, k, b, l) for (x, y) in Any[[1,2], [1,3], [1,4], [2,3], [2,4], [3, 4]]]
    positions = [Point(0., 0, 0.2), Point(0.05, 0, 0.2), Point(0., 0.05, 0.2), Point(0.02, 0.02, 0.25)]
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
        # for k = 1:edges_per_node
        #     connected_nodes = []
        #     j = rand(1:num_nodes)
        #     while j == i || j in connected_nodes
        #         j = rand(1:num_nodes)
        #     end
        #     push!(connected_nodes, j)
        #     rest_length = norm(positions[i] - positions[j])
        #     push!(edges, DampedSpring(i, j, k, b, rest_length))
        # end
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

function mesh{T}(robot::SoftRobot, state::SoftRobotState{T})
    HomogenousMesh(state.positions, robot.faces)
end

function draw(mesh::AbstractMesh)
    LCMGLClient("soft_robot") do lcmgl
        color(lcmgl, .7, .7, .2, 0.4)
        begin_mode(lcmgl, LCMGL.TRIANGLES)
        verts = vertices(mesh)
        for face in faces(mesh)
            for point in verts[face][end:-1:1]
                vertex(lcmgl, point[1], point[2], point[3])
            end
        end
        end_mode(lcmgl)
        color(lcmgl, 0., 0., 0., 1.0)
        begin_mode(lcmgl, LCMGL.LINES)
        verts = vertices(mesh)
        for face in faces(mesh)
            for i = 1:length(verts[face])
                vertex(lcmgl, verts[face][i]...)
                j = i + 1
                if i == length(verts[face])
                    j = 1
                end
                vertex(lcmgl, verts[face][j]...)
            end
        end
        end_mode(lcmgl)
        switch_buffer(lcmgl)
    end
end

draw(robot::SoftRobot, state::SoftRobotState) = draw(mesh(robot, state))

# include("pyplot_visualizer.jl")

# type Range
#     min
#     max
# end

# function draw(barrier::ScalarField, variable_ranges::Array{Range, 1}, level=0.0)
#     @assert length(variable_ranges) == 2
#     X = linspace(variable_ranges[1].min, variable_ranges[1].max)
#     Y = linspace(variable_ranges[2].min, variable_ranges[2].max)

#     Z = zeros(length(Y), length(X))
#     for j = 1:length(X)
#         for k = 1:length(Y)
#             Z[k,j] = evaluate(barrier, [X[j], Y[k]])
#         end
#     end
#     PyPlot.contour(X, Y, Z, [level, level])
# end

# function draw(snake::SoftRobot, state::SoftRobotState)
#     draw(state.barrier, [Range(0,1), Range(0,1)])

#     for edge in snake.edges
#         p1 = state.positions[edge.parent]
#         p2 = state.positions[edge.child]
#         PyPlot.plot([p1[1], p2[1]], [p1[2], p2[2]], "bo-")
#     end
# end

# function draw(object::FixedObject, state::FixedObjectState)
#     draw(state.barrier, [Range(0,1), Range(0,1)])
# end

# function draw(states::Dict{Object, ObjectState})
#     for (object, state) in states
#         draw(object, state)
#     end
# end

end
