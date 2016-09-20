module SoftRobots

import DataStructures: OrderedDict
import Iterators: product
import Base: convert
import MultiPoly: MPoly

using SpatialFields
using Meshing
using GeometryTypes
using PyLCM
using PyCall
global spatial
global lcmdrake

function __init__()
    const global spatial = pyimport("scipy.spatial")
    const global lcmdrake = pyimport("drake")
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

type SoftRobot{T, IndexOffset} <: DynamicObject
    nodes::Vector{PointMass}
    edges::Vector{DampedSpring}
    faces::Vector{Face{3, T, IndexOffset}}
end

function SoftRobot{T, IndexOffset}(nodes::Vector{PointMass}, edges::Vector{DampedSpring}, faces::Vector{Face{3, T, IndexOffset}})
  SoftRobot{T, IndexOffset}(nodes, edges, faces)
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
        if isa(object, SoftRobot)
            if any(x -> any(isnan, x), states[object].positions)
                @show states[object]
                error("nans in positions after collision")
            end
        end
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

function update_barrier!{T}(robot::SoftRobot, state::SoftRobotState{T})
    node_is_on_face = falses(length(robot.nodes))
    for face in robot.faces
        for i in 1:length(face)
            node_is_on_face[onebased(face, i)] = true
        end
    end
    mesh = HomogenousMesh(state.positions, robot.faces)
    normals = decompose(GeometryTypes.Normal{3, T}, mesh)
    if any(x -> any(isnan, x), normals[node_is_on_face])
        @show robot.faces
        @show state.positions
        @show state.positions[node_is_on_face]
        @show normals[node_is_on_face]
    end
    # state.barrier = HermiteRadialField(state.positions[node_is_on_face], normals[node_is_on_face])
    values = zeros(length(state.positions))
    values[!node_is_on_face] = -1.0
    state.barrier = InterpolatingSurface(state.positions, values, SpatialFields.XSquaredLogX())
    # state.barrier = InterpolatingSurface(state.positions, values)
end

function update_velocity_field!{T}(robot::SoftRobot, state::SoftRobotState{T})
    # mean_velocity = mean(state.velocities)
    # dimension = length(state.positions[1])
    # state.velocity_field = PolynomialVectorField([MPoly{T}(OrderedDict(zeros(dimension) => v), [:x, :y, :z][1:dimension]) for v in mean_velocity])
    state.velocity_field = PolynomialVectorField(SpatialFields.linear_fit(state.positions, state.velocities))
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
        barrier_velocity::Point{3, T} = evaluate(other_state.velocity_field, my_state.positions[i])
        relative_velocity::Point{3, T} = my_state.velocities[i] - barrier_velocity
        predicted_relative_position::Point{3, T} = my_state.positions[i] + relative_velocity * dt
        if evaluate(other_state.barrier, predicted_relative_position) <= 0
            gradient::Point{3, T} = evaluate(barrier_graident, predicted_relative_position)
            normal = gradient ./ norm(gradient)
            if any(isnan, normal)
                @show gradient
                @show normal
                error("nan in normal")
            end
            penetration_velocity = -dot(relative_velocity, normal)
            if penetration_velocity > 0
                # my_state.velocities[i] += normal * penetration_velocity
                my_state.velocities[i] = barrier_velocity
            end
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

function euler_integrate!{T}(robot::SoftRobot, state::SoftRobotState{T}, total_forces, dt::Real)
    dim = length(state.positions[1])
    for i = 1:length(robot.nodes)
        delta_v = Point{3, T}(total_forces[i]) / robot.nodes[i].mass * dt
        # state.velocities[i] += half_delta_v
        # state.velocities[i] += Point(total_forces[i]) / robot.nodes[i].mass * dt
        state.positions[i] += state.velocities[i] * dt
        state.velocities[i] += delta_v
    end
end

function convex_hull{T}(nodes::Vector{Point{3, T}})
    hull = spatial[:ConvexHull](hcat(map(x -> convert(Vector, x), nodes)...)')
    simplices = hull[:simplices]
    simplices += 1

    faces = Face{3, Int, 0}[]
    # Reorient simplices so that normals always point out from the hull
    for i = 1:size(simplices, 1)
        v = Vector{T}[nodes[simplices[i, j]] for j in 1:size(simplices, 2)]

        # Given a face of a convex hull, all of the points in the body must be on the same side of that face. So in theory, we just need to pick one point not on the face and check which side it's on. But to attempt to be robust to numerical issues, we'll actually sum the dot product with the normal for every point in the body, and check whether that sum is positive or negative. If this becomes a performance bottleneck, we can revisit it later.
        n = Point{3, T}(cross(v[2] - v[1], v[3] - v[1]))
        b = dot(n, nodes[simplices[i, 1]])
        total = zero(T)
        for j = 1:length(nodes)
            total += dot(n, nodes[j]) - b
        end
        if total < 0
            # Then the normal is pointing the right way
            push!(faces, Face{3, Int, 0}((simplices[i,:])...))
        else
            # Otherwise the normal is pointing the wrong way, so we flip the face
            push!(faces, Face{3, Int, 0}((simplices[i,end:-1:1])...))
        end
    end
    return faces
end

function barrier_mesh(bounds::HyperRectangle, barrier::ScalarField)
    ub = maximum(bounds)
    lb = minimum(bounds)
    widths = maximum(bounds) - minimum(bounds)
    barrier_field = SignedDistanceField(x -> SoftRobots.evaluate(barrier, x), bounds, minimum(widths)/20);
    barrier_mesh = HomogenousMesh(barrier_field, 0.0)
    rescaled_points = Point{3,Float64}[Vec(v-1) ./ (Vec(size(barrier_field))-1) .* (ub - lb) + lb for v in vertices(barrier_mesh)]
     HomogenousMesh(rescaled_points, barrier_mesh.faces)
 end

 function barrier_mesh(field::HermiteRadialField)
    center = Vec(mean(field.points))
    widths = Vec(2*(maximum(field.points) - minimum(field.points)))
    lb = Vec(center - 0.5*widths)
    ub = Vec(center + 0.5*widths)
    bounds = HyperRectangle(lb, widths)
    barrier_mesh(bounds, field)
end

function barrier_mesh(field::InterpolatingSurface)
    center = Vec(mean(field.points))
    widths = Vec(2*(maximum(field.points) - minimum(field.points)))
    lb = Vec(center - 0.5*widths)
    ub = Vec(center + 0.5*widths)
    bounds = HyperRectangle(lb, widths)
    barrier_mesh(bounds, field)
end

barrier_mesh(state::SoftRobotState) = barrier_mesh(state.barrier)

include("drawing.jl")

function copy_dict_values(dict::Dict)
    Dict(zip(keys(dict), deepcopy(values(dict))))
end

function simulate{StateType <: WorldState}(world::World, initial_world_state::StateType, times; draw_callback::Function=(args...)->nothing, history_downsampling=1)
    state = copy_dict_values(initial_world_state)
    history = Tuple{Float64, StateType}[]
    dts = diff(times)

    for i in 2:length(times)
        SoftRobots.update!(world, state, times[i] - times[i-1])
        if mod(i, history_downsampling) == 0
            push!(history, (times[i], copy_dict_values(state)))
        end
        draw_callback(world, state)
    end
    history
end

include("robots.jl")


end
