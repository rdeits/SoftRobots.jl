using ForwardDiff
using PyPlot
using NLopt
using Ipopt
using JuMP
import ForwardDiff: value
using FixedSizeArrays
import SoftRobots
reload("SoftRobots")

function potential_energy{T}(robot, state::SoftRobots.SoftRobotState{T})
    energy = zero(T)
    for edge in robot.edges
        energy += 0.5 * edge.stiffness * 
         (norm(state.positions[edge.parent] - state.positions[edge.child]) - 
          edge.rest_length)^2
    end
    energy
end

Base.convert{T}(::Type{SoftRobots.SoftRobotState{T}}, x::Vector{T}) = SoftRobots.SoftRobotState(SoftRobots.Point{3, T}[x[i:i+2]
    for i in 1:3:(3*length(robot.nodes)-2)],
    [SoftRobots.Point{3, T}(0) for i in 1:length(robot.nodes)])

potential_energy{T}(robot, x::Vector{T}) = potential_energy(robot, convert(SoftRobots.SoftRobotState{T}, x))


draw{T}(robot, x::Vector{T}) = SoftRobots.draw(robot, convert(SoftRobots.SoftRobotState{T}, x))

value(x::Real) = x

function tower()
    k = 10.
    b = 4.0
    l = 1.0
    m = 0.1
    nodes = SoftRobots.PointMass[]
    positions = SoftRobots.Point{3, Float64}[]
    for z = 0:1:3
        bulge = 1 - 0.1*(z - 1.5)^2
        for x = -0.5:0.5
            for y = -0.5:0.5
                push!(nodes, SoftRobots.PointMass(m))
                push!(positions, SoftRobots.Point{3, Float64}(x*bulge, y*bulge, z))
            end
        end
    end
    num_nodes = length(nodes)
    edges = SoftRobots.DampedSpring[]
    for i = 1:num_nodes
        for j = (i+1):num_nodes
            rest_length = norm(positions[i] - positions[j])
            push!(edges, SoftRobots.DampedSpring(i, j, k, b, rest_length))
        end
    end
    faces = SoftRobots.convex_hull(positions)
    robot = SoftRobots.SoftRobot(nodes, edges, faces)
    velocities = [SoftRobots.Point{3, Float64}(0) for n in robot.nodes]
    state = SoftRobots.SoftRobotState(positions, velocities)
    robot, state
end
    
        
robot, state = tower()
SoftRobots.draw(robot, state)
point = [0.25; 0.0]

function potential_energy_fixed_base{T}(robot, state::SoftRobots.SoftRobotState{T})
    energy = zero(T)
    for edge in robot.edges
        energy += 0.5 * edge.stiffness * 
         (norm(state.positions[edge.parent] - state.positions[edge.child]) - 
          edge.rest_length)^2
    end
    for position in state.positions[1:4]
        energy += 0.5 * 1000 * position[3]^2
    end
    energy
end

function cost{T}(x::Vector{T})
    state = convert(SoftRobots.SoftRobotState{T}, x)
    energy = potential_energy_fixed_base(robot, state)
    SoftRobots.update_barrier!(robot, state)
#     v = SoftRobots.evaluate(state.barrier, point)
#     if v > 0
#         energy += 0.5 * 1000 * v^2
#     end
    energy
end

x0 = collect(destructure(state.positions))
f = cost
g = ForwardDiff.gradient(f)
x = copy(x0)
energies = []
step = 1e-4
@time for i = 1:1
    gi = g(x)
    x -= step * gi
    draw(robot, map(value, x))
    push!(energies, f(x))
end
