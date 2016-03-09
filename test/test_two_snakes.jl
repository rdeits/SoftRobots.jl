using SoftRobots
using PyPlot
using MultiPoly

function copy_dict_values(dict::Dict)
    Dict(zip(keys(dict), deepcopy(values(dict))))
end

function test_two_snakes()

    world = SoftRobots.World2D()
    world.friction = 1.5
    world_state = SoftRobots.WorldState()

    snake1, snake1_state = SoftRobots.snake()
    push!(world.objects, snake1)
    world_state[snake1] = snake1_state

    snake2, snake2_state = SoftRobots.snake()
    for i in 1:length(snake2_state.positions)
        snake2_state.positions[i] -= [0, 0.15]
    #     snake2_state.velocities[i] += [-0.2, 4]
    end
    push!(world.objects, snake2)
    world_state[snake2] = snake2_state

    x, y = generators(MPoly{Float64}, :x, :y)
    terrain = SoftRobots.FixedObject()
    # terrain_state = SoftRobots.FixedObjectState(-(0.1 + 0.5x^2 - 1y))
    terrain_state = SoftRobots.FixedObjectState(-0.1 + y)
    push!(world.objects, terrain)
    world_state[terrain] = terrain_state

    vis = SoftRobots.PyPlotVisualizer(world)

    history = typeof(world_state)[]

    function playback(history::Array{SoftRobots.WorldState}, dt=1)
        for i = 1:dt:length(history)
            SoftRobots.draw(vis, history[i])
            xlim([0,1])
            ylim([0,1])
            PyPlot.draw()
        end
    end

    pygui(true)
    PyPlot.ion()

    for j = 0:2000
        SoftRobots.update!(world, world_state, 0.001)
        push!(history, copy_dict_values(world_state))
        if mod(j, 100) == 0
            SoftRobots.draw(vis, world_state)
            xlim([0,1])
            ylim([0,1])
        end
    end

    playback(history, 10)
end

test_two_snakes()
