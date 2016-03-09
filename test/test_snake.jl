using SoftRobots
using MultiPoly
using ProfileView
# using PyPlot

function test_snake()
    world = SoftRobots.World2D()
    world_state = Dict{SoftRobots.Object, SoftRobots.ObjectState}()

    r, snake_state = SoftRobots.snake()
    push!(world.objects, r)
    world_state[r] = snake_state

    x, y = generators(MPoly{Float64}, :x, :y)
    terrain = SoftRobots.FixedObject()
    terrain_state = SoftRobots.FixedObjectState(-(0.1 + 0.5x^2 - 1y))
    push!(world.objects, terrain)
    world_state[terrain] = terrain_state

    # vis = SoftRobots.PyPlotVisualizer(world)

    # pygui(true)
    # PyPlot.ion()

    Profile.clear()
    SoftRobots.update!(world, world_state, 0.001)
    @time @profile for j = 1:1000
        # if mod(j, 100) == 0
        #     SoftRobots.draw(vis, world_state)
        #     xlim([0,1])
        #     ylim([0,1])
        # end
        SoftRobots.update!(world, world_state, 0.001)
    end
    ProfileView.view()
end

test_snake()
