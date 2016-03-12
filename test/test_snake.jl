using SoftRobots
using MultiPoly
using ProfileView

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

    Profile.clear()
    SoftRobots.update!(world, world_state, 0.001)
    @time @profile for j = 1:1000
        SoftRobots.update!(world, world_state, 0.001)
    end
    ProfileView.view()
end

test_snake()
