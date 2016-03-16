using SoftRobots
using MultiPoly
# using ProfileView

function test_snake()
    world = SoftRobots.World3D()
    world_state = Dict{SoftRobots.Object, SoftRobots.ObjectState}()

    robot, robot_state = SoftRobots.tetrahedron()
    push!(world.objects, robot)
    world_state[robot] = robot_state

    x, y, z = generators(MPoly{Float64}, :x, :y, :z)
    terrain = SoftRobots.FixedObject()
    terrain_state = SoftRobots.FixedObjectState(-(0.1 + 0.5x^2  + 0.5y^2 - 1z))
    push!(world.objects, terrain)
    world_state[terrain] = terrain_state

    # Profile.clear()
    SoftRobots.update!(world, world_state, 0.001)
    @time @profile for j = 1:1000
        SoftRobots.update!(world, world_state, 0.001)
        SoftRobots.draw(robot, world_state[robot])
    end
    # ProfileView.view()
end

test_snake()
