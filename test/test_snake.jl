using SoftRobots
using MultiPoly
using PyPlot

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

	pygui(true)
	PyPlot.ion()

	for j = 1:2000
	    if mod(j, 5) == 0
	    	cla()
	        SoftRobots.draw(world_state)
	        xlim([0,1])
	        ylim([0,1])
	        PyPlot.draw()
	    end
	    SoftRobots.update!(world, world_state, 0.001)
	end
end

test_snake()
