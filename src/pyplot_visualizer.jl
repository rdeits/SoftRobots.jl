abstract Visualizer

type PyPlotVisualizer <: Visualizer
	world::World
	fig
	ax
	actors
	drawn_static_objects
end

type Actors
	edges
	barrier
end

function PyPlotVisualizer(world::World)
	fig = PyPlot.figure()
	ax = fig[:add_subplot](1,1,1)
	actors = Dict([(obj, construct_actors(ax, obj)) for obj in world.objects])
	drawn_static_objects = Set{Object}()
	PyPlotVisualizer(world, fig, ax, actors, drawn_static_objects)
end

construct_actors(ax, object::StaticObject) = Actors([], ax[:contour](zeros(25, 25), [0.0, 0.0]))

function construct_actors(ax, object::SoftRobot)
	edge_actors = []
	for edge in object.edges
		push!(edge_actors, ax[:plot]([0,0], [0,0], "bo-"))
	end
	Actors(edge_actors, ax[:contour](zeros(25, 25), [0.0, 0.0]))
end

function draw(vis::PyPlotVisualizer, state::WorldState)
	for object in vis.world.objects
		draw(vis, object, state[object])
	end
	PyPlot.draw()
end

function draw(vis::PyPlotVisualizer, object::StaticObject, state::StaticObjectState)
	if !(object in vis.drawn_static_objects)
		draw(vis.actors[object], object, state)
	end
end

function draw(vis::PyPlotVisualizer, object::SoftRobot, state::SoftRobotState)
	draw(vis.actors[object], object, state)
end

type Range
    min
    max
end

function draw(actors::Actors, object::StaticObject, state::StaticObjectState)
	draw(actors.barrier, state.barrier, [Range(0, 1), Range(0, 1)])
end

function draw(actor, barrier::ScalarField, variable_ranges::Array{Range, 1}, level=0.0)
    @assert length(variable_ranges) == 2
    X = linspace(variable_ranges[1].min, variable_ranges[1].max, 25)
    Y = linspace(variable_ranges[2].min, variable_ranges[2].max, 25)
    
    Z = zeros(length(Y), length(X))
    for j = 1:length(X)
        for k = 1:length(Y)
            Z[k,j] = evaluate(barrier, [X[j], Y[k]])
        end
    end
    PyCall.pyeval("[actor.ax.collections.remove(c) for c in actor.collections]", actor=actor)
    actor[:__init__](actor[:ax], X, Y, Z, [level, level])
end

function draw(actors::Actors, robot::SoftRobot, state::SoftRobotState)
	draw(actors.barrier, state.barrier, [Range(0, 1), Range(0, 1)])

	for (i, edge) in enumerate(robot.edges)
		actors.edges[i][1][:set_xdata]([state.positions[edge.parent][1], state.positions[edge.child][1]])
		actors.edges[i][1][:set_ydata]([state.positions[edge.parent][2], state.positions[edge.child][2]])
	end
end

